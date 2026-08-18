#!/usr/bin/env python3
"""Aggregate per-marker accuracy across folds and across runs into one reliability table.

compute_accuracy_metrics.py writes metrics_by_snp.tsv for every imputer and
every fold, and nothing has ever read it -- aggregate_cv_metrics.py consumes
only the run-level summary. That per-marker file is what answers the question
this benchmark exists for: which markers impute consistently.

Two rules make the answer usable rather than flattering:

  * A marker is filtered on its *minimum* R2 across folds, not its mean. A
    marker that collapses in one fold out of five is not a marker you want in a
    reference panel; averaging hides exactly that.
  * Passing --root more than once intersects runs. A V4 marker earns its place
    only by holding up in both the real V3 -> V4 array transition and the
    masked 50K -> V4 test, which stress different things.

Marker IDs are recovered from a .bim by chromosome and position, because
metrics_by_snp.tsv keys on CHROM:POS:REF:ALT and a --extract list needs IDs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

BIM_COLUMNS = ["chrom", "snp", "cm", "pos", "a1", "a2"]


def load_fold_metrics(root: Path, imputers: list[str], folds: list[int]) -> pd.DataFrame:
    """Stack every metrics_by_snp.tsv under one run directory."""
    frames = []
    for imputer in imputers:
        for fold in folds:
            path = root / imputer / f"fold{fold}" / "metrics_by_snp.tsv"
            if not path.exists():
                continue
            frame = pd.read_csv(path, sep="\t")
            frame["imputer"] = imputer
            frame["fold"] = fold
            frame["run"] = root.name
            frames.append(frame)
    if not frames:
        raise ValueError(
            f"No metrics_by_snp.tsv found under {root} for imputers "
            f"{', '.join(imputers)} and folds {', '.join(str(f) for f in folds)}"
        )
    return pd.concat(frames, ignore_index=True)


def summarise_markers(metrics: pd.DataFrame) -> pd.DataFrame:
    """One row per marker per imputer, carrying the worst fold as well as the mean."""
    grouped = metrics.groupby(["run", "imputer", "variant"], as_index=False).agg(
        n_folds_scored=("allelic_r2", "count"),
        mean_r2=("allelic_r2", "mean"),
        min_r2=("allelic_r2", "min"),
        sd_r2=("allelic_r2", "std"),
        mean_concordance=("concordance", "mean"),
        min_concordance=("concordance", "min"),
        maf=("maf", "mean"),
        n_evaluated=("n_evaluated", "sum"),
    )
    return grouped


def attach_marker_ids(markers: pd.DataFrame, bim_path: Path) -> pd.DataFrame:
    """Map CHROM:POS:REF:ALT back to the marker ID a --extract list needs."""
    bim = pd.read_csv(bim_path, sep=r"\s+", header=None, names=BIM_COLUMNS, dtype=str)
    lookup = bim[["chrom", "pos", "snp"]].drop_duplicates(subset=["chrom", "pos"])

    parts = markers["variant"].str.split(":", n=3, expand=True)
    markers = markers.copy()
    markers["chrom"] = parts[0]
    markers["pos"] = parts[1]
    merged = markers.merge(lookup, on=["chrom", "pos"], how="left")
    unmapped = merged["snp"].isna().sum()
    if unmapped:
        print(f"WARNING: {unmapped} scored markers had no ID in {bim_path}", file=sys.stderr)
    return merged


def combine_runs(per_run: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """Collapse runs to one row per marker, keeping the worst run's numbers."""
    return per_run.groupby(keys, as_index=False).agg(
        n_runs=("run", "nunique"),
        n_folds_scored=("n_folds_scored", "sum"),
        mean_r2=("mean_r2", "mean"),
        min_r2=("min_r2", "min"),
        mean_concordance=("mean_concordance", "mean"),
        min_concordance=("min_concordance", "min"),
        maf=("maf", "mean"),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, nargs="+", type=Path,
                        help="One or more accuracy output directories; several intersect")
    parser.add_argument("--imputers", required=True, nargs="+")
    parser.add_argument("--folds", required=True, nargs="+", type=int)
    parser.add_argument("--bim", type=Path,
                        help="PLINK .bim used to recover marker IDs from positions")
    parser.add_argument("--r2-threshold", type=float, default=0.90,
                        help="A marker is reliable when its worst fold clears this")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--reliable-out", type=Path,
                        help="Marker ID list of everything that cleared the threshold")
    args = parser.parse_args()

    per_run = pd.concat(
        [
            summarise_markers(load_fold_metrics(root, args.imputers, args.folds))
            for root in args.root
        ],
        ignore_index=True,
    )

    keys = ["variant"]
    if args.bim:
        per_run = attach_marker_ids(per_run, args.bim)
        keys = ["variant", "snp"]

    combined = combine_runs(per_run, keys)
    # A marker missing from one run was never stress-tested there, so it is not
    # eligible however well it did in the runs that did score it.
    n_runs_expected = len(args.root)
    combined["reliable"] = (
        (combined["n_runs"] == n_runs_expected)
        & (combined["min_r2"] >= args.r2_threshold)
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    combined.sort_values("min_r2", ascending=False).to_csv(args.out, sep="\t", index=False)

    reliable = combined.loc[combined["reliable"]]
    print(f"runs={n_runs_expected}")
    print(f"markers_scored={len(combined)}")
    print(f"markers_reliable={len(reliable)}")
    print(f"r2_threshold={args.r2_threshold}")

    if args.reliable_out:
        args.reliable_out.parent.mkdir(parents=True, exist_ok=True)
        if "snp" in reliable.columns:
            ids = reliable["snp"].dropna().tolist()
        else:
            ids = reliable["variant"].tolist()
        args.reliable_out.write_text("\n".join(ids) + ("\n" if ids else ""))


if __name__ == "__main__":
    main()
