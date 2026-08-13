#!/usr/bin/env python3
"""Aggregate per-fold imputation accuracy summaries."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _to_long(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalise a per-fold summary.tsv to long (metric, value) form.

    compute_accuracy_metrics.build_summary writes one wide row with a column
    per metric, so that is the expected shape. A file that is already long is
    passed through, which keeps hand-written summaries usable.
    """
    if set(df.columns) >= {"metric", "value"}:
        return df[["metric", "value"]].copy()
    return df.melt(var_name="metric", value_name="value")


def aggregate_summaries(
    root: str | Path,
    imputers: list[str],
    folds: list[int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return long per-fold metrics and imputer-level mean/sd metrics."""
    root = Path(root)
    frames = []
    for imputer in imputers:
        for fold in folds:
            path = root / imputer / f"fold{fold}" / "summary.tsv"
            if not path.exists():
                continue
            df = _to_long(pd.read_csv(path, sep="\t"))
            df.insert(0, "fold", fold)
            df.insert(0, "imputer", imputer)
            frames.append(df)

    if not frames:
        raise ValueError(f"No per-fold summary.tsv files found under {root}")

    long_df = pd.concat(frames, ignore_index=True)
    long_df["value"] = pd.to_numeric(long_df["value"], errors="coerce")
    imputer_df = (
        long_df.groupby(["imputer", "metric"], as_index=False)["value"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    return long_df, imputer_df


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--imputers", nargs="+", required=True)
    parser.add_argument("--folds", nargs="+", required=True, type=int)
    parser.add_argument("--summary-out", required=True, type=Path)
    parser.add_argument("--imputer-summary-out", required=True, type=Path)
    args = parser.parse_args()

    long_df, imputer_df = aggregate_summaries(args.root, args.imputers, args.folds)
    args.summary_out.parent.mkdir(parents=True, exist_ok=True)
    args.imputer_summary_out.parent.mkdir(parents=True, exist_ok=True)
    long_df.to_csv(args.summary_out, sep="\t", index=False)
    imputer_df.to_csv(args.imputer_summary_out, sep="\t", index=False)


if __name__ == "__main__":
    main()
