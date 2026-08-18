#!/usr/bin/env python3
"""Decide which paired animals really are one physical fish, from LD-vs-HD concordance.

A cross-array test pairs an animal's LD and HD genotypes on its ID and assumes
both typings came from the same fish. When an ID covers two physical samples --
a mislabelled tube, a reused ID -- the pair is two unrelated fish, and the run
scores an imputation that never had a chance.

The two typings share the markers both arrays carry, so the check needs no
imputation at all: genotype concordance between LD and HD at those markers,
per animal. The result is sharply bimodal -- the same animal lands at 0.95-1.0,
two different salmon land near 0.5, and nothing sits in between -- so any
threshold in the gap works.

On the 2026-08 salmon benchmark this separated 974 of 3,881 animals (25%) in
one step, which had scored mean R2 0.53 purely because a quarter of its truth
was a different fish.

Reads the per-individual concordance written by compute_accuracy_metrics.py and
writes the animal list every downstream rule is built from.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def read_pairs(path: Path) -> pd.DataFrame:
    pairs = pd.read_csv(path, sep=r"\s+", header=None, names=["fid", "iid"], dtype=str)
    if pairs.empty:
        raise ValueError(f"No paired animals found in {path}")
    return pairs


def classify(pairs: pd.DataFrame, metrics: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Attach each animal's LD-vs-HD concordance and its pass/fail verdict.

    An animal the metrics never scored -- no overlapping non-missing calls --
    cannot be shown to be the same fish, so it fails rather than defaulting to
    pass. Silently keeping it is how an unverified pair reaches the results.
    """
    scored = metrics.rename(columns={"sample": "iid"})[["iid", "concordance", "n_evaluated"]]
    merged = pairs.merge(scored, on="iid", how="left")
    merged["concordance"] = pd.to_numeric(merged["concordance"], errors="coerce")
    merged["passed"] = merged["concordance"].ge(threshold).fillna(False)
    return merged


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", required=True, type=Path,
                        help="FID/IID list of animals present on both arrays")
    parser.add_argument("--metrics", required=True, type=Path,
                        help="metrics_by_individual.tsv from the LD-vs-HD comparison")
    parser.add_argument("--threshold", type=float, default=0.90,
                        help="Minimum LD-vs-HD concordance for one physical fish")
    parser.add_argument("--min-pass-rate", type=float, default=0.5,
                        help="Abort below this pass rate; that is a data problem, not a result")
    parser.add_argument("--identity-out", required=True, type=Path)
    parser.add_argument("--fail-out", required=True, type=Path)
    parser.add_argument("--keep-out", required=True, type=Path)
    args = parser.parse_args()

    pairs = read_pairs(args.pairs)
    metrics = pd.read_csv(args.metrics, sep="\t")
    verdicts = classify(pairs, metrics, args.threshold)

    for path in (args.identity_out, args.fail_out, args.keep_out):
        path.parent.mkdir(parents=True, exist_ok=True)

    verdicts.to_csv(args.identity_out, sep="\t", index=False)
    failed = verdicts.loc[~verdicts["passed"], ["fid", "iid"]]
    kept = verdicts.loc[verdicts["passed"], ["fid", "iid"]]
    failed.to_csv(args.fail_out, sep="\t", header=False, index=False)
    kept.to_csv(args.keep_out, sep="\t", header=False, index=False)

    n_total = len(verdicts)
    pass_rate = len(kept) / n_total if n_total else 0.0
    print(f"paired_animals={n_total}")
    print(f"passed={len(kept)}")
    print(f"failed={len(failed)}")
    print(f"pass_rate={pass_rate:.4f}")
    print(f"threshold={args.threshold}")

    if kept.empty:
        print("ERROR: no animal passed the identity check.", file=sys.stderr)
        sys.exit(1)

    if pass_rate < args.min_pass_rate:
        print(
            f"ERROR: only {pass_rate:.1%} of paired animals are the same fish on both "
            f"arrays, below --min-pass-rate {args.min_pass_rate:.1%}. Fix the pairing "
            "upstream; an accuracy number from this cohort would not mean anything.",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
