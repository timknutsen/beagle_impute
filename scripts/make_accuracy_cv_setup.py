#!/usr/bin/env python3
"""Create deterministic fold assignments and a shared LD SNP panel."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def assign_folds(samples: pd.DataFrame, n_folds: int, seed: int) -> pd.DataFrame:
    """Return samples with a balanced 1-based fold assignment."""
    if n_folds < 2:
        raise ValueError("n_folds must be at least 2")
    if samples.empty:
        raise ValueError("No samples available for CV fold assignment")

    shuffled = samples[["fid", "iid"]].sample(frac=1.0, random_state=seed).reset_index(drop=True)
    shuffled["fold"] = [(idx % n_folds) + 1 for idx in range(len(shuffled))]
    return shuffled.sort_values(["fold", "fid", "iid"]).reset_index(drop=True)


def choose_ld_panel(
    bim: pd.DataFrame,
    n_snps: int,
    seed: int,
    snp_list: str | Path | None = None,
) -> list[str]:
    """Return a deterministic LD panel, or copy the order from an explicit SNP list."""
    if snp_list:
        path = Path(snp_list)
        snps = [line.strip() for line in path.read_text().splitlines() if line.strip()]
        if not snps:
            raise ValueError(f"Explicit SNP list is empty: {path}")
        return snps

    if n_snps < 1:
        raise ValueError("n_snps must be at least 1")
    n = min(n_snps, len(bim))
    return bim[1].sample(n=n, random_state=seed).sort_values().astype(str).tolist()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fam", required=True, type=Path)
    parser.add_argument("--bim", required=True, type=Path)
    parser.add_argument("--n-folds", required=True, type=int)
    parser.add_argument("--target-n-snps", required=True, type=int)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--target-snp-list", default="", type=Path)
    parser.add_argument("--folds-out", required=True, type=Path)
    parser.add_argument("--snps-out", required=True, type=Path)
    args = parser.parse_args()

    fam = pd.read_csv(
        args.fam,
        sep=r"\s+",
        header=None,
        names=["fid", "iid", "pat", "mat", "sex", "phen"],
    )
    bim = pd.read_csv(args.bim, sep=r"\s+", header=None)

    folds = assign_folds(fam[["fid", "iid"]], args.n_folds, args.seed)
    snps = choose_ld_panel(
        bim,
        n_snps=args.target_n_snps,
        seed=args.seed,
        snp_list=args.target_snp_list if str(args.target_snp_list) else None,
    )

    args.folds_out.parent.mkdir(parents=True, exist_ok=True)
    args.snps_out.parent.mkdir(parents=True, exist_ok=True)
    folds.to_csv(args.folds_out, sep="\t", index=False)
    args.snps_out.write_text("\n".join(snps) + "\n")


if __name__ == "__main__":
    main()
