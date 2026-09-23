#!/usr/bin/env python3
"""Intersect two or more PLINK .bim files and write the shared marker IDs.

Used for two things that must agree exactly:

  * the LD panel of the cross_array mode -- the markers both arrays carry, and
    therefore the only ones the imputer may be handed as observed genotypes;
  * the `cv.target_snp_list` of a masked run, e.g. the 50K markers a V4 fish is
    downsampled to.

Defining it once means the masked test and the real cross-array test mask to
the same set, which is the only way their per-marker numbers can be compared.

A marker ID present on both arrays but sitting at a different physical position,
or carrying a different allele pair, is not the same marker. Those are dropped
and counted rather than silently trusted.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

BIM_COLUMNS = ["chrom", "snp", "cm", "pos", "a1", "a2"]


def read_bim(path: Path) -> pd.DataFrame:
    bim = pd.read_csv(path, sep=r"\s+", header=None, names=BIM_COLUMNS, dtype=str)
    if bim.empty:
        raise ValueError(f"No variants found in {path}")
    duplicated = bim["snp"].duplicated()
    if duplicated.any():
        examples = ", ".join(bim.loc[duplicated, "snp"].head(5))
        raise ValueError(f"{duplicated.sum()} duplicate marker IDs in {path}: {examples}")
    return bim


def _allele_pair(frame: pd.DataFrame) -> pd.Series:
    """Order-independent allele pair, so an A1/A2 swap is not a mismatch."""
    lo = frame[["a1", "a2"]].min(axis=1)
    hi = frame[["a1", "a2"]].max(axis=1)
    return lo + "/" + hi


def shared_markers(bims: list[pd.DataFrame]) -> tuple[list[str], dict[str, int]]:
    """Marker IDs carried by every panel at the same position and allele pair."""
    merged = None
    for index, bim in enumerate(bims):
        frame = bim[["snp", "chrom", "pos"]].copy()
        frame[f"alleles{index}"] = _allele_pair(bim)
        frame = frame.rename(columns={"chrom": f"chrom{index}", "pos": f"pos{index}"})
        merged = frame if merged is None else merged.merge(frame, on="snp", how="inner")

    stats = {"shared_ids": len(merged)}

    position_ok = pd.Series(True, index=merged.index)
    allele_ok = pd.Series(True, index=merged.index)
    for index in range(1, len(bims)):
        position_ok &= (merged["chrom0"] == merged[f"chrom{index}"]) & (
            merged["pos0"] == merged[f"pos{index}"]
        )
        allele_ok &= merged["alleles0"] == merged[f"alleles{index}"]

    stats["position_mismatch"] = int((~position_ok).sum())
    stats["allele_mismatch"] = int((position_ok & ~allele_ok).sum())

    kept = merged.loc[position_ok & allele_ok].copy()
    kept["_chrom_order"] = kept["chrom0"].map(_chrom_key)
    kept["_pos_int"] = kept["pos0"].astype(int)
    kept = kept.sort_values(["_chrom_order", "_pos_int"])
    stats["kept"] = len(kept)
    return kept["snp"].tolist(), stats


def _chrom_key(value: str):
    return (0, int(value)) if str(value).isdigit() else (1, 0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bim", required=True, nargs="+", type=Path,
                        help="Two or more PLINK .bim files to intersect")
    parser.add_argument("--out", required=True, type=Path,
                        help="Destination for the shared marker ID list")
    parser.add_argument("--min-shared", type=int, default=1,
                        help="Abort if fewer than this many markers survive")
    args = parser.parse_args()

    if len(args.bim) < 2:
        parser.error("--bim needs at least two .bim files")

    bims = [read_bim(path) for path in args.bim]
    markers, stats = shared_markers(bims)

    for path, bim in zip(args.bim, bims):
        print(f"{path}: {len(bim)} markers")
    for key, value in stats.items():
        print(f"{key}={value}")

    if len(markers) < args.min_shared:
        print(
            f"ERROR: only {len(markers)} shared markers, below --min-shared "
            f"{args.min_shared}. The two panels do not describe the same assay.",
            file=sys.stderr,
        )
        sys.exit(1)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(markers) + "\n")


if __name__ == "__main__":
    main()
