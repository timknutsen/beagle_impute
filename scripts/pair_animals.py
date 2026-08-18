#!/usr/bin/env python3
"""Pair animals across two PLINK filesets on individual ID.

Per-array PLINK exports normally carry FID=0 for everyone, and the two arrays
are exported in separate runs, so matching on the (FID, IID) pair is a trap: a
difference in FID alone yields an empty intersection and the run fails with no
hint as to why. The individual ID is the only column that identifies the fish,
so that is what is matched, and the FID is taken from the HD side so every
downstream --keep file agrees with the truth fileset.

Pairing on ID says only that the two rows claim to be the same animal. Whether
they are the same physical fish is what the identity gate decides next; see
scripts/identity_gate.py.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

FAM_COLUMNS = ["fid", "iid", "sire", "dam", "sex", "pheno"]


def read_fam(path: Path) -> pd.DataFrame:
    fam = pd.read_csv(path, sep=r"\s+", header=None, names=FAM_COLUMNS, dtype=str)
    if fam.empty:
        raise ValueError(f"No samples found in {path}")
    duplicated = fam["iid"].duplicated()
    if duplicated.any():
        examples = ", ".join(fam.loc[duplicated, "iid"].head(5))
        raise ValueError(
            f"{duplicated.sum()} duplicate individual IDs in {path}: {examples}. "
            "One row per animal per array is required; deduplicate upstream."
        )
    return fam


def pair(ld_fam: pd.DataFrame, hd_fam: pd.DataFrame) -> pd.DataFrame:
    shared = hd_fam.loc[hd_fam["iid"].isin(set(ld_fam["iid"])), ["fid", "iid"]]
    return shared.sort_values("iid").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ld-fam", required=True, type=Path)
    parser.add_argument("--hd-fam", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    ld_fam = read_fam(args.ld_fam)
    hd_fam = read_fam(args.hd_fam)
    paired = pair(ld_fam, hd_fam)

    print(f"ld_animals={len(ld_fam)}")
    print(f"hd_animals={len(hd_fam)}")
    print(f"paired={len(paired)}")

    if paired.empty:
        print(
            "ERROR: no individual ID appears on both arrays. Check that the two "
            "filesets were exported with the same ID convention.",
            file=sys.stderr,
        )
        sys.exit(1)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    paired.to_csv(args.out, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
