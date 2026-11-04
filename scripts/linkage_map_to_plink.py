#!/usr/bin/env python3
"""Write a Beagle `map=` file (PLINK .map: chrom, id, cM, bp) from a linkage map table.

Built for the Brekke et al. 2023 Atlantic salmon map
(viktige_filer/linkage_maps/AtlanticSalmon_LinkageMap_Brekke_2021_affyID.txt;
columns chr snpid Mb male_cM female_cM avg_cM, where `Mb` holds base pairs).
Beagle phases one population, so the sex-averaged column is the default.

Beagle interpolates between map positions and needs, per chromosome, positions
in increasing order and cM that never decreases. Rows sharing a position are
collapsed to the first, and any cM dip is raised to the running maximum; both
are counted rather than silently fixed.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def build_map(table: pd.DataFrame, chrom: str, snp: str, bp: str, cm: str) -> tuple[pd.DataFrame, dict]:
    d = table[[chrom, snp, cm, bp]].dropna().copy()
    d.columns = ["chrom", "snp", "cm", "bp"]
    d["chrom"] = d["chrom"].astype(str)
    d["bp"] = d["bp"].astype("int64")
    key = d["chrom"].map(lambda c: (0, int(c), "") if c.isdigit() else (1, 0, c))
    d = d.assign(_k=key).sort_values(["_k", "bp"], kind="stable").drop(columns="_k")
    n_dup = int(d.duplicated(["chrom", "bp"]).sum())
    d = d.drop_duplicates(["chrom", "bp"])
    raised = d.groupby("chrom", sort=False)["cm"].cummax()
    n_dip = int((raised > d["cm"]).sum())
    d["cm"] = raised
    return d, {"markers": len(d), "duplicate_positions_dropped": n_dup, "cm_dips_raised": n_dip}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", type=Path, help="whitespace-separated linkage map with a header row")
    ap.add_argument("out", type=Path, help="PLINK .map to write")
    ap.add_argument("--chrom-col", default="chr")
    ap.add_argument("--id-col", default="snpid")
    ap.add_argument("--bp-col", default="Mb", help="physical position in bp (Brekke's file calls it Mb)")
    ap.add_argument("--cm-col", default="avg_cM", help="avg_cM (sex-averaged), male_cM or female_cM")
    a = ap.parse_args(argv)

    table = pd.read_csv(a.table, sep=r"\s+")
    missing = [c for c in (a.chrom_col, a.id_col, a.bp_col, a.cm_col) if c not in table.columns]
    if missing:
        sys.exit(f"{a.table}: missing column(s) {missing}; found {list(table.columns)}")
    d, stats = build_map(table, a.chrom_col, a.id_col, a.bp_col, a.cm_col)
    d.to_csv(a.out, sep="\t", header=False, index=False, float_format="%.6f")
    print(f"{a.out}: {stats['markers']} markers on {d['chrom'].nunique()} chromosomes; "
          f"{stats['duplicate_positions_dropped']} duplicate positions dropped; "
          f"{stats['cm_dips_raised']} cM dips raised", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
