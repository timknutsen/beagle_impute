#!/usr/bin/env python3
"""Choose the members of a reference panel, and record why each animal is in or out.

A breeding cohort is not a diverse sample. Full sibs share long haplotype
blocks, so the tenth sib in a family adds almost nothing a reference panel can
use while costing the same phasing time as an unrelated fish -- and it drags the
panel's allele frequencies toward whichever families happened to be bred most.
Capping per sire x dam pair buys diversity per animal phased.

Exclusions are applied here, before phasing, not by subsetting a phased panel
afterwards. Removing samples from an already-phased VCF is cheaper but leaks:
the excluded animals' genotypes informed the phase estimates of the relatives
that remain, so an accuracy benchmark built that way scores against a panel its
own test animals helped build.

Writes the ID list plus a manifest naming, for every animal, the reason it was
kept or dropped -- so a panel can be explained a year later.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

FAM_COLUMNS = ["fid", "iid", "sire", "dam", "sex", "phen"]


def read_fam(path: Path) -> pd.DataFrame:
    fam = pd.read_csv(path, sep=r"\s+", header=None, names=FAM_COLUMNS, dtype=str)
    if fam.empty:
        raise ValueError(f"No samples found in {path}")
    if fam["iid"].duplicated().any():
        dupes = fam.loc[fam["iid"].duplicated(), "iid"].head(5).tolist()
        raise ValueError(f"Duplicate individual IDs in {path}: {', '.join(dupes)}")
    return fam


def read_pedigree(path: Path) -> pd.DataFrame:
    """Three or more columns; the first three are read as iid, sire, dam."""
    ped = pd.read_csv(path, sep=r"\s+", header=None, dtype=str, comment="#")
    if ped.shape[1] < 3:
        raise ValueError(f"{path} needs at least three columns: iid sire dam")
    ped = ped.iloc[:, :3]
    ped.columns = ["iid", "sire", "dam"]
    # A header line survives the read as a row; drop it if it self-describes.
    if ped.iloc[0]["iid"].lower() in {"iid", "id", "individ_id", "animal"}:
        ped = ped.iloc[1:]
    return ped.drop_duplicates(subset="iid")


def read_id_file(path: Path) -> set[str]:
    """Accept either a bare ID per line or PLINK's FID/IID pair."""
    ids: set[str] = set()
    for line in Path(path).read_text().splitlines():
        fields = line.split()
        if not fields:
            continue
        ids.add(fields[1] if len(fields) >= 2 else fields[0])
    return ids


def family_key(row) -> str:
    """The sire x dam pair, or a per-animal key when the parents are unknown.

    Animals with no recorded parents must not collapse into one enormous
    pseudo-family: as far as the data says they are unrelated, and capping them
    as a group would throw away most of the founders.
    """
    sire, dam = row.get("sire", "0"), row.get("dam", "0")
    if not sire or not dam or sire in {"0", "NA", ""} or dam in {"0", "NA", ""}:
        return f"__unknown__{row['iid']}"
    return f"{sire}\x00{dam}"


def select(
    fam: pd.DataFrame,
    pedigree: pd.DataFrame | None,
    excluded: set[str],
    max_per_family: int,
    seed: int,
) -> pd.DataFrame:
    """Return the fam with `kept` and `reason` columns attached."""
    panel = fam[["fid", "iid"]].copy()

    if pedigree is not None:
        panel = panel.merge(pedigree, on="iid", how="left")
    else:
        panel["sire"] = "0"
        panel["dam"] = "0"
    panel[["sire", "dam"]] = panel[["sire", "dam"]].fillna("0")

    panel["family"] = panel.apply(family_key, axis=1)
    panel["excluded"] = panel["iid"].isin(excluded)

    panel["kept"] = ~panel["excluded"]
    panel["reason"] = ""
    panel.loc[panel["excluded"], "reason"] = "held out by exclusion list"

    if max_per_family and max_per_family > 0:
        eligible = panel.loc[panel["kept"]]
        # Deterministic: shuffle once at the seed, then take the first N of
        # each family in that order.
        shuffled = eligible.sample(frac=1.0, random_state=seed)
        rank = shuffled.groupby("family").cumcount()
        over = shuffled.loc[rank >= max_per_family, "iid"]
        panel.loc[panel["iid"].isin(set(over)), ["kept", "reason"]] = [
            False,
            f"over the cap of {max_per_family} per full-sib family",
        ]

    panel.loc[panel["kept"], "reason"] = "selected"
    return panel


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fam", required=True, type=Path)
    parser.add_argument("--pedigree", default="", help="TSV of iid sire dam; enables family capping")
    parser.add_argument("--exclude", nargs="*", default=[],
                        help="ID files whose animals are held out before phasing")
    parser.add_argument("--max-per-family", type=int, default=0,
                        help="Cap per sire x dam pair; 0 keeps everyone")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--ids-out", required=True, type=Path)
    parser.add_argument("--manifest-out", required=True, type=Path)
    args = parser.parse_args()

    fam = read_fam(args.fam)
    pedigree = read_pedigree(Path(args.pedigree)) if str(args.pedigree).strip() else None
    if args.max_per_family and pedigree is None:
        sys.exit(
            "ERROR: --max-per-family needs --pedigree. A PLINK .fam straight off an "
            "array export carries no parents, so there are no families to cap; "
            "resolve them from GPA first."
        )

    excluded: set[str] = set()
    for path in args.exclude:
        if str(path).strip():
            excluded |= read_id_file(Path(path))

    panel = select(fam, pedigree, excluded, args.max_per_family, args.seed)
    kept = panel.loc[panel["kept"]]

    if kept.empty:
        sys.exit("ERROR: the selection kept no animals.")

    for path in (args.ids_out, args.manifest_out):
        path.parent.mkdir(parents=True, exist_ok=True)
    kept[["fid", "iid"]].to_csv(args.ids_out, sep="\t", header=False, index=False)
    panel.drop(columns=["family"]).to_csv(args.manifest_out, sep="\t", index=False)

    named = panel.loc[~panel["family"].str.startswith("__unknown__")]
    print(f"animals_in_fileset={len(panel)}")
    print(f"excluded_by_list={int(panel['excluded'].sum())}")
    print(f"dropped_over_family_cap={int((~panel['kept'] & ~panel['excluded']).sum())}")
    print(f"selected={len(kept)}")
    print(f"families_with_known_parents={named['family'].nunique()}")
    print(f"animals_without_known_parents={len(panel) - len(named)}")


if __name__ == "__main__":
    main()
