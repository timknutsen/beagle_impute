#!/usr/bin/env python3
"""FImpute text input/output helpers for accuracy benchmarking."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd


def _open_text(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def _sex_to_fimpute(value: object) -> str:
    text = str(value)
    if text == "1":
        return "M"
    if text == "2":
        return "F"
    return "U"


FIMPUTE_MISSING = 5


def _raw_value_to_fimpute(value: object) -> str:
    """Single-value form of the raw → FImpute mapping (kept for readability/tests)."""
    if pd.isna(value) or str(value) == "NA":
        return str(FIMPUTE_MISSING)
    return str(int(float(value)))


def _raw_block_to_fimpute_strings(block: pd.DataFrame) -> list[str]:
    """
    Vectorised equivalent of joining _raw_value_to_fimpute over every column.

    The per-row/per-column form is O(n_animals x n_markers) Python-level
    operations, which on a real dataset (10k animals, 50k markers) is hundreds
    of millions of pandas lookups. Here the whole block is converted once in
    numpy and sliced into fixed-width ASCII rows instead.
    """
    values = block.to_numpy(dtype="float64", na_value=np.nan)
    codes = np.where(np.isnan(values), FIMPUTE_MISSING, values).astype(np.uint8)
    # FImpute codes are single digits, so the byte value is just the ASCII digit.
    ascii_rows = (codes + ord("0")).tobytes()
    width = codes.shape[1]
    return [
        ascii_rows[i * width : (i + 1) * width].decode("ascii")
        for i in range(codes.shape[0])
    ]


def build_panel_map(bims: list[pd.DataFrame]) -> pd.DataFrame:
    """
    Merge one .bim per chip into FImpute's SNP info table.

    Returns SNP_ID / Chr / Pos plus one Chip<N> column per panel, holding that
    marker's 1-based index *within that panel* and 0 when the panel does not
    carry it. The union is ordered by position, which is the order every
    genotype string must follow.
    """
    union = (
        pd.concat([bim[["snp", "chrom", "pos"]] for bim in bims])
        .drop_duplicates(subset="snp")
        .sort_values("pos")
        .reset_index(drop=True)
    )

    out = pd.DataFrame(
        {
            "SNP_ID": union["snp"].astype(str),
            "Chr": union["chrom"].astype(str),
            "Pos": union["pos"].astype(int),
        }
    )

    for idx, bim in enumerate(bims, start=1):
        # Index within this panel follows the panel's own position order.
        panel_order = bim.sort_values("pos")["snp"].astype(str).tolist()
        rank = {snp: i + 1 for i, snp in enumerate(panel_order)}
        out[f"Chip{idx}"] = [rank.get(snp, 0) for snp in out["SNP_ID"]]

    return out


def _panel_genotype_strings(raw: pd.DataFrame, bim: pd.DataFrame) -> list[str]:
    """Encode one panel's animals, ordered by that panel's position order."""
    # .raw columns are "<snp>_<allele>", so match on the leading SNP id.
    geno_cols = list(raw.columns[6:])
    by_snp = {col.rsplit("_", 1)[0]: col for col in geno_cols}
    ordered = [
        by_snp[snp]
        for snp in bim.sort_values("pos")["snp"].astype(str)
        if snp in by_snp
    ]
    if len(ordered) != len(bim):
        raise ValueError(
            f"RAW carries {len(ordered)} of the panel's {len(bim)} markers; "
            "the .raw and .bim do not describe the same panel"
        )
    return _raw_block_to_fimpute_strings(raw[ordered])


def write_fimpute_inputs_from_raw(
    raw_path: str | Path,
    bim_path: str | Path,
    fam_path: str | Path,
    out_dir: str | Path,
    chrom: str,
    nthreads: int,
    ref_raw_path: str | Path | None = None,
    ref_bim_path: str | Path | None = None,
    ref_fam_path: str | Path | None = None,
    keep_observed: bool = True,
    phase_output: bool = True,
) -> dict[str, Path]:
    """
    Write FImpute .genos, .snps, .ped, .ctrl and short ID map files.

    With ref_* supplied the run becomes a two-panel imputation: chip 1 is the
    reference (typically the denser array, declared to FImpute via ref_chip=1)
    and chip 2 is the target. Without them a single panel is written, which
    still fills sporadic missing calls and phases.
    """
    raw_path = Path(raw_path)
    bim_path = Path(bim_path)
    fam_path = Path(fam_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(raw_path, sep=r"\s+")
    bim = pd.read_csv(
        bim_path,
        sep=r"\s+",
        header=None,
        names=["chrom", "snp", "cm", "pos", "a1", "a2"],
    )
    fam = pd.read_csv(
        fam_path,
        sep=r"\s+",
        header=None,
        names=["fid", "iid", "pat", "mat", "sex", "phen"],
    )

    use_reference = ref_raw_path is not None
    if use_reference:
        ref_raw = pd.read_csv(ref_raw_path, sep=r"\s+")
        ref_bim = pd.read_csv(
            ref_bim_path,
            sep=r"\s+",
            header=None,
            names=["chrom", "snp", "cm", "pos", "a1", "a2"],
        )
        ref_fam = pd.read_csv(
            ref_fam_path if ref_fam_path is not None else fam_path,
            sep=r"\s+",
            header=None,
            names=["fid", "iid", "pat", "mat", "sex", "phen"],
        )
    else:
        ref_raw = ref_bim = ref_fam = None

    # Chip 1 is the reference when one is supplied, so it is listed first and
    # named as ref_chip below. Animal order follows the same convention.
    if use_reference:
        all_fam = pd.concat([ref_fam, fam]).drop_duplicates(subset="iid")
    else:
        all_fam = fam

    short_ids = {iid: str(idx + 1) for idx, iid in enumerate(all_fam["iid"].astype(str))}
    id_map = all_fam[["fid", "iid"]].copy()
    id_map.insert(0, "short_id", [short_ids[iid] for iid in all_fam["iid"].astype(str)])

    if use_reference:
        snps = build_panel_map([ref_bim, bim])
        genos = pd.DataFrame(
            {
                "ID": (
                    [short_ids[i] for i in ref_raw["IID"].astype(str)]
                    + [short_ids[i] for i in raw["IID"].astype(str)]
                ),
                "Chip": ["1"] * len(ref_raw) + ["2"] * len(raw),
                "Genotypes": (
                    _panel_genotype_strings(ref_raw, ref_bim)
                    + _panel_genotype_strings(raw, bim)
                ),
            }
        )
    else:
        snps = build_panel_map([bim])
        genos = pd.DataFrame(
            {
                "ID": [short_ids[iid] for iid in raw["IID"].astype(str)],
                "Chip": "1",
                "Genotypes": _panel_genotype_strings(raw, bim),
            }
        )

    parent_to_short = lambda value: short_ids.get(str(value), "0")
    ped = pd.DataFrame(
        {
            "ID": [short_ids[iid] for iid in all_fam["iid"].astype(str)],
            "sire ID": [parent_to_short(value) for value in all_fam["pat"]],
            "dam ID": [parent_to_short(value) for value in all_fam["mat"]],
            "sex": [_sex_to_fimpute(value) for value in all_fam["sex"]],
        }
    )

    genos_path = out_dir / f"chr{chrom}.genos"
    snps_path = out_dir / f"chr{chrom}.snps"
    ped_path = out_dir / f"chr{chrom}.ped"
    ctrl_path = out_dir / f"chr{chrom}.ctrl"
    id_map_path = out_dir / f"chr{chrom}.id_map.tsv"

    genos.to_csv(genos_path, sep="\t", index=False)
    snps.to_csv(snps_path, sep="\t", index=False)
    ped.to_csv(ped_path, sep="\t", index=False)
    id_map.to_csv(id_map_path, sep="\t", index=False)

    output_folder = out_dir / f"chr{chrom}"
    output_folder.mkdir(parents=True, exist_ok=True)

    lines = [
        f'title="chr{chrom}";',
        f'genotype_file="{genos_path}";',
        f'snp_info_file="{snps_path}";',
        f'ped_file="{ped_path}";',
        f'output_folder="{output_folder}";',
        # Parentage testing runs by default and costs a pass over every animal
        # for a result this pipeline does not consume.
        "parentage_test /off;",
        "njob=1;",
        f"nthread={nthreads};",
    ]
    if not phase_output:
        # save_genotype collapses the output to plain 0/1/2/5 calls. That also
        # discards phase: with it set FImpute reports heterozygotes as code 1,
        # without it as the resolved codes 3 and 4. So it is only requested
        # when the caller has said it does not want phase.
        lines.insert(-2, "save_genotype;")
    if keep_observed:
        # Leave genotypes that were actually measured untouched.
        lines.insert(-2, "keep_og;")
    if use_reference:
        lines.insert(-2, "ref_chip=1;")

    ctrl_path.write_text("\n".join(lines) + "\n")

    return {
        "genos": genos_path,
        "snps": snps_path,
        "ped": ped_path,
        "ctrl": ctrl_path,
        "id_map": id_map_path,
    }


# FImpute reports heterozygotes as 3 or 4 once it has resolved phase: 3 is
# "reference allele on the paternal strand", 4 is the reverse. Code 1 is an
# unphased heterozygote. Writing 3/4 with "/" would silently discard the phase
# FImpute computed -- in practice it resolves nearly every het, so that is most
# of the phasing information in the run.
_GT_UNPHASED = {
    "0": "0/0",
    "1": "0/1",
    "2": "1/1",
    "3": "0/1",
    "4": "0/1",
    "5": "./.",
}

_GT_PHASED = {
    "0": "0|0",
    "1": "0/1",   # genuinely unphased; keep it that way
    "2": "1|1",
    "3": "0|1",
    "4": "1|0",
    "5": "./.",
}


def fimpute_calls_to_vcf_gt(code: str, phased: bool = False) -> str:
    """
    Map an FImpute single-character call to a VCF GT string.

    With phased=True the resolved heterozygote codes keep their phase
    ("0|1" / "1|0"); otherwise every heterozygote collapses to "0/1".
    """
    mapping = _GT_PHASED if phased else _GT_UNPHASED
    try:
        return mapping[str(code)]
    except KeyError as exc:
        raise ValueError(f"Unsupported FImpute genotype code: {code!r}") from exc


def write_fimpute_vcf(
    imputed_path: str | Path,
    bim_path: str | Path | list,
    id_map_path: str | Path,
    out_vcf: str | Path,
    snp_info_path: str | Path | None = None,
    phased: bool = True,
) -> None:
    """
    Convert FImpute genotypes_imp.txt to a VCF with IID sample names.

    bim_path may be several .bim files when more than one chip was used; they
    supply REF/ALT. With snp_info_path given, marker order follows that SNP
    info table (the union across chips) rather than a single .bim, which is
    what a two-panel run returns.
    """
    imputed_path = Path(imputed_path)
    bim_paths = [Path(p) for p in (bim_path if isinstance(bim_path, list) else [bim_path])]
    id_map_path = Path(id_map_path)
    out_vcf = Path(out_vcf)
    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    calls = pd.read_csv(imputed_path, sep="\t", dtype={"ID": str, "Calls...": str})
    bims = [
        pd.read_csv(
            p, sep=r"\s+", header=None,
            names=["chrom", "snp", "cm", "pos", "a1", "a2"],
        )
        for p in bim_paths
    ]
    # Alleles come from whichever chip carries the marker; identical SNP ids
    # across chips must agree, so the first one wins.
    alleles = pd.concat(bims).drop_duplicates(subset=["snp"]).set_index("snp")

    if snp_info_path is not None:
        snp_info = pd.read_csv(snp_info_path, sep="\t")
        order = snp_info["SNP_ID"].astype(str).tolist()
    else:
        order = bims[0]["snp"].astype(str).tolist()

    id_map = pd.read_csv(id_map_path, sep="\t", dtype={"short_id": str, "fid": str, "iid": str})
    sample_map = dict(zip(id_map["short_id"], id_map["iid"]))

    samples = [sample_map[str(short_id)] for short_id in calls["ID"].astype(str)]
    call_strings = calls["Calls..."].astype(str).tolist()

    for short_id, string in zip(calls["ID"].astype(str), call_strings):
        if len(string) != len(order):
            raise ValueError(
                f"Animal {short_id} has {len(string)} calls but the marker map "
                f"lists {len(order)}; genotype string and SNP map disagree"
            )

    with _open_text(out_vcf, "wt") as out:
        contigs = sorted(
            {str(c) for c in alleles.loc[order, "chrom"]},
            key=lambda value: int(value) if str(value).isdigit() else 1 << 30,
        )
        out.write("##fileformat=VCFv4.2\n")
        for chrom in contigs:
            out.write(f"##contig=<ID={chrom}>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        if samples:
            out.write("\t" + "\t".join(samples))
        out.write("\n")

        for idx, snp in enumerate(order):
            row = alleles.loc[snp]
            gts = [
                fimpute_calls_to_vcf_gt(call_string[idx], phased=phased)
                for call_string in call_strings
            ]
            out.write(
                f"{row.chrom}\t{int(row.pos)}\t{snp}\t{row.a2}\t{row.a1}"
                f"\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prep = subparsers.add_parser("prepare-inputs")
    prep.add_argument("--raw", required=True, type=Path)
    prep.add_argument("--bim", required=True, type=Path)
    prep.add_argument("--fam", required=True, type=Path)
    prep.add_argument("--out-dir", required=True, type=Path)
    prep.add_argument("--chrom", required=True)
    prep.add_argument("--nthreads", required=True, type=int)
    prep.add_argument("--ref-raw", type=Path, default=None,
                      help="Reference panel .raw; enables two-chip imputation")
    prep.add_argument("--ref-bim", type=Path, default=None)
    prep.add_argument("--ref-fam", type=Path, default=None)
    prep.add_argument("--no-phase", action="store_true",
                      help="Ask FImpute for plain 0/1/2/5 calls (drops phase)")
    prep.add_argument("--no-keep-observed", action="store_true",
                      help="Let FImpute overwrite genotypes that were measured")

    to_vcf = subparsers.add_parser("to-vcf")
    to_vcf.add_argument("--imputed", required=True, type=Path)
    to_vcf.add_argument("--bim", required=True, type=Path, nargs="+")
    to_vcf.add_argument("--snp-info", type=Path, default=None,
                        help="SNP info table giving marker order across chips")
    to_vcf.add_argument("--unphased", action="store_true",
                        help="Collapse resolved heterozygotes to 0/1")
    to_vcf.add_argument("--id-map", required=True, type=Path)
    to_vcf.add_argument("--out-vcf", required=True, type=Path)

    args = parser.parse_args()
    if args.command == "prepare-inputs":
        write_fimpute_inputs_from_raw(
            args.raw, args.bim, args.fam, args.out_dir, args.chrom, args.nthreads,
            ref_raw_path=args.ref_raw,
            ref_bim_path=args.ref_bim,
            ref_fam_path=args.ref_fam,
            keep_observed=not args.no_keep_observed,
            phase_output=not args.no_phase,
        )
    elif args.command == "to-vcf":
        write_fimpute_vcf(
            args.imputed, args.bim, args.id_map, args.out_vcf,
            snp_info_path=args.snp_info,
            phased=not args.unphased,
        )


if __name__ == "__main__":
    main()
