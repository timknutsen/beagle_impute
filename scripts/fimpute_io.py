#!/usr/bin/env python3
"""FImpute text input/output helpers for accuracy benchmarking."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

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


def _raw_value_to_fimpute(value: object) -> str:
    if pd.isna(value) or str(value) == "NA":
        return "5"
    return str(int(float(value)))


def write_fimpute_inputs_from_raw(
    raw_path: str | Path,
    bim_path: str | Path,
    fam_path: str | Path,
    out_dir: str | Path,
    chrom: str,
    nthreads: int,
) -> dict[str, Path]:
    """Write FImpute .genos, .snps, .ped, .ctrl and short ID map files."""
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

    short_ids = {iid: str(idx + 1) for idx, iid in enumerate(fam["iid"].astype(str))}
    id_map = fam[["fid", "iid"]].copy()
    id_map.insert(0, "short_id", [short_ids[iid] for iid in fam["iid"].astype(str)])

    geno_cols = list(raw.columns[6:])
    if len(geno_cols) != len(bim):
        raise ValueError(
            f"RAW marker count ({len(geno_cols)}) does not match BIM marker count ({len(bim)})"
        )

    genos = pd.DataFrame(
        {
            "ID": [short_ids[iid] for iid in raw["IID"].astype(str)],
            "Chip": "1",
            "Genotypes": [
                "".join(_raw_value_to_fimpute(row[col]) for col in geno_cols)
                for _, row in raw.iterrows()
            ],
        }
    )

    parent_to_short = lambda value: short_ids.get(str(value), "0")
    ped = pd.DataFrame(
        {
            "ID": [short_ids[iid] for iid in fam["iid"].astype(str)],
            "sire ID": [parent_to_short(value) for value in fam["pat"]],
            "dam ID": [parent_to_short(value) for value in fam["mat"]],
            "sex": [_sex_to_fimpute(value) for value in fam["sex"]],
        }
    )
    snps = pd.DataFrame(
        {
            "SNP_ID": bim["snp"].astype(str),
            "Chr": str(chrom),
            "Pos": bim["pos"].astype(int),
            "Chip1": list(range(1, len(bim) + 1)),
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
    ctrl_path.write_text(
        f'genotype_file="{genos_path}";\n'
        f'snp_info_file="{snps_path}";\n'
        f'output_folder="{output_folder}";\n'
        "njob=1;\n"
        f"nthread={nthreads};\n"
        f'ped_file="{ped_path}";\n'
    )

    return {
        "genos": genos_path,
        "snps": snps_path,
        "ped": ped_path,
        "ctrl": ctrl_path,
        "id_map": id_map_path,
    }


def fimpute_calls_to_vcf_gt(code: str) -> str:
    """Map FImpute single-character calls to VCF GT strings."""
    mapping = {
        "0": "0/0",
        "1": "0/1",
        "2": "1/1",
        "3": "0/1",
        "4": "1/0",
        "5": "./.",
    }
    try:
        return mapping[str(code)]
    except KeyError as exc:
        raise ValueError(f"Unsupported FImpute genotype code: {code!r}") from exc


def write_fimpute_vcf(
    imputed_path: str | Path,
    bim_path: str | Path,
    id_map_path: str | Path,
    out_vcf: str | Path,
) -> None:
    """Convert FImpute genotypes_imp.txt to a VCF with IID sample names."""
    imputed_path = Path(imputed_path)
    bim_path = Path(bim_path)
    id_map_path = Path(id_map_path)
    out_vcf = Path(out_vcf)
    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    calls = pd.read_csv(imputed_path, sep="\t", dtype={"ID": str, "Calls...": str})
    bim = pd.read_csv(
        bim_path,
        sep=r"\s+",
        header=None,
        names=["chrom", "snp", "cm", "pos", "a1", "a2"],
    )
    id_map = pd.read_csv(id_map_path, sep="\t", dtype={"short_id": str, "fid": str, "iid": str})
    sample_map = dict(zip(id_map["short_id"], id_map["iid"]))

    samples = [sample_map[str(short_id)] for short_id in calls["ID"].astype(str)]

    with _open_text(out_vcf, "wt") as out:
        contigs = sorted({str(chrom) for chrom in bim["chrom"]}, key=lambda value: int(value))
        out.write("##fileformat=VCFv4.2\n")
        for chrom in contigs:
            out.write(f"##contig=<ID={chrom}>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        if samples:
            out.write("\t" + "\t".join(samples))
        out.write("\n")

        call_strings = calls["Calls..."].astype(str).tolist()
        for idx, row in bim.reset_index(drop=True).iterrows():
            gts = [fimpute_calls_to_vcf_gt(call_string[idx]) for call_string in call_strings]
            out.write(
                f"{row.chrom}\t{int(row.pos)}\t{row.snp}\t{row.a2}\t{row.a1}"
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

    to_vcf = subparsers.add_parser("to-vcf")
    to_vcf.add_argument("--imputed", required=True, type=Path)
    to_vcf.add_argument("--bim", required=True, type=Path)
    to_vcf.add_argument("--id-map", required=True, type=Path)
    to_vcf.add_argument("--out-vcf", required=True, type=Path)

    args = parser.parse_args()
    if args.command == "prepare-inputs":
        write_fimpute_inputs_from_raw(args.raw, args.bim, args.fam, args.out_dir, args.chrom, args.nthreads)
    elif args.command == "to-vcf":
        write_fimpute_vcf(args.imputed, args.bim, args.id_map, args.out_vcf)


if __name__ == "__main__":
    main()
