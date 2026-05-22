#!/usr/bin/env python3
"""Mask validation genotypes outside a low-density SNP panel in a PLINK BED set."""

from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path


BED_MAGIC = b"\x6c\x1b\x01"
MISSING_BITS = 0b01


def read_ids(path: Path) -> set[tuple[str, str]]:
    ids: set[tuple[str, str]] = set()
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"Expected FID IID columns in {path}, got: {line!r}")
            ids.add((fields[0], fields[1]))
    if not ids:
        raise ValueError(f"No validation IDs found in {path}")
    return ids


def read_snp_ids(path: Path) -> set[str]:
    snps = {line.split()[0] for line in path.read_text().splitlines() if line.strip()}
    if not snps:
        raise ValueError(f"No SNP IDs found in {path}")
    return snps


def read_fam_ids(path: Path) -> list[tuple[str, str]]:
    ids: list[tuple[str, str]] = []
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"Expected PLINK .fam with at least two columns, got: {line!r}")
            ids.append((fields[0], fields[1]))
    if not ids:
        raise ValueError(f"No samples found in {path}")
    return ids


def read_bim_snp_ids(path: Path) -> list[str]:
    snps: list[str] = []
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"Expected PLINK .bim with at least two columns, got: {line!r}")
            snps.append(fields[1])
    if not snps:
        raise ValueError(f"No variants found in {path}")
    return snps


def set_missing(byte_value: int, sample_offset: int) -> int:
    shift = 2 * sample_offset
    return (byte_value & ~(0b11 << shift)) | (MISSING_BITS << shift)


def mask_bed(
    in_bed: Path,
    out_bed: Path,
    bim_snps: list[str],
    ld_snps: set[str],
    validation_sample_indexes: list[int],
    n_samples: int,
) -> tuple[int, int]:
    bytes_per_variant = math.ceil(n_samples / 4)
    masked_variants = 0
    masked_genotypes = 0

    with in_bed.open("rb") as src, out_bed.open("wb") as dst:
        magic = src.read(3)
        if magic != BED_MAGIC:
            raise ValueError(
                f"{in_bed} is not a SNP-major PLINK BED file; expected magic {BED_MAGIC!r}, got {magic!r}"
            )
        dst.write(magic)

        for snp_id in bim_snps:
            block = bytearray(src.read(bytes_per_variant))
            if len(block) != bytes_per_variant:
                raise ValueError(f"Unexpected end of BED file while reading SNP {snp_id}")

            if snp_id not in ld_snps:
                masked_variants += 1
                for sample_index in validation_sample_indexes:
                    byte_index, sample_offset = divmod(sample_index, 4)
                    block[byte_index] = set_missing(block[byte_index], sample_offset)
                    masked_genotypes += 1
            dst.write(block)

        extra = src.read(1)
        if extra:
            raise ValueError(f"{in_bed} has trailing bytes after {len(bim_snps)} variants")

    return masked_variants, masked_genotypes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bfile", required=True, help="Input PLINK BED prefix")
    parser.add_argument("--validation-ids", required=True, type=Path, help="FID/IID validation animals")
    parser.add_argument("--ld-snps", required=True, type=Path, help="Low-density SNP IDs to preserve")
    parser.add_argument("--out", required=True, help="Output PLINK BED prefix")
    args = parser.parse_args()

    in_prefix = Path(args.bfile)
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    fam_path = Path(f"{in_prefix}.fam")
    bim_path = Path(f"{in_prefix}.bim")
    bed_path = Path(f"{in_prefix}.bed")

    fam_ids = read_fam_ids(fam_path)
    validation_ids = read_ids(args.validation_ids)
    validation_indexes = [i for i, sample_id in enumerate(fam_ids) if sample_id in validation_ids]
    missing_validation = validation_ids.difference(fam_ids)
    if missing_validation:
        examples = ", ".join(f"{fid}:{iid}" for fid, iid in sorted(missing_validation)[:5])
        raise ValueError(f"{len(missing_validation)} validation IDs were not found in {fam_path}: {examples}")

    ld_snps = read_snp_ids(args.ld_snps)
    bim_snps = read_bim_snp_ids(bim_path)
    missing_snps = ld_snps.difference(bim_snps)
    if missing_snps:
        examples = ", ".join(sorted(missing_snps)[:5])
        raise ValueError(f"{len(missing_snps)} LD SNP IDs were not found in {bim_path}: {examples}")

    shutil.copy2(bim_path, Path(f"{out_prefix}.bim"))
    shutil.copy2(fam_path, Path(f"{out_prefix}.fam"))
    masked_variants, masked_genotypes = mask_bed(
        bed_path,
        Path(f"{out_prefix}.bed"),
        bim_snps,
        ld_snps,
        validation_indexes,
        len(fam_ids),
    )

    print(f"samples={len(fam_ids)}")
    print(f"validation_samples={len(validation_indexes)}")
    print(f"variants={len(bim_snps)}")
    print(f"ld_variants={len(ld_snps)}")
    print(f"masked_variants={masked_variants}")
    print(f"masked_genotypes={masked_genotypes}")


if __name__ == "__main__":
    main()
