#!/usr/bin/env python3
"""Mask validation genotypes outside a low-density SNP panel in a PLINK BED set.

With --replace-from the validation animals additionally get their genotypes at
the panel markers taken from a second fileset. That is what a cross-array test
needs: fold k's calls must come from the low-density array it was really typed
on, not from the high-density array's calls at the same positions. Reusing the
HD calls would remove the array-transition effect the test exists to measure --
probe differences and between-array genotyping error -- and quietly turn a
cross-array benchmark into a masking benchmark.
"""

from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path


BED_MAGIC = b"\x6c\x1b\x01"
MISSING_BITS = 0b01

# PLINK1 .bed two-bit codes, sample-major within each variant block.
HOM_A1 = 0b00
MISSING = 0b01
HET = 0b10
HOM_A2 = 0b11
_FLIP = {HOM_A1: HOM_A2, HOM_A2: HOM_A1, HET: HET, MISSING: MISSING}


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



def read_bim_alleles(path: Path) -> dict[str, tuple[int, str, str]]:
    """Marker ID -> (row index, a1, a2), for orienting a second fileset."""
    alleles: dict[str, tuple[int, str, str]] = {}
    with path.open() as handle:
        for index, line in enumerate(handle):
            fields = line.split()
            if len(fields) < 6:
                raise ValueError(f"Expected a six-column PLINK .bim, got: {line!r}")
            alleles[fields[1]] = (index, fields[4], fields[5])
    if not alleles:
        raise ValueError(f"No variants found in {path}")
    return alleles


def flip_needed(snp_id: str, source: tuple[str, str], target: tuple[str, str]) -> bool:
    """Whether the source fileset counts the opposite allele to the target.

    An undetected swap mirrors every call at that marker, and allelic r2 cannot
    see it -- squaring the correlation hides the sign -- so it surfaces only as
    a collapsed concordance much later. Refuse anything that is not cleanly the
    same pair or the same pair reversed.
    """
    src_a1, src_a2 = source
    tgt_a1, tgt_a2 = target
    if (src_a1, src_a2) == (tgt_a1, tgt_a2):
        return False
    if (src_a1, src_a2) == (tgt_a2, tgt_a1):
        return True
    # PLINK writes "0" for the absent allele of a monomorphic marker, so only
    # the observed one can be matched.
    if src_a1 == "0":
        if src_a2 == tgt_a2:
            return False
        if src_a2 == tgt_a1:
            return True
    if tgt_a1 == "0":
        if src_a2 == tgt_a2:
            return False
        if src_a1 == tgt_a2:
            return True
    raise ValueError(
        f"Marker {snp_id} carries alleles {src_a1}/{src_a2} in the replacement fileset "
        f"but {tgt_a1}/{tgt_a2} in the target. These are not the same marker; "
        "harmonise the two arrays upstream (plink2 --ref-allele force) before splicing."
    )


class ReplacementSource:
    """Random-access reader for one variant of a second PLINK fileset.

    Only the validation animals and only the panel markers are ever read, so
    the file is seeked rather than streamed: the panel is a small slice of a
    large .bed and the animals are a small slice of each variant block.
    """

    def __init__(self, prefix: Path, validation_iids: list[str], target_alleles: dict[str, tuple[str, str]]):
        self.bed_path = Path(f"{prefix}.bed")
        fam = read_fam_ids(Path(f"{prefix}.fam"))
        self.n_samples = len(fam)
        self.bytes_per_variant = math.ceil(self.n_samples / 4)

        index_by_iid: dict[str, int] = {}
        for position, (_fid, iid) in enumerate(fam):
            if iid in index_by_iid:
                raise ValueError(f"Duplicate individual ID {iid} in {prefix}.fam")
            index_by_iid[iid] = position
        missing = [iid for iid in validation_iids if iid not in index_by_iid]
        if missing:
            raise ValueError(
                f"{len(missing)} validation animals are absent from {prefix}.fam: "
                + ", ".join(missing[:5])
            )
        self.sample_positions = [index_by_iid[iid] for iid in validation_iids]

        self.variants = read_bim_alleles(Path(f"{prefix}.bim"))
        self.flip = {
            snp_id: flip_needed(snp_id, (a1, a2), target_alleles[snp_id])
            for snp_id, (_index, a1, a2) in self.variants.items()
            if snp_id in target_alleles
        }
        self._handle = None

    def __enter__(self) -> "ReplacementSource":
        self._handle = self.bed_path.open("rb")
        magic = self._handle.read(3)
        if magic != BED_MAGIC:
            raise ValueError(f"{self.bed_path} is not a SNP-major PLINK BED file")
        return self

    def __exit__(self, *exc) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None

    def has(self, snp_id: str) -> bool:
        return snp_id in self.variants

    def codes(self, snp_id: str) -> list[int]:
        """Two-bit codes for the validation animals, oriented to the target."""
        variant_index = self.variants[snp_id][0]
        self._handle.seek(3 + variant_index * self.bytes_per_variant)
        block = self._handle.read(self.bytes_per_variant)
        if len(block) != self.bytes_per_variant:
            raise ValueError(f"Unexpected end of {self.bed_path} at marker {snp_id}")
        flip = self.flip.get(snp_id, False)
        out = []
        for position in self.sample_positions:
            byte_index, offset = divmod(position, 4)
            code = (block[byte_index] >> (2 * offset)) & 0b11
            out.append(_FLIP[code] if flip else code)
        return out


def set_code(byte_value: int, sample_offset: int, code: int) -> int:
    shift = 2 * sample_offset
    return (byte_value & ~(0b11 << shift)) | (code << shift)


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
    replacement: "ReplacementSource | None" = None,
) -> dict[str, int]:
    bytes_per_variant = math.ceil(n_samples / 4)
    counts = {
        "masked_variants": 0,
        "masked_genotypes": 0,
        "replaced_variants": 0,
        "replaced_genotypes": 0,
        "panel_variants_without_replacement": 0,
    }

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
                counts["masked_variants"] += 1
                for sample_index in validation_sample_indexes:
                    byte_index, sample_offset = divmod(sample_index, 4)
                    block[byte_index] = set_missing(block[byte_index], sample_offset)
                    counts["masked_genotypes"] += 1
            elif replacement is not None:
                # A panel marker the replacement fileset does not carry has no
                # observed call to offer, so it is masked rather than left as
                # the target array's own genotype -- which would hand the
                # imputer an answer it was never given in the real transition.
                if not replacement.has(snp_id):
                    counts["panel_variants_without_replacement"] += 1
                    for sample_index in validation_sample_indexes:
                        byte_index, sample_offset = divmod(sample_index, 4)
                        block[byte_index] = set_missing(block[byte_index], sample_offset)
                        counts["masked_genotypes"] += 1
                else:
                    counts["replaced_variants"] += 1
                    for sample_index, code in zip(
                        validation_sample_indexes, replacement.codes(snp_id)
                    ):
                        byte_index, sample_offset = divmod(sample_index, 4)
                        block[byte_index] = set_code(block[byte_index], sample_offset, code)
                        counts["replaced_genotypes"] += 1
            dst.write(block)

        extra = src.read(1)
        if extra:
            raise ValueError(f"{in_bed} has trailing bytes after {len(bim_snps)} variants")

    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bfile", required=True, help="Input PLINK BED prefix")
    parser.add_argument("--validation-ids", required=True, type=Path, help="FID/IID validation animals")
    parser.add_argument("--ld-snps", required=True, type=Path, help="Low-density SNP IDs to preserve")
    parser.add_argument("--out", required=True, help="Output PLINK BED prefix")
    parser.add_argument(
        "--replace-from",
        default=None,
        help=(
            "PLINK prefix of the array the validation animals were really typed on. "
            "Their calls at the panel markers are taken from there instead of being "
            "left as the target array's own."
        ),
    )
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

    if args.replace_from:
        target_alleles = {
            snp_id: (a1, a2)
            for snp_id, (_index, a1, a2) in read_bim_alleles(bim_path).items()
            if snp_id in ld_snps
        }
        validation_iids = [fam_ids[i][1] for i in validation_indexes]
        with ReplacementSource(
            Path(args.replace_from), validation_iids, target_alleles
        ) as source:
            counts = mask_bed(
                bed_path,
                Path(f"{out_prefix}.bed"),
                bim_snps,
                ld_snps,
                validation_indexes,
                len(fam_ids),
                replacement=source,
            )
    else:
        counts = mask_bed(
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
    for key, value in counts.items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()
