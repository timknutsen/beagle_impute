"""
Shared pytest fixtures for beagle_impute tests.

Provides:
  synth_plink   — tiny synthetic PLINK binary dataset with pedigree structure,
                  generated in-memory (no binary files stored in git).
  synth_ai2     — corresponding AlphaImpute2 genotype + pedigree text files.

Synthetic dataset layout
------------------------
Chromosomes : 2  (chr1, chr2)
SNPs        : 40 per chromosome = 80 total
Individuals : 12
  sire1, sire2  — founder males   (fully genotyped)
  dam1,  dam2   — founder females (fully genotyped)
  off01–off08   — offspring with known parents, 25 % missingness introduced
"""

import random
import sys
from pathlib import Path

import pytest

# Make scripts/ importable so tests can import conversion helpers directly.
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

REPO_ROOT = Path(__file__).parent.parent

# ── Synthetic dataset parameters ─────────────────────────────────────────────
_N_SNPS_PER_CHROM = 40
_N_CHROMS = 2
_MISSING_RATE = 0.25
_SEED = 42

# PLINK .bed genotype encoding (SNP-major, 2 bits per individual, LSB first):
#   00 = hom A1 (ALT in VCF)  ← dosage 2 from plink2 --export A
#   01 = missing
#   10 = het                   ← dosage 1
#   11 = hom A2 (REF in VCF)  ← dosage 0
_BED_ENC = {0: 0b11, 1: 0b10, 2: 0b00, 9: 0b01}


# ── Helpers ───────────────────────────────────────────────────────────────────

def _build_pedigree():
    """Return list of (FID, IID, PAT, MAT, SEX, PHENO) tuples."""
    founders = [
        ("FAM1", "sire1", "0", "0", "1", "-9"),
        ("FAM1", "sire2", "0", "0", "1", "-9"),
        ("FAM1", "dam1",  "0", "0", "2", "-9"),
        ("FAM1", "dam2",  "0", "0", "2", "-9"),
    ]
    offspring = []
    for i in range(1, 9):
        sire = "sire1" if i <= 4 else "sire2"
        dam  = "dam1"  if i <= 4 else "dam2"
        offspring.append(("FAM1", f"off{i:02d}", sire, dam, "1", "-9"))
    return founders + offspring


def _build_geno_matrix(indivs, n_snps, rng):
    """
    Return geno[snp_idx][sample_idx] as integers in {0,1,2,9}.

    Founders are fully genotyped; offspring have _MISSING_RATE missingness.
    """
    n_founders = sum(1 for ind in indivs if ind[2] == "0")  # PAT == "0"
    genos = []
    for _ in range(n_snps):
        row = []
        for s_idx in range(len(indivs)):
            g = rng.choice([0, 0, 1, 1, 2])          # slight ref bias
            if s_idx >= n_founders and rng.random() < _MISSING_RATE:
                g = 9
            row.append(g)
        genos.append(row)
    return genos


def _write_bed(path, genos, n_samples):
    """Write PLINK .bed (SNP-major mode)."""
    with open(path, "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))            # magic bytes
        for row in genos:
            n_bytes = (n_samples + 3) // 4
            for byte_idx in range(n_bytes):
                byte_val = 0
                for bit_pos in range(4):
                    s_idx = byte_idx * 4 + bit_pos
                    g = row[s_idx] if s_idx < n_samples else 0
                    byte_val |= _BED_ENC.get(g, 0b01) << (bit_pos * 2)
                fh.write(bytes([byte_val]))


def _write_bim(path, n_chroms, n_snps_per_chrom):
    """Write PLINK .bim; A1 = ALT ('A'), A2 = REF ('G')."""
    records = []
    with open(path, "w") as fh:
        for chrom in range(1, n_chroms + 1):
            for i in range(1, n_snps_per_chrom + 1):
                snp_id = f"snp_chr{chrom}_{i:03d}"
                pos = i * 10_000
                fh.write(f"{chrom}\t{snp_id}\t0\t{pos}\tA\tG\n")
                records.append((str(chrom), snp_id, pos, "A", "G"))
    return records


def _write_fam(path, indivs):
    with open(path, "w") as fh:
        for ind in indivs:
            fh.write(" ".join(ind) + "\n")


# ── Session-scoped fixtures ───────────────────────────────────────────────────

@pytest.fixture(scope="session")
def synth_plink(tmp_path_factory):
    """
    Tiny synthetic PLINK dataset with pedigree structure.

    Returns a dict with keys:
      bed, bim, fam       — file paths (str)
      bfile               — common prefix (str)
      genos               — raw genotype matrix [snp][sample], values in {0,1,2,9}
      snp_records         — list of (chrom, snp_id, pos, a1, a2)
      indivs              — list of (FID, IID, PAT, MAT, SEX, PHENO) tuples
    """
    tmp = tmp_path_factory.mktemp("synth_plink")
    rng = random.Random(_SEED)

    indivs = _build_pedigree()
    n_samples = len(indivs)
    n_snps = _N_SNPS_PER_CHROM * _N_CHROMS

    genos = _build_geno_matrix(indivs, n_snps, rng)

    snp_records = _write_bim(tmp / "synth.bim", _N_CHROMS, _N_SNPS_PER_CHROM)
    _write_fam(tmp / "synth.fam", indivs)
    _write_bed(tmp / "synth.bed", genos, n_samples)

    return {
        "bed": str(tmp / "synth.bed"),
        "bim": str(tmp / "synth.bim"),
        "fam": str(tmp / "synth.fam"),
        "bfile": str(tmp / "synth"),
        "genos": genos,
        "snp_records": snp_records,
        "indivs": indivs,
    }


@pytest.fixture(scope="session")
def synth_ai2(synth_plink, tmp_path_factory):
    """
    AlphaImpute2-format genotype and pedigree files derived from synth_plink.

    Genotypes are drawn from the synthetic matrix (missing → '9').
    Returns a dict with keys: genotypes, pedigree (file paths).
    """
    tmp = tmp_path_factory.mktemp("synth_ai2")
    genos = synth_plink["genos"]          # [snp][sample]
    indivs = synth_plink["indivs"]

    n_samples = len(indivs)
    n_snps = len(genos)

    # Transpose: individual × SNP
    geno_by_indiv = [
        [genos[snp_idx][s_idx] for snp_idx in range(n_snps)]
        for s_idx in range(n_samples)
    ]

    geno_path = tmp / "genotypes.txt"
    with open(geno_path, "w") as fh:
        for ind, row in zip(indivs, geno_by_indiv):
            iid = ind[1]
            fh.write(iid + " " + " ".join(str(g) for g in row) + "\n")

    ped_path = tmp / "pedigree.txt"
    with open(ped_path, "w") as fh:
        for ind in indivs:
            _fid, iid, pat, mat, *_ = ind
            fh.write(f"{iid} {pat} {mat}\n")

    return {
        "genotypes": str(geno_path),
        "pedigree": str(ped_path),
        "geno_by_indiv": geno_by_indiv,
        "indivs": indivs,
    }


# ── Two-chip fixture ─────────────────────────────────────────────────────────
#
# synth_plink cannot express any of the failure modes the cross-array path
# actually hits: it is one chip, has no duplicate physical positions, no animal
# whose two typings disagree, and a .raw that encodes the same wrong assumption
# about which allele --export A counts that the production code once had.
# Fixture and code agreed with each other and both disagreed with reality,
# which is how a green suite coexisted with a module that could not complete a
# single real run.
#
# This fixture is built around those failures instead:
#   * HD carries markers LD does not, and LD carries markers HD does not
#   * one pair of HD markers sits at the same physical position
#   * some markers have A1/A2 swapped between the two chips
#   * two animals are mispaired -- their LD genotypes belong to another fish

_TC_N_ANIMALS = 20
_TC_N_SHARED = 24          # markers on both chips, per chromosome pair below
_TC_N_HD_ONLY = 16         # markers only the high-density chip carries
_TC_N_LD_ONLY = 3          # markers only the low-density chip carries
_TC_SWAPPED = {2, 7, 15}   # indexes into the shared markers, A1/A2 reversed on LD
_TC_MISPAIRED = ("off03", "off09")
_TC_SEED = 20260818


def _tc_animals():
    founders = [
        ("FAM1", "sire1", "0", "0", "1", "-9"),
        ("FAM1", "sire2", "0", "0", "1", "-9"),
        ("FAM1", "dam1", "0", "0", "2", "-9"),
        ("FAM1", "dam2", "0", "0", "2", "-9"),
    ]
    offspring = [
        (
            "FAM1",
            f"off{i:02d}",
            "sire1" if i <= 8 else "sire2",
            "dam1" if i <= 8 else "dam2",
            "1" if i % 2 else "2",
            "-9",
        )
        for i in range(1, _TC_N_ANIMALS - 3)
    ]
    return founders + offspring


def _tc_markers():
    """Return (shared, hd_only, ld_only) marker records as (chrom, id, pos, a1, a2)."""
    shared, hd_only, ld_only = [], [], []
    pos = 10_000
    for i in range(_TC_N_SHARED):
        chrom = "1" if i < _TC_N_SHARED // 2 else "2"
        shared.append((chrom, f"shared_{i:03d}", pos + i * 1_000, "A", "G"))
    for i in range(_TC_N_HD_ONLY):
        chrom = "1" if i < _TC_N_HD_ONLY // 2 else "2"
        hd_only.append((chrom, f"hdonly_{i:03d}", 500_000 + i * 1_000, "C", "T"))
    for i in range(_TC_N_LD_ONLY):
        ld_only.append(("1", f"ldonly_{i:03d}", 900_000 + i * 1_000, "A", "T"))
    # One duplicate physical position on the HD chip. FImpute aborts the whole
    # run on these ("SNPs with the same physical position found"), so the code
    # that drops them has to see one.
    hd_only.append(("1", "hdonly_dup", hd_only[0][2], "C", "T"))
    return shared, hd_only, ld_only


@pytest.fixture(scope="session")
def synth_two_chip(tmp_path_factory):
    """
    A low-density and a high-density PLINK fileset over the same animals.

    Returns a dict with the two prefixes plus the facts a test needs to assert
    against: which markers are shared, which are chip-only, which had their
    alleles swapped, and which animals are deliberately not the same fish on
    both chips.
    """
    tmp = tmp_path_factory.mktemp("synth_two_chip")
    rng = random.Random(_TC_SEED)

    animals = _tc_animals()
    n_samples = len(animals)
    shared, hd_only, ld_only = _tc_markers()

    # Truth genotypes, indexed by marker ID so both chips can draw from them.
    truth = {
        record[1]: [rng.choice([0, 0, 1, 1, 2]) for _ in range(n_samples)]
        for record in shared + hd_only + ld_only
    }

    mispaired_indexes = [
        idx for idx, animal in enumerate(animals) if animal[1] in _TC_MISPAIRED
    ]

    # LD chip: the shared markers plus its own, with a little between-array
    # noise, the swapped markers mirrored, and the mispaired animals carrying
    # genotypes drawn independently of the truth.
    ld_records = []
    ld_genos = []
    for index, record in enumerate(shared):
        chrom, snp_id, pos, a1, a2 = record
        swapped = index in _TC_SWAPPED
        ld_records.append((chrom, snp_id, pos, a2, a1) if swapped else record)
        row = []
        for sample_index, call in enumerate(truth[snp_id]):
            if sample_index in mispaired_indexes:
                call = rng.choice([0, 1, 2])
            elif rng.random() < 0.02:
                call = rng.choice([0, 1, 2])
            row.append(2 - call if swapped else call)
        ld_genos.append(row)
    for record in ld_only:
        ld_records.append(record)
        ld_genos.append(list(truth[record[1]]))

    # plink2 refuses a .bim whose chromosomes are interleaved ("has a split
    # chromosome"), so both chips are emitted in coordinate order -- which is
    # what a real export looks like anyway.
    ld_order = sorted(
        range(len(ld_records)), key=lambda i: (int(ld_records[i][0]), ld_records[i][2])
    )
    ld_records = [ld_records[i] for i in ld_order]
    ld_genos = [ld_genos[i] for i in ld_order]

    # HD chip: shared plus HD-only, in coordinate order, unmodified truth.
    hd_records = sorted(shared + hd_only, key=lambda r: (int(r[0]), r[2]))
    hd_genos = [list(truth[record[1]]) for record in hd_records]

    def _write(prefix, records, genos):
        with open(f"{prefix}.bim", "w") as handle:
            for chrom, snp_id, pos, a1, a2 in records:
                handle.write(f"{chrom}\t{snp_id}\t0\t{pos}\t{a1}\t{a2}\n")
        _write_fam(f"{prefix}.fam", animals)
        _write_bed(f"{prefix}.bed", genos, n_samples)

    ld_prefix = str(tmp / "ld")
    hd_prefix = str(tmp / "hd")
    _write(ld_prefix, ld_records, ld_genos)
    _write(hd_prefix, hd_records, hd_genos)

    return {
        "ld_bfile": ld_prefix,
        "hd_bfile": hd_prefix,
        "animals": animals,
        "shared_ids": [record[1] for record in shared],
        "hd_only_ids": [record[1] for record in hd_only],
        "ld_only_ids": [record[1] for record in ld_only],
        "swapped_ids": [shared[i][1] for i in sorted(_TC_SWAPPED)],
        "mispaired_iids": list(_TC_MISPAIRED),
        "duplicate_position_ids": ["hdonly_000", "hdonly_dup"],
        "truth": truth,
        "ld_records": ld_records,
        "hd_records": hd_records,
    }


@pytest.fixture(scope="session")
def synth_raw_counting_a2(synth_two_chip, tmp_path_factory):
    """
    A plink2 --export A .raw whose column suffix names a2, not a1.

    This is what our real exports produce, and assuming otherwise mirrors every
    genotype. Allelic r2 cannot see it -- squaring the correlation hides the
    sign -- so it surfaces only as a collapsed concordance, long after the run.
    """
    tmp = tmp_path_factory.mktemp("synth_raw_a2")
    records = synth_two_chip["hd_records"]
    animals = synth_two_chip["animals"]
    truth = synth_two_chip["truth"]

    path = tmp / "counts_a2.raw"
    header = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"] + [
        f"{snp_id}_{a2}" for _chrom, snp_id, _pos, _a1, a2 in records
    ]
    with open(path, "w") as handle:
        handle.write(" ".join(header) + "\n")
        for index, (fid, iid, pat, mat, sex, pheno) in enumerate(animals):
            # The .raw counts a2, so it reports the complement of the a1 dosage.
            counts = [str(2 - truth[record[1]][index]) for record in records]
            handle.write(" ".join([fid, iid, pat, mat, sex, pheno] + counts) + "\n")

    return {"raw": str(path), "records": records, "animals": animals}
