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
