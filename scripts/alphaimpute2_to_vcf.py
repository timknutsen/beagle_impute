"""
Convert AlphaImpute2 imputed genotype file to bgzipped VCF.

AlphaImpute2 genotype format (input):
    individual_id g1 g2 ... gN   (space-separated; codes 0/1/2/9)

Output: bgzipped + tabix-indexed VCF using SNP metadata from the PLINK .bim file.

Allele conventions (consistent with plink2 --export A):
    0 = 0 copies of A1 (ALT) → REF/REF  → VCF 0/0
    1 = 1 copy  of A1        → REF/ALT  → VCF 0/1
    2 = 2 copies of A1       → ALT/ALT  → VCF 1/1
    9 = missing              → VCF ./.

Pure functions (read_bim, read_ai2_genotypes, convert) are importable so that
unit tests can exercise the conversion logic without a Snakemake context.
"""

import subprocess
from datetime import date

# Genotype code → VCF GT string
GT_MAP = {"0": "0/0", "1": "0/1", "2": "1/1", "9": "./."}


def read_bim(bim_path: str) -> list:
    """
    Parse a PLINK .bim file.

    Returns a list of (chrom, snp_id, pos, a1, a2) tuples where:
      a1 = ALT allele (PLINK column 5, counted by plink2 --export A)
      a2 = REF allele (PLINK column 6)
    """
    snps = []
    with open(bim_path) as fh:
        for line in fh:
            parts = line.strip().split()
            chrom, snp_id, _cm, pos, a1, a2 = parts
            snps.append((chrom, snp_id, int(pos), a1, a2))
    return snps


def read_ai2_genotypes(geno_path: str) -> tuple:
    """
    Parse an AlphaImpute2 genotype file.

    Returns (indiv_ids, rows) where rows[i] is a list of string genotype
    codes ("0"/"1"/"2"/"9") for individual i.
    """
    indivs, rows = [], []
    with open(geno_path) as fh:
        for line in fh:
            parts = line.strip().split()
            indivs.append(parts[0])
            rows.append(parts[1:])
    return indivs, rows


def convert(genotypes_path: str, bim_path: str, vcf_out_path: str) -> None:
    """
    Convert AlphaImpute2 imputed genotypes to a bgzipped, tabix-indexed VCF.

    Requires ``bgzip`` and ``tabix`` (from htslib) to be on PATH.
    The .tbi index is written alongside the output VCF.
    """
    snps = read_bim(bim_path)
    indivs, rows = read_ai2_genotypes(genotypes_path)

    n_snps = len(snps)
    if rows and len(rows[0]) != n_snps:
        raise ValueError(
            f"Genotype column count mismatch: BIM has {n_snps} SNPs "
            f"but imputed file has {len(rows[0])} columns"
        )

    bgzip_proc = subprocess.Popen(
        ["bgzip", "-c"],
        stdin=subprocess.PIPE,
        stdout=open(vcf_out_path, "wb"),
        stderr=subprocess.PIPE,
    )

    def emit(line: str) -> None:
        bgzip_proc.stdin.write((line + "\n").encode())

    emit("##fileformat=VCFv4.2")
    emit(f"##fileDate={date.today().strftime('%Y%m%d')}")
    emit("##source=AlphaImpute2")
    emit('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    emit("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(indivs))

    for snp_idx, (chrom, snp_id, pos, a1, a2) in enumerate(snps):
        # a2 = REF (PLINK A2), a1 = ALT (PLINK A1)
        gts = "\t".join(GT_MAP.get(row[snp_idx], "./.") for row in rows)
        emit(f"{chrom}\t{pos}\t{snp_id}\t{a2}\t{a1}\t.\tPASS\t.\tGT\t{gts}")

    bgzip_proc.stdin.close()
    rc = bgzip_proc.wait()
    if rc != 0:
        raise RuntimeError(
            f"bgzip exited with code {rc}: {bgzip_proc.stderr.read().decode()}"
        )

    subprocess.run(["tabix", "-f", "-p", "vcf", vcf_out_path], check=True)


# ---------------------------------------------------------------------------
# Snakemake entry point
# ---------------------------------------------------------------------------
if "snakemake" in dir():
    convert(
        snakemake.input.genotypes,
        snakemake.input.bim,
        snakemake.output.vcf,
    )
