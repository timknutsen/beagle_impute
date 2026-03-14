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
"""

import subprocess
from datetime import date

GT_MAP = {"0": "0/0", "1": "0/1", "2": "1/1", "9": "./."}


# ---------------------------------------------------------------------------
# Read BIM: each row → (chrom, id, pos, a1=ALT, a2=REF)
# ---------------------------------------------------------------------------
snps = []
with open(snakemake.input.bim) as fh:
    for line in fh:
        parts = line.strip().split()
        chrom, snp_id, _cm, pos, a1, a2 = parts
        snps.append((chrom, snp_id, int(pos), a1, a2))

n_snps = len(snps)

# ---------------------------------------------------------------------------
# Read imputed genotypes (individual × SNP matrix)
# ---------------------------------------------------------------------------
indivs = []
rows = []
with open(snakemake.input.genotypes) as fh:
    for line in fh:
        parts = line.strip().split()
        indivs.append(parts[0])
        rows.append(parts[1:])

if rows and len(rows[0]) != n_snps:
    raise ValueError(
        f"Genotype column count mismatch: BIM has {n_snps} SNPs "
        f"but imputed file has {len(rows[0])} columns"
    )

# ---------------------------------------------------------------------------
# Write bgzipped VCF (pipe through bgzip to get bgzip format for tabix)
# ---------------------------------------------------------------------------
bgzip_proc = subprocess.Popen(
    ["bgzip", "-c"],
    stdin=subprocess.PIPE,
    stdout=open(snakemake.output.vcf, "wb"),
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
_rc = bgzip_proc.wait()
if _rc != 0:
    raise RuntimeError(f"bgzip exited with code {_rc}: {bgzip_proc.stderr.read().decode()}")

# Index with tabix
subprocess.run(
    ["tabix", "-f", "-p", "vcf", snakemake.output.vcf],
    check=True,
)
