"""
Compute imputation accuracy metrics by comparing an imputed VCF against a
truth VCF.

Metrics
-------
allelic_r2
    Squared Pearson correlation between imputed dosage (DS tag from Beagle,
    or GT converted to 0/1/2) and true genotype dosage (0/1/2).
    This is the standard imputation quality metric used in the literature.

concordance
    Fraction of samples where round(imputed_dosage) matches true_genotype.
    Sensitive to rare-allele errors because a heterozygote called as
    homozygous counts as a mistake regardless of MAF.

Outputs (all tab-separated)
---------------------------
metrics_by_snp.tsv
    Per-SNP: variant, maf, allelic_r2, concordance, n_evaluated
metrics_by_maf_bin.tsv
    Aggregated by MAF bin: n_snps, mean_r2, median_r2, mean_concordance
metrics_by_individual.tsv
    Per-animal: sample, concordance, n_evaluated
summary.tsv
    Overall: mean/median allelic R², % SNPs with R²≥0.8 / ≥0.9, mean concordance
"""

import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Compute imputation accuracy metrics (allelic R², concordance)."
    )
    p.add_argument("--imputed",   required=True, help="Imputed VCF (.vcf.gz, tabix-indexed)")
    p.add_argument("--truth",     required=True, help="Truth VCF (.vcf.gz, tabix-indexed)")
    p.add_argument("--out-dir",   required=True, help="Directory for output TSV files")
    p.add_argument(
        "--maf-bins", nargs="+", type=float,
        default=[0.01, 0.05, 0.1, 0.2, 0.5],
        help="MAF bin boundaries (default: 0.01 0.05 0.1 0.2 0.5)",
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# VCF helpers (thin wrappers around bcftools)
# ---------------------------------------------------------------------------

def get_samples(vcf: str) -> list:
    """Return list of sample IDs in a VCF."""
    r = subprocess.run(
        ["bcftools", "query", "-l", vcf],
        capture_output=True, text=True, check=True,
    )
    return r.stdout.strip().split("\n") if r.stdout.strip() else []


def get_variant_ids(vcf: str) -> list:
    """Return list of 'CHROM:POS:REF:ALT' strings for every variant."""
    r = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM:%POS:%REF:%ALT\n", vcf],
        capture_output=True, text=True, check=True,
    )
    return r.stdout.strip().split("\n") if r.stdout.strip() else []


def has_ds_tag(vcf: str) -> bool:
    """Return True if the VCF header declares a DS FORMAT field."""
    r = subprocess.run(
        ["bcftools", "view", "-h", vcf],
        capture_output=True, text=True, check=True,
    )
    return "##FORMAT=<ID=DS," in r.stdout


def _pipe_query(vcf: str, fmt: str, samples_file: str, regions_file: str) -> str:
    """
    Run: bcftools view -S samples -R regions vcf | bcftools query -f fmt
    Returns decoded stdout.
    """
    view = subprocess.Popen(
        ["bcftools", "view", "-S", samples_file, "-R", regions_file, vcf],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    query = subprocess.Popen(
        ["bcftools", "query", "-f", fmt],
        stdin=view.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    view.stdout.close()
    out, _ = query.communicate()
    view.wait()
    return out.decode()


# ---------------------------------------------------------------------------
# Matrix parsers
# ---------------------------------------------------------------------------

def _parse_gt(raw: str) -> np.ndarray:
    """
    Parse bcftools query '[%GT\\t]\\n' output into a float32 dosage matrix
    (n_variants × n_samples).  Missing calls (./. or .) become NaN.
    """
    rows = []
    for line in raw.strip().split("\n"):
        vals = line.rstrip("\t").split("\t")
        dosages = []
        for v in vals:
            v = v.replace("|", "/")
            alleles = v.split("/")
            try:
                dosages.append(float(sum(int(a) for a in alleles)))
            except ValueError:
                dosages.append(np.nan)
        rows.append(dosages)
    return np.array(rows, dtype=np.float32)


def _parse_ds(raw: str) -> np.ndarray:
    """
    Parse bcftools query '[%DS\\t]\\n' output into a float32 dosage matrix
    (n_variants × n_samples).  Missing ('.') becomes NaN.
    """
    rows = []
    for line in raw.strip().split("\n"):
        vals = line.rstrip("\t").split("\t")
        rows.append(
            [float(v) if v not in (".", "") else np.nan for v in vals]
        )
    return np.array(rows, dtype=np.float32)


def load_imputed_matrix(vcf: str, samples_file: str, regions_file: str) -> np.ndarray:
    """
    Load imputed dosages.  Prefers the DS tag (Beagle output) for dosage
    precision; falls back to GT-derived dosage (AlphaImpute2 or plain VCF).
    """
    if has_ds_tag(vcf):
        raw = _pipe_query(vcf, "[%DS\t]\n", samples_file, regions_file)
        return _parse_ds(raw)
    raw = _pipe_query(vcf, "[%GT\t]\n", samples_file, regions_file)
    return _parse_gt(raw)


def load_truth_matrix(vcf: str, samples_file: str, regions_file: str) -> np.ndarray:
    """Load truth genotypes as dosage matrix (GT only)."""
    raw = _pipe_query(vcf, "[%GT\t]\n", samples_file, regions_file)
    return _parse_gt(raw)


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def pearsonr_rowwise(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Compute per-row Pearson r between two (n_variants × n_samples) arrays,
    skipping NaN pairs.  Returns a 1-D array of length n_variants.
    """
    r = np.full(x.shape[0], np.nan, dtype=np.float64)
    for i in range(x.shape[0]):
        mask = ~(np.isnan(x[i]) | np.isnan(y[i]))
        if mask.sum() < 2:
            continue
        xi = x[i, mask].astype(np.float64)
        yi = y[i, mask].astype(np.float64)
        xi -= xi.mean()
        yi -= yi.mean()
        denom = np.sqrt((xi ** 2).sum() * (yi ** 2).sum())
        if denom > 0:
            r[i] = (xi * yi).sum() / denom
    return r


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # ── Find common samples and variants ─────────────────────────────────────
    imp_samples = get_samples(args.imputed)
    tru_samples = get_samples(args.truth)
    common_samples = [s for s in imp_samples if s in set(tru_samples)]

    if not common_samples:
        sys.exit(
            "ERROR: No common samples found between imputed and truth VCFs.\n"
            f"  Imputed samples ({len(imp_samples)}): {imp_samples[:5]} ...\n"
            f"  Truth   samples ({len(tru_samples)}): {tru_samples[:5]} ..."
        )

    imp_variants = get_variant_ids(args.imputed)
    tru_variants = get_variant_ids(args.truth)
    imp_set      = set(imp_variants)
    eval_variants = [v for v in tru_variants if v in imp_set]

    if not eval_variants:
        sys.exit(
            "ERROR: No common variants found between imputed and truth VCFs.\n"
            "This can happen if allele coding differs; check that the same\n"
            "reference genome and allele orientation were used."
        )

    print(f"Common samples     : {len(common_samples)}", flush=True)
    print(f"Truth variants     : {len(tru_variants)}", flush=True)
    print(f"Evaluated variants : {len(eval_variants)}", flush=True)

    # Write temp subsetting files
    tmp = tempfile.mkdtemp()
    samples_file = os.path.join(tmp, "samples.txt")
    regions_file = os.path.join(tmp, "regions.txt")

    with open(samples_file, "w") as fh:
        fh.write("\n".join(common_samples) + "\n")

    with open(regions_file, "w") as fh:
        for v in eval_variants:
            chrom, pos, *_ = v.split(":", 3)
            fh.write(f"{chrom}\t{pos}\n")

    # ── Load genotype matrices ────────────────────────────────────────────────
    print("Loading imputed data ...", flush=True)
    imp_matrix = load_imputed_matrix(args.imputed, samples_file, regions_file)

    print("Loading truth data ...", flush=True)
    tru_matrix = load_truth_matrix(args.truth, samples_file, regions_file)

    if imp_matrix.shape != tru_matrix.shape:
        sys.exit(
            f"ERROR: matrix shape mismatch — "
            f"imputed {imp_matrix.shape} vs truth {tru_matrix.shape}.\n"
            "Check that both VCFs were produced from the same set of animals."
        )

    # ── Per-SNP metrics ───────────────────────────────────────────────────────
    missing_mask = np.isnan(imp_matrix) | np.isnan(tru_matrix)
    match        = np.where(
        missing_mask, np.nan,
        (np.round(imp_matrix) == tru_matrix).astype(np.float32),
    )
    concordance = np.nanmean(match, axis=1)
    r2          = pearsonr_rowwise(imp_matrix, tru_matrix) ** 2

    # MAF from truth genotypes
    alt_counts = np.nansum(tru_matrix, axis=1)
    n_alleles  = np.sum(~np.isnan(tru_matrix), axis=1) * 2
    freq       = np.where(n_alleles > 0, alt_counts / n_alleles, np.nan)
    maf        = np.minimum(freq, 1.0 - freq)

    snp_df = pd.DataFrame({
        "variant"     : eval_variants,
        "maf"         : maf,
        "allelic_r2"  : r2,
        "concordance" : concordance,
        "n_evaluated" : np.sum(~missing_mask, axis=1),
    })
    snp_df.to_csv(
        os.path.join(args.out_dir, "metrics_by_snp.tsv"), sep="\t", index=False
    )

    # ── Per-individual metrics ────────────────────────────────────────────────
    indiv_miss    = np.isnan(imp_matrix) | np.isnan(tru_matrix)
    indiv_match   = np.where(
        indiv_miss, np.nan,
        (np.round(imp_matrix) == tru_matrix).astype(np.float32),
    )
    indiv_concord = np.nanmean(indiv_match, axis=0)

    indiv_df = pd.DataFrame({
        "sample"      : common_samples,
        "concordance" : indiv_concord,
        "n_evaluated" : np.sum(~indiv_miss, axis=0),
    })
    indiv_df.to_csv(
        os.path.join(args.out_dir, "metrics_by_individual.tsv"), sep="\t", index=False
    )

    # ── Per-MAF-bin metrics ───────────────────────────────────────────────────
    bins   = sorted(set([0.0] + list(args.maf_bins) + [0.5]))
    labels = [f"{bins[i]:.2f}-{bins[i+1]:.2f}" for i in range(len(bins) - 1)]
    snp_df["maf_bin"] = pd.cut(
        snp_df["maf"], bins=bins, labels=labels, include_lowest=True
    )
    maf_df = (
        snp_df.groupby("maf_bin", observed=True)
        .agg(
            n_snps           = ("allelic_r2",  "count"),
            mean_r2          = ("allelic_r2",  "mean"),
            median_r2        = ("allelic_r2",  "median"),
            mean_concordance = ("concordance", "mean"),
        )
        .reset_index()
    )
    maf_df.to_csv(
        os.path.join(args.out_dir, "metrics_by_maf_bin.tsv"), sep="\t", index=False
    )

    # ── Overall summary ───────────────────────────────────────────────────────
    valid_r2 = snp_df["allelic_r2"].dropna()
    summary = pd.DataFrame([{
        "n_samples"            : len(common_samples),
        "n_variants_truth"     : len(tru_variants),
        "n_variants_evaluated" : len(eval_variants),
        "mean_allelic_r2"      : valid_r2.mean(),
        "median_allelic_r2"    : valid_r2.median(),
        "pct_r2_ge_0.8"        : (valid_r2 >= 0.8).mean() * 100,
        "pct_r2_ge_0.9"        : (valid_r2 >= 0.9).mean() * 100,
        "mean_concordance"     : snp_df["concordance"].mean(),
    }])
    summary.to_csv(
        os.path.join(args.out_dir, "summary.tsv"), sep="\t", index=False
    )

    # Print summary to stdout (captured in Snakemake log)
    print("\n=== Imputation Accuracy Summary ===")
    print(f"Samples evaluated    : {len(common_samples)}")
    print(f"Variants evaluated   : {len(eval_variants)}")
    print(f"Mean allelic R\u00b2      : {valid_r2.mean():.4f}")
    print(f"Median allelic R\u00b2    : {valid_r2.median():.4f}")
    print(f"% SNPs with R\u00b2 \u2265 0.8 : {(valid_r2 >= 0.8).mean() * 100:.1f}%")
    print(f"% SNPs with R\u00b2 \u2265 0.9 : {(valid_r2 >= 0.9).mean() * 100:.1f}%")
    print(f"Mean concordance     : {snp_df['concordance'].mean():.4f}")
    print(f"\nMAF-bin breakdown:")
    print(maf_df.to_string(index=False))


if __name__ == "__main__":
    main()
