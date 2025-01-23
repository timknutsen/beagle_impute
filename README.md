
# Chromosome-wise Imputation Pipeline

This pipeline converts PLINK files to VCF per chromosome, optionally harmonizes against a reference VCF, and imputes with Beagle.

## 1. Requirements

- **Conda** and **Snakemake** installed
- Everything else (PLINK, bcftools, Java, etc.) is automatically handled by Snakemake via `envs/workflow_env.yaml`.

## 2. Configuration

Inside the **Snakefile**, we have the following default `config`:

```python
localrules: make_per_chrom_vcf, normalize_vcf, concat_chromosomes, vcf_to_plink

config = {
    # Paths that typically need changing
    "bfile": "tests/data/test_salmon",
    "reference_vcf": "tests/data/test_salmon_ref.PHASED.vcf.gz",  # (Optional) path to reference VCF

    # Default paths and parameters
    "plink_path": "/mnt/efshome/home/timknu/bioinf_tools/plink2.0/plink2",
    "beagle_jar": "/mnt/efshome/home/timknu/bioinf_tools/beagle5.5/beagle.17Dec24.224.jar",
    "output_dir": "vcf_output",
    "beagle_params": {
        "window": 80,
        "overlap": 10,
        "ne": 500,
        "nthreads": 1
    },
    "conform_gt_jar": "/mnt/efshome/home/timknu/bioinf_tools/conform-gt.24May16.cee.jar"
}
```

### Key Parameters
- **`bfile`**: Prefix of your PLINK `.bed/.bim/.fam`
- **`reference_vcf`**: (Optional) Reference panel VCF for harmonization
- **`plink_path`**: Location of PLINK2 executable
- **`beagle_jar`**: Path to Beagle JAR file
- **`conform_gt_jar`**: Path to conform-gt JAR (only used if `reference_vcf` is set)

You can override these in your terminal:
```bash
snakemake --use-conda \
          --cores 8 \
          --config bfile=/path/to/mydata \
                   reference_vcf=/path/to/ref.vcf.gz \
                   plink_path=/path/to/plink2
```

## 3. Running the Pipeline

1. **Clone** this repository and ensure `Snakefile` + `envs/workflow_env.yaml` are present.
2. **Check** or **override** the `config` paths above.
3. **Run**:
   ```bash
   snakemake --use-conda --cores 8
   ```
   (Add `--executor slurm`, resource settings, or job grouping as needed.)

## 4. Outputs

- Per-chromosome imputed VCFs: `vcf_output/imputed/chr{chrom}.vcf.gz`
- Single combined VCF: `vcf_output/all_chromosomes.vcf.gz`
- Final PLINK files: `plink_binary/imputed_data.bed/.bim/.fam`

## 5. Notes

- **Harmonization** (bcftools + conform-gt) only runs if `reference_vcf` is specified.
- **Logs** are written to `logs/` per step and per chromosome.
- The `beagle_params` (e.g. `nthreads`) should match the resources allocated on your cluster or local machine.