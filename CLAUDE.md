# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this pipeline does

Chromosome-wise genotype imputation pipeline. Converts PLINK binary files to per-chromosome VCFs, optionally harmonizes against a phased reference panel, imputes with **Beagle 5.5** or **AlphaImpute2**, and outputs a merged VCF + PLINK binary.

## Running the pipeline

```bash
# Local run
snakemake --use-conda --cores 8

# SLURM cluster (see snakemake_slurm_example.sh for a real example)
snakemake --use-conda --cores 48 --executor slurm --jobs 35 \
    --group-components intersect=5 conform=5 beagle=5

# Dry-run to validate DAG
snakemake --use-conda --cores 8 -n

# Override config values without editing config.yaml
snakemake --use-conda --cores 8 \
    --config bfile=/path/to/data beagle_jar=/path/to/beagle.jar
```

## Tests

```bash
# Run all tests (from repo root, with workflow_env active)
conda activate workflow_env
pytest

# Skip tests requiring external tools (bgzip, tabix, plink2, snakemake)
pytest -m "not slow"

# Run a single test file
pytest tests/test_convert.py -v
```

Test fixtures are generated synthetically at runtime in a temp dir — no binary files in git. Tests that require missing binaries auto-skip.

## Accuracy evaluation (separate Snakefile)

```bash
snakemake --snakefile Snakefile_accuracy --use-conda --cores 8
```

Reads `config.yaml` + `config_accuracy.yaml`. Set `accuracy_mode: mask_and_impute` or `cross_array`.

## Architecture

### Mode flags (set at parse time in `Snakefile`)

Three boolean flags control which rules and files are included:
- `_use_ref` — `reference_vcf` is set → includes `rules/intersect_and_conform.smk` and routes Beagle input through harmonization
- `_use_bref3` — `bref3_jar` is set → adds a bref3 conversion step before Beagle for faster reference loading
- `_use_alphaimpute2` — `imputer: "alphaimpute2"` → includes `rules/alphaimpute2.smk` instead of Beagle rules

### Main pipeline DAG (Beagle mode, no reference)

```
PLINK .bed → make_per_chrom_vcf → normalize_vcf → run_beagle → concat_chromosomes → vcf_to_plink
```

### With reference panel (`reference_vcf` set)

```
normalize_vcf → bcftools_isec → conform_gt → run_beagle → merge_imputed_with_target_only → concat_chromosomes
```

`bcftools_isec` outputs `0002.vcf.gz` (target ∩ ref, fed to Beagle) and `0000.vcf.gz` (target-only markers). After Beagle, `merge_imputed_with_target_only` adds the target-only markers back so no chip data is lost.

### Key config parameters to tune per species

| Parameter | Salmon/small livestock | Cattle/dogs | Humans |
|-----------|----------------------|-------------|--------|
| `ne` | ~500 | ~1000 | ~1,000,000 |
| `window` | 80 (default) | 80 | 40 |

### Rule files

- `Snakefile` — main entry point; all Beagle rules + concat + vcf_to_plink
- `rules/intersect_and_conform.smk` — `bcftools_isec`, `conform_gt`, `convert_ref_to_bref3` (only loaded when `reference_vcf` is set)
- `rules/alphaimpute2.smk` — AlphaImpute2 mode rules (only loaded when `imputer: "alphaimpute2"`)
- `rules/accuracy.smk` — imputation accuracy evaluation (only used via `Snakefile_accuracy`)
- `scripts/alphaimpute2_to_vcf.py` — converts AlphaImpute2 output to VCF
- `scripts/compute_accuracy_metrics.py` — concordance/r² metrics for accuracy evaluation

### Conda environments

- `envs/workflow_env.yaml` — main env: openjdk, bcftools, htslib, pandas, pytest
- `envs/alphaimpute2_env.yaml` — AlphaImpute2 env (Python 3.10 required; 3.12/3.14 incompatible)
- `envs/accuracy_env.yaml` — accuracy evaluation env

### SLURM resource groups

Rules in the same group (`intersect`, `conform`, `beagle`) are submitted together. The `--group-components` flag controls how many chromosomes per job. `run_beagle` uses 70 GB RAM; `conform_gt` and `bcftools_isec` use 32 GB.
