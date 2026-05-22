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

`workflow_env` is **not** a pre-existing top-level conda env on most machines — it's
defined by `envs/workflow_env.yaml` and Snakemake builds it on-demand under
`.snakemake/conda/<hash>_` when running with `--use-conda`. To run pytest, either
create it once explicitly or use any env that already has pytest + bcftools + htslib + plink2 (e.g. an aquagen `base_env`).

```bash
# One-time setup (only if you want `conda activate workflow_env` to work)
conda env create -f envs/workflow_env.yaml

# Run tests
pytest                           # all tests
pytest -m "not slow"             # skip tests needing bgzip/tabix/plink2/snakemake
pytest tests/test_convert.py -v  # single file
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

### SLURM resource groups and partitions

Rules in the same group (`intersect`, `conform`, `beagle`) are submitted together. The `--group-components` flag controls how many chromosomes per job.

This cluster's partitions are split by node size, so every heavy rule must
declare `slurm_partition` in `resources:` — without it Snakemake submits to
the default `r6i-ondemand-large` (15 GiB / 2 CPU) and sbatch rejects the job
("CPU count per node can not be satisfied").

| Rule | `mem_mb` | `slurm_partition` |
|------|---------:|-------------------|
| `run_beagle` | 70000 | `r6i-ondemand-4xlarge` (124 GiB / 16 CPU) |
| `concat_chromosomes` | 64000 | `r6i-ondemand-4xlarge` |
| `bcftools_isec` / `conform_gt` | 32000 | `r6i-ondemand-2xlarge` (62 GiB / 8 CPU) |
| `merge_imputed_with_target_only` / `convert_ref_to_bref3` | 16000 | `r6i-ondemand-2xlarge` |
| `make_per_chrom_vcf` / `normalize_vcf` / `vcf_to_plink` | (none) | default `r6i-ondemand-large` |

When adding a new rule that needs >15 GiB RAM, pick the smallest partition
whose `RealMemory` (see `sinfo -o "%P %m %c"`) is ≥ the rule's `mem_mb`, and
add `slurm_partition = "<name>"` to its `resources:`.

## Conventions when editing rules

- **All rule outputs must live under `${output_dir}`.** Hard-coded paths
  (e.g. `"plink_binary/imputed_data.bed"`) cause parallel runs that share a
  cwd to clobber each other. Use `config["output_dir"] + "/..."` for every
  `output:` field.
- **plink2 invocations must always pass `--dog`.** The Snakefile prepends
  `--dog` to `config["plink_extra_flags"]` at parse time so every existing
  `params.extra_flags` / `params.extra` slot inherits it. When adding a new
  plink2 rule, route through `config["plink_extra_flags"]` (or pass the same
  rewritten value) rather than hardcoding flags — this guarantees non-human
  chromosome codes (1–38) work for every species we run (salmon, trout,
  livestock, dogs, humans). Do **not** add a separate `--chr-set N` flag.
- New per-chromosome rules should use `temp(...)` for intermediate VCFs so
  Snakemake cleans them up after `concat_chromosomes`.
- Heavy rules need an explicit `resources: mem_mb=...` so the SLURM executor
  asks for the right partition (`run_beagle`=70000, `concat_chromosomes`=64000,
  `bcftools_isec`/`conform_gt`=32000).

## Working with the GitHub remote

`gh` CLI is authenticated via HTTPS token. If `git fetch`/`push` fails with
`Permission denied (publickey)`, the origin is on SSH — switch it:

```bash
git remote set-url origin https://github.com/timknutsen/beagle_impute.git
```
