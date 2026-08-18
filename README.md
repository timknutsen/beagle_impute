# Chromosome-wise Imputation Pipeline

This pipeline imputes chromosome-wise PLINK data with Beagle 5.5, FImpute3, or
AlphaImpute2. Beagle accepts a phased VCF reference; FImpute accepts a PLINK
reference panel; both engines also run without an external reference.

## 1. Requirements

- **Conda** and **Snakemake** installed
- **PLINK2**, plus the selected engine: Beagle 5.5 JAR, licensed FImpute3
  binary, or AlphaImpute2 (see paths in `config.yaml`)
- Everything else (bcftools, pandas, Java) is installed automatically via `envs/workflow_env.yaml`

## 2. Choosing an Imputer

The pipeline supports three imputation engines. Set `imputer` in `config.yaml`:

| Setting | Engine | Best for |
|---------|--------|----------|
| `imputer: "beagle"` (default) | Beagle 5.5 | LD-based imputation; no pedigree needed |
| `imputer: "alphaimpute2"` | AlphaImpute2 | Pedigree + population; structured livestock/aquaculture |
| `imputer: "fimpute"` | FImpute3 | Fast pedigree + population imputation; licensed binary required |

AlphaImpute2 is installed automatically via pip inside `envs/alphaimpute2_env.yaml` — no JAR file needed. It requires Python 3.10 (3.12 and 3.14 are incompatible).

## 3. Configuration

Edit `config.yaml` in the repository root before running. The key parameters to set:

```yaml
bfile: "path/to/your/data"          # PLINK prefix (without .bed/.bim/.fam)
reference_vcf: ""                    # Optional phased reference VCF; leave empty to skip harmonization
output_dir: "vcf_output"

plink_path: "plink2"                 # Path to PLINK2 executable
beagle_jar: "bin/beagle.jar"         # Auto-downloaded on first run (Beagle 5.5)
conform_gt_jar: ""                   # Only needed if reference_vcf is set
bref3_jar: ""                        # Optional: path to bref3.jar (same download page as Beagle)
                                     # Converts reference VCF to binary bref3 format once —
                                     # Beagle then loads it 3–43× faster per run

# Extra PLINK flags. `--dog` is ALWAYS prepended automatically so plink2
# accepts non-human chromosome codes (1–38: salmon, trout, livestock, human).
# Use this field only for additional dataset-specific flags (e.g. "--maf 0.01").
plink_extra_flags: ""

beagle_params:
  window: 80       # Window size in cM (Beagle default is 40; 80 gives more LD context)
  overlap: 20      # Overlap between windows in cM (~25% of window is recommended)
  ne: 500          # Effective population size — tune to your species:
                   #   ~500 for salmon/small livestock, ~1000 for cattle/dogs,
                   #   ~1,000,000 for large outbred human populations (Beagle default)
  nthreads: 4      # Beagle threads — match to your CPU/cluster allocation
```

For AlphaImpute2 mode, use `alphaimpute2_params`; for FImpute, use
`fimpute_params` and set `reference_bfile` to a PLINK prefix or leave it empty.
`reference_vcf` is Beagle-only, and `fimpute_params.reference_bfile` is
FImpute-only: stale settings for another engine are ignored.

You can also override individual values on the command line without editing the file:

```bash
snakemake --use-conda --cores 8 \
          --config bfile=/path/to/mydata \
                   beagle_jar=/path/to/beagle.17Dec24.jar \
                   plink_path=/path/to/plink2
```

## 4. Running the Pipeline

```bash
snakemake --use-conda --cores 8
```

For SLURM cluster execution, see `snakemake_slurm_example.sh` and add `--set-resources` flags for partition names as needed.

### Real-world example: Rainbow trout with reference panel

Impute a small trout cohort (105 fish, 47k SNPs) using a large phased reference panel (4,439 fish, 42k SNPs):

```bash
snakemake --use-conda --cores 8 \
    --config \
        bfile="/mnt/efshome/aquagen/projects/AG_global_breeding/Rainbowtrout/issue_191_trout_parent_control/test_output.filtered_miss_hwe" \
        reference_vcf="/mnt/efshome/aquagen/projects/AG_global_breeding/Rainbowtrout/latemat_trout_2023/6159_regnbueorret_nye_arter/all_batches_clean_good_probes_merged.PHASED.vcf.gz" \
        output_dir="issue_191_ref_imputed"
```

For a large cohort (25k fish, 55k SNPs) without a reference panel, run on SLURM:

```bash
snakemake --use-conda \
    --cores 48 \
    --executor slurm \
    --jobs 35 \
    --group-components beagle=5 \
    --config \
        bfile="/mnt/efshome/aquagen/projects/AG_global_breeding/Rainbowtrout/RT24/plink/RT24_fish" \
        output_dir="RT24_imputed"
```

## 5. Outputs

- Per-chromosome imputed VCFs: `vcf_output/imputed/chr{chrom}.vcf.gz`
- Single combined VCF + index: `vcf_output/all_chromosomes.vcf.gz` (+ `.tbi`)
- Final PLINK files: `plink_binary/imputed_data.bed/.bim/.fam`
- Logs per step and chromosome: `logs/`

**AlphaImpute2 mode outputs** (`imputer: "alphaimpute2"`):
- Genome-wide imputed VCF + index: `vcf_output/alphaimpute2/all_chromosomes.vcf.gz` (+ `.tbi`)
- Final PLINK files: `plink_binary/imputed_data.bed/.bim/.fam`
- Intermediate AlphaImpute2 files: `vcf_output/alphaimpute2_input/` and `vcf_output/alphaimpute2_output/`
- Logs: `logs/`

**FImpute mode outputs** (`imputer: "fimpute"`):
- Genome-wide imputed VCF + index: `vcf_output/fimpute/all_chromosomes.vcf.gz` (+ `.tbi`)
- Final PLINK files: `vcf_output/plink_binary/imputed_data.bed/.bim/.fam`
- With `fimpute_params.reference_bfile`, the final files contain target
  (`bfile`) animals only; reference-chip animals are not output samples.

## 6. Accuracy Evaluation

Run the accuracy workflow separately from the main imputation workflow:

```bash
snakemake --snakefile Snakefile_accuracy --use-conda --cores 8
```

`config_accuracy.yaml` supports three modes:

| Mode | Purpose |
|------|---------|
| `mask_and_impute` | Hold out one validation set, mask it to LD density, impute, and compare to truth |
| `cross_array` | Compare a real low-density array against a real high-density array for the same animals, scored K-fold and identity-gated |
| `kfold_mask_and_impute` | Run K-fold animal CV from a shared LD marker panel to full density |

For the 10k-to-full benchmark requested for real data, use:

```yaml
accuracy_mode: "kfold_mask_and_impute"
accuracy_output_dir: "accuracy_cv"

cv:
  n_folds: 10
  target_n_snps: 10000
  target_snp_list: ""
  random_seed: 42
  imputers: ["beagle", "alphaimpute2", "fimpute"]

fimpute_params:
  executable: "/mnt/efshome/applications/FImpute3/2026/FImpute3"
  nthreads: 1
```

In this mode each fold masks the held-out animals to the same 10k SNP panel.
The other nine folds remain full density as the reference set.  Beagle uses
the held-out LD animals as `gt=` and the non-held-out full-density animals as
`ref=`; AlphaImpute2 and FImpute use a combined dataset where only the held-out
animals are masked outside the LD panel.

**Accuracy is measured on imputed markers only.** The truth set excludes the LD
panel, since those markers were given to the imputer as observed genotypes and
are returned unchanged — including them would inflate the reported accuracy by
roughly `n_panel / n_total` (on a 50K chip with a 10k panel, about a fifth of
all markers). The same exclusion applies to `mask_and_impute` and, for the
overlap between the two arrays, to `cross_array`.

### `cross_array`

The same animals on two arrays, scored K-fold: each fold contributes its real
low-density typings, and the other folds' high-density genotypes are the
reference panel. A held-out fold rather than an external panel is what makes
the mode work at all when nearly every high-density animal is already in the
cohort — on the salmon arrays exactly one Ssa70kv4 fish exists outside it.

The LD panel is not sampled here: it is exactly the markers both arrays carry.

```bash
snakemake --snakefile Snakefile_accuracy --use-conda --cores 8 \
  --config accuracy_mode=cross_array \
           accuracy_output_dir=accuracy_v3_to_v4 \
           cross_array_ld_bfile=/path/to/Ssa70kv3 \
           cross_array_hd_bfile=/path/to/Ssa70kv4 \
           cv_n_folds=5 \
           cv_imputers="beagle fimpute"
```

**An identity gate runs before anything is scored, and is not optional.**
Pairing an animal's two typings on its ID assumes both came from one fish; when
an ID covers two physical samples the pair is two unrelated animals and the run
scores an imputation that never had a chance. Concordance between the arrays at
their shared markers separates the cases cleanly — the same animal lands at
0.95–1.0, two different fish near 0.5 — so `cross_array.identity_threshold`
(default 0.90) sits in an empty gap. Failures land in
`setup/identity_fail.ids`; the run aborts entirely below
`cross_array.min_identity_pass_rate`.

To mask a high-density cohort down to a real lower-density array instead of a
random thinning — e.g. V4 fish downsampled to the 50K markers — build the panel
first and hand it to `kfold_mask_and_impute`:

```bash
python scripts/shared_marker_list.py \
  --bim /path/to/Ssa50kv1.bim /path/to/Ssa70kv4.bim \
  --out 50k_v4_shared.txt

snakemake --snakefile Snakefile_accuracy --use-conda --cores 8 \
  --config accuracy_mode=kfold_mask_and_impute \
           bfile=/path/to/Ssa70kv4 \
           cv_target_snp_list=50k_v4_shared.txt
```

### Which markers impute reliably

Both CV modes write `snp_reliability.tsv` and `reliable_markers.txt`. A marker
qualifies on its **worst** fold, not its mean — one collapsed fold out of five
is exactly what a mean hides. Passing `--root` more than once intersects runs,
so a marker earns its place only by holding up in every test it was part of:

```bash
python scripts/aggregate_snp_reliability.py \
  --root accuracy_v3_to_v4 accuracy_50k_masked_v4 \
  --imputers beagle fimpute --folds 1 2 3 4 5 \
  --bim /path/to/Ssa70kv4.bim \
  --out v4_snp_reliability.tsv --reliable-out v4_reliable_markers.txt
```

### QC is expected upstream

The pipeline does not filter genotypes. Bfiles handed to any accuracy mode
should already be QC'd **per chip, before merging or pairing** — genotype status
`OK`, `--geno`, `--maf`, and `--hwe <p> <k> midp keep-fewhet`. Run that QC
before writing the pedigree into `.fam` (plink2's `--hwe` considers founders
only), and apply sample-missingness exclusions to both arrays so the pairing
survives.

Command-line overrides use flat names because Snakemake does not accept dotted
keys in `--config`:

```bash
snakemake --snakefile Snakefile_accuracy --use-conda --cores 8 \
  --config accuracy_mode=kfold_mask_and_impute \
           accuracy_output_dir=accuracy_cv \
           bfile=/path/to/full_density_plink \
           cv_n_folds=10 \
           cv_target_n_snps=10000 \
           cv_imputers="beagle alphaimpute2 fimpute"
```

CV outputs:

- `accuracy_cv/cv_summary.tsv`: long table with one row per imputer/fold/metric
- `accuracy_cv/cv_imputer_summary.tsv`: mean, SD, and count by imputer/metric
- `accuracy_cv/{imputer}/fold{n}/summary.tsv`: per-fold summary
- `accuracy_cv/{imputer}/fold{n}/metrics_by_snp.tsv`
- `accuracy_cv/{imputer}/fold{n}/metrics_by_maf_bin.tsv`
- `accuracy_cv/{imputer}/fold{n}/metrics_by_individual.tsv`
- `accuracy_cv/snp_reliability.tsv` and `accuracy_cv/reliable_markers.txt`
- `cross_array` also writes `setup/identity.tsv` and `setup/identity_fail.ids`

## 6b. Building a reference panel

Beagle's `reference_vcf` must be phased, and until now nothing in this repo
produced one. `Snakefile_refpanel` does:

```bash
snakemake --snakefile Snakefile_refpanel --use-conda --executor slurm --jobs 30 \
  --config refpanel_bfile=/path/to/qc_Ssa70kv1 \
           refpanel_name=Ssa70kv1 \
           refpanel_pedigree=/path/to/pedigree.tsv \
           refpanel_max_per_family=4 \
           output_dir=refpanel_output
```

Per chromosome: export → `bcftools norm -d snps` → **phase with Beagle** (`gt=`,
no `ref=`, which also fills missing calls) → **bref3**. Outputs per variant:

| File | What |
|---|---|
| `phased/chr<N>.vcf.gz` | phased, no missing genotypes |
| `bref3/chr<N>.bref3` | what Beagle should actually read — smaller and far faster to load |
| `manifest.tsv` | every animal, kept or dropped, with the reason |
| `panel_report.tsv` | source bfile, pedigree, cap, seed, exclusions, counts |

It is a one-time cost of roughly 25–30 node-minutes genome-wide; every later
imputation loads the bref3 instead of re-phasing. The array is a parameter, so
the same rules build the v1, v3 and v4 panels.

**Diversity.** `refpanel_max_per_family` caps animals per sire × dam pair,
because a breeding cohort is not a diverse sample — full sibs share long
haplotype blocks, so extra sibs cost phasing time without adding much a panel
can use. It needs a pedigree; an array export has `0` in both parent columns, so
resolve parents from GPA first. Founders are each treated as their own family
rather than capped as one group.

**Holding animals out.** `refpanel.variants` maps a name to ID files to exclude,
so one selection pass produces both the production panel and benchmark-safe
variants. Exclusions are applied *before* phasing on purpose: subsetting an
already-phased panel is cheaper but leaks, since the excluded animals' genotypes
would already have shaped the phase of the relatives that remain.

Then use it:

```bash
snakemake --use-conda --cores 8 \
  --config bfile=/path/to/target \
           reference_vcf=refpanel_output/full/phased/chr{chrom}.vcf.gz \
           bref3_jar=bin/bref3.jar
```

## 7. Testing

The test suite lives in `tests/` and uses [pytest](https://docs.pytest.org). No pipeline tools (plink2, Beagle, AlphaImpute2) are required for the default test run.

### Running the tests

```bash
# Option 1: create the workflow_env once (matches envs/workflow_env.yaml) and use it.
# Note: Snakemake's --use-conda builds this env under .snakemake/conda/ on demand,
# so you only need to create the top-level env if you want to `conda activate` it.
conda env create -f envs/workflow_env.yaml
conda activate workflow_env
pytest

# Option 2: install pytest into any existing env that already has bcftools/htslib/plink2.
pip install pytest && pytest
```

### What gets tested

| Test file | What it covers | External tools needed |
|-----------|---------------|----------------------|
| `test_convert.py` (pure-Python group) | `GT_MAP` encoding, `read_bim`, `read_ai2_genotypes`, pedigree column format, genotype roundtrip against fixture matrix | None |
| `test_convert.py` (`TestConvertRoundtrip`) | Full `convert()` call: bgzipped VCF content, REF/ALT allele direction, GT strings, tabix index creation, mismatch error | `bgzip`, `tabix` (htslib) |
| `test_dryrun.py` | Snakemake DAG validation for Beagle and FImpute, each with/without a reference, plus AlphaImpute2; also checks engine isolation | `snakemake` |
| `test_accuracy_cv.py` | Ten-fold CV setup, FImpute I/O helpers, CV aggregation, and CV Snakemake dry-run | `snakemake` for dry-run only |
| `test_cross_array.py` | Cross-array pairing, the shared marker panel, the identity gate, allele orientation, the LD splice, marker reliability, and the `cross_array` dry-run | `snakemake` for dry-run only |

Tests that require missing tools are **automatically skipped** with a clear reason — they never fail due to a missing binary.

### Synthetic test dataset

The test fixtures in `tests/conftest.py` generate a small synthetic PLINK binary **at test time in a temp directory** — no binary files are stored in git. The dataset is:

- **12 individuals**: 4 fully-genotyped founders (`sire1`, `sire2`, `dam1`, `dam2`) + 8 offspring with known parents and 25% missingness
- **80 SNPs** across 2 chromosomes (40 each), positions spaced 10 kbp apart
- **Deterministic** (fixed random seed), so failures are reproducible

This pedigree structure allows meaningful tests of both the genotype encoding logic and the pedigree file format required by AlphaImpute2.

## 8. Notes

### AlphaImpute2 key parameters

- **`cycles`**: number of peeling rounds for pedigree imputation (default 4). Increasing to 6–8 can improve accuracy in deep pedigrees at the cost of runtime.
- **`length`**: effective chromosome length in Morgans (default 1.0 = 100 cM). Adjust for species with unusually short or long chromosomes.
- **`final_peeling_threshold`**: confidence cutoff for calling a genotype (0.1 = best-guess, 0.9+ = high-confidence only). Use higher values when downstream analyses require reliable homozygosity calls.
- **`hd_threshold`**: proportion of non-missing genotypes required to use an individual as a high-density reference (default 0.95).
- **`ped_only` / `pop_only`**: force pedigree-only or population-only imputation. Default uses both jointly.

AlphaImpute2 processes all chromosomes in a single run (unlike Beagle, which is run per chromosome). The pipeline therefore produces a single genome-wide VCF rather than per-chromosome files in this mode.

### Reference panel behaviour

For Beagle, `reference_vcf` must be phased and complete. The pipeline runs:
1. `bcftools isec` — splits target markers into two groups: those **in the reference** and those **only in the target**
2. `conform-gt` — harmonises strand/allele order of the intersection to match the reference
3. `Beagle` with `ref=` — phases and imputes only the intersection markers (Beagle drops target-only markers by design)
4. **Merge** — target-only markers are merged back with the Beagle output so the final VCF contains all original chip markers plus any markers imputed from the reference

This means you can use a reference panel to improve phasing and imputation of overlapping markers without losing any of your original chip data.

For FImpute, `fimpute_params.reference_bfile` activates a two-chip run:
reference is chip 1 and target `bfile` is chip 2. The final VCF/PLINK contains
only target animals. Leave the setting empty for one-chip phasing and sporadic
missing-call fill.

### bref3: faster reference loading
Setting `bref3_jar` causes the pipeline to convert the reference VCF to bref3 binary format before running Beagle. This conversion runs once; subsequent per-chromosome Beagle jobs load it 3–43× faster depending on reference panel size. The `bref3.jar` file is available on the same download page as the Beagle JAR.

### Other notes
- **`ne`** has a large effect on imputation accuracy for non-human populations. Use a value appropriate for your species — the Beagle default of 1,000,000 is designed for large outbred human populations and performs poorly for livestock and aquaculture species.
- **`nthreads`** should match the number of CPUs you allocate. Beagle scales well with more threads — the default of 4 is a reasonable starting point; increase for cluster jobs.
- **`plink_extra_flags`**: `--dog` is always prepended automatically so non-human chromosome codes (1–38) work for every species we run (salmon, trout, livestock, dogs, humans). Use this field only for *additional* dataset-specific flags (e.g. `"--maf 0.01"`). Empty by default.

## Project Structure

- **Snakefile**: Main imputation workflow (Beagle/FImpute/AlphaImpute2)
- **Snakefile_accuracy**: Imputation accuracy evaluation workflow
- **config.yaml**: Main pipeline configuration
- **config_accuracy.yaml**: Accuracy evaluation configuration
- **rules/**: Modular Snakemake rule files (per-imputer, reference panel, accuracy)
- **scripts/**: Python utilities for VCF conversion and accuracy metrics
- **envs/**: Conda environment YAMLs for reproducible setups
- **tests/**: Pytest-based test suite with synthetic data fixtures

## Workflows

### 1. Imputation Pipeline
- Converts PLINK data and imputes with Beagle, FImpute, or AlphaImpute2
- Select imputer in `config.yaml` (`beagle`, `fimpute`, or `alphaimpute2`)
- Modular rules: per-chromosome, reference panel, and imputer-specific logic

### 2. Accuracy Evaluation
- Run with `snakemake --snakefile Snakefile_accuracy`
- Configure in `config_accuracy.yaml`
- Modes: `mask_and_impute`, `cross_array`, and `kfold_mask_and_impute`
- Produces detailed accuracy metrics and summaries

## Scripts
- `scripts/alphaimpute2_to_vcf.py`: Converts AlphaImpute2 output to VCF
- `scripts/compute_accuracy_metrics.py`: Computes R², concordance, and summary stats
- `scripts/make_accuracy_cv_setup.py`: Creates deterministic CV folds and LD SNP panels
- `scripts/fimpute_io.py`: Converts PLINK raw exports to FImpute input files and FImpute output to VCF
- `scripts/aggregate_cv_metrics.py`: Aggregates per-fold CV summary metrics
- `scripts/select_refpanel.py`: Reference panel membership — family capping, exclusions, and the manifest
- `scripts/aggregate_snp_reliability.py`: Per-marker R² across folds and runs; writes the reliable-marker list
- `scripts/shared_marker_list.py`: Intersects two or more `.bim` files into a marker list
- `scripts/pair_animals.py`: Matches animals across two arrays on individual ID
- `scripts/identity_gate.py`: Turns LD-vs-HD concordance into the animal list a cross-array run uses
- `scripts/mask_validation_genotypes.py`: Masks genotypes outside an LD panel; `--replace-from` splices a second array's calls in

## Environments
- `envs/workflow_env.yaml`: Main pipeline dependencies
- `envs/alphaimpute2_env.yaml`: AlphaImpute2-specific (Python 3.10)
- `envs/accuracy_env.yaml`: Accuracy evaluation tools
