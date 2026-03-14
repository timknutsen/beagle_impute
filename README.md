
# Chromosome-wise Imputation Pipeline

This pipeline converts PLINK files to VCF per chromosome, optionally harmonizes against a reference VCF, and imputes with Beagle 5.5.

## 1. Requirements

- **Conda** and **Snakemake** installed
- **PLINK2** and **Beagle 5.5 JAR** (see paths in `config.yaml`)
- Everything else (bcftools, pandas, Java) is installed automatically via `envs/workflow_env.yaml`

## 2. Choosing an Imputer

The pipeline supports two imputation engines. Set `imputer` in `config.yaml`:

| Setting | Engine | Best for |
|---------|--------|----------|
| `imputer: "beagle"` (default) | Beagle 5.5 | LD-based imputation; no pedigree needed |
| `imputer: "alphaimpute2"` | AlphaImpute2 | Pedigree + population; structured livestock/aquaculture |

AlphaImpute2 is installed automatically via pip inside `envs/alphaimpute2_env.yaml` — no JAR file needed. It requires Python 3.10 (3.12 and 3.14 are incompatible).

## 3. Configuration

Edit `config.yaml` in the repository root before running. The key parameters to set:

```yaml
bfile: "path/to/your/data"          # PLINK prefix (without .bed/.bim/.fam)
reference_vcf: ""                    # Optional phased reference VCF; leave empty to skip harmonization
output_dir: "vcf_output"

plink_path: "plink2"                 # Path to PLINK2 executable
beagle_jar: "/path/to/beagle.17Dec24.jar"   # Latest Beagle 5.5 JAR
conform_gt_jar: ""                   # Only needed if reference_vcf is set
bref3_jar: ""                        # Optional: path to bref3.jar (same download page as Beagle)
                                     # Converts reference VCF to binary bref3 format once —
                                     # Beagle then loads it 3–43× faster per run

# Species-specific PLINK flags (default: none)
# Examples:
#   "--dog --aec"   for canine data with non-standard chromosome names
#   ""              for human or most other species
plink_extra_flags: ""

beagle_params:
  window: 80       # Window size in cM (Beagle default is 40; 80 gives more LD context)
  overlap: 20      # Overlap between windows in cM (~25% of window is recommended)
  ne: 500          # Effective population size — tune to your species:
                   #   ~500 for salmon/small livestock, ~1000 for cattle/dogs,
                   #   ~1,000,000 for large outbred human populations (Beagle default)
  nthreads: 4      # Beagle threads — match to your CPU/cluster allocation
```

For AlphaImpute2 mode, replace `beagle_jar` with the `alphaimpute2_params` block (see `config.yaml` comments) and set `imputer: "alphaimpute2"`. The `beagle_jar`, `conform_gt_jar`, and `bref3_jar` fields are ignored in this mode.

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

## 6. Testing

The test suite lives in `tests/` and uses [pytest](https://docs.pytest.org). No pipeline tools (plink2, Beagle, AlphaImpute2) are required for the default test run.

### Running the tests

```bash
# Activate the workflow conda env (has pytest) and run from the repo root:
conda activate workflow_env
pytest

# Or install pytest standalone:
pip install pytest && pytest
```

### What gets tested

| Test file | What it covers | External tools needed |
|-----------|---------------|----------------------|
| `test_convert.py` (pure-Python group) | `GT_MAP` encoding, `read_bim`, `read_ai2_genotypes`, pedigree column format, genotype roundtrip against fixture matrix | None |
| `test_convert.py` (`TestConvertRoundtrip`) | Full `convert()` call: bgzipped VCF content, REF/ALT allele direction, GT strings, tabix index creation, mismatch error | `bgzip`, `tabix` (htslib) |
| `test_dryrun.py` | Snakemake DAG validation for Beagle (no ref), Beagle (with ref), and AlphaImpute2 modes — confirms rules resolve without executing anything | `snakemake` |

Tests that require missing tools are **automatically skipped** with a clear reason — they never fail due to a missing binary.

### Synthetic test dataset

The test fixtures in `tests/conftest.py` generate a small synthetic PLINK binary **at test time in a temp directory** — no binary files are stored in git. The dataset is:

- **12 individuals**: 4 fully-genotyped founders (`sire1`, `sire2`, `dam1`, `dam2`) + 8 offspring with known parents and 25% missingness
- **80 SNPs** across 2 chromosomes (40 each), positions spaced 10 kbp apart
- **Deterministic** (fixed random seed), so failures are reproducible

This pedigree structure allows meaningful tests of both the genotype encoding logic and the pedigree file format required by AlphaImpute2.

## 7. Notes

### AlphaImpute2 key parameters

- **`cycles`**: number of peeling rounds for pedigree imputation (default 4). Increasing to 6–8 can improve accuracy in deep pedigrees at the cost of runtime.
- **`length`**: effective chromosome length in Morgans (default 1.0 = 100 cM). Adjust for species with unusually short or long chromosomes.
- **`final_peeling_threshold`**: confidence cutoff for calling a genotype (0.1 = best-guess, 0.9+ = high-confidence only). Use higher values when downstream analyses require reliable homozygosity calls.
- **`hd_threshold`**: proportion of non-missing genotypes required to use an individual as a high-density reference (default 0.95).
- **`ped_only` / `pop_only`**: force pedigree-only or population-only imputation. Default uses both jointly.

AlphaImpute2 processes all chromosomes in a single run (unlike Beagle, which is run per chromosome). The pipeline therefore produces a single genome-wide VCF rather than per-chromosome files in this mode.

### Reference panel behaviour
When `reference_vcf` is set, the pipeline runs:
1. `bcftools isec` — splits target markers into two groups: those **in the reference** and those **only in the target**
2. `conform-gt` — harmonises strand/allele order of the intersection to match the reference
3. `Beagle` with `ref=` — phases and imputes only the intersection markers (Beagle drops target-only markers by design)
4. **Merge** — target-only markers are merged back with the Beagle output so the final VCF contains all original chip markers plus any markers imputed from the reference

This means you can use a reference panel to improve phasing and imputation of overlapping markers without losing any of your original chip data.

### bref3: faster reference loading
Setting `bref3_jar` causes the pipeline to convert the reference VCF to bref3 binary format before running Beagle. This conversion runs once; subsequent per-chromosome Beagle jobs load it 3–43× faster depending on reference panel size. The `bref3.jar` file is available on the same download page as the Beagle JAR.

### Other notes
- **`ne`** has a large effect on imputation accuracy for non-human populations. Use a value appropriate for your species — the Beagle default of 1,000,000 is designed for large outbred human populations and performs poorly for livestock and aquaculture species.
- **`nthreads`** should match the number of CPUs you allocate. Beagle scales well with more threads — the default of 4 is a reasonable starting point; increase for cluster jobs.
- **`plink_extra_flags`**: Set species-specific PLINK flags here, e.g. `"--dog --aec"` for canine data. Empty by default.
