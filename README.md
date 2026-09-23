# beagle_impute

Chromosome-wise genotype imputation for PLINK data with **Beagle 5.5**,
**FImpute3** or **AlphaImpute2**, plus K-fold accuracy benchmarking and
reference-panel construction. Built for salmon and trout; works for any species
with chromosome codes 1–38.

| Workflow | Snakefile | Config | Produces |
|---|---|---|---|
| Impute | `Snakefile` | `config.yaml` | imputed VCF + PLINK |
| Measure accuracy | `Snakefile_accuracy` | `config_accuracy.yaml` | r² / concordance per imputer, reliable-marker list |
| Build a reference panel | `Snakefile_refpanel` | `config_refpanel.yaml` | phased VCF + bref3 per chromosome |

Every setting is documented in its config file. Architecture and the traps
behind each design choice are in `CLAUDE.md`; current project state is in
`ONBOARDING.md`.

## Requirements

- Conda and Snakemake 8+
- The imputer you pick: Beagle and conform-gt are auto-downloaded to `bin/`;
  FImpute3 is licensed and must already exist (`fimpute_params.executable`);
  AlphaImpute2 is installed by conda
- Everything else (plink2, bcftools, htslib, Java, pandas) comes from `envs/`

## 1. Impute

```bash
snakemake --use-conda --cores 8 \
  --config bfile=/path/to/data output_dir=my_run
```

`bfile` defaults to the small test fixture, so always set it for a real run.

Choose the engine with `imputer=`:

| `imputer` | Uses | Optional reference |
|---|---|---|
| `beagle` (default) | LD only, no pedigree | `reference_vcf=` phased VCF (a file, or `.../chr{chrom}.vcf.gz`) |
| `fimpute` | pedigree + LD, fast | `fimpute_params.reference_bfile` (second chip, PLINK) |
| `alphaimpute2` | pedigree + LD, slow | none |

Without a reference, Beagle and FImpute phase and fill sporadic missing calls.
With one, they also impute the reference's extra markers. Add `bref3_jar=bin/bref3.jar`
to make Beagle load the reference faster.

Set `beagle_params.ne` for your species (~500 salmon/small livestock, ~1000
cattle/dogs); Beagle's default of 1,000,000 is for humans and hurts accuracy here.

**Outputs** (under `output_dir`, default `vcf_output/`):

- `all_chromosomes.vcf.gz` (Beagle), `fimpute/all_chromosomes.vcf.gz` or
  `alphaimpute2/all_chromosomes.vcf.gz`, each with `.tbi`
- `plink_binary/imputed_data.{bed,bim,fam}`
- `logs/`

**On SLURM** (see `snakemake_slurm_example.sh`):

```bash
snakemake --use-conda --executor slurm --jobs 35 --cores 48 \
  --group-components beagle=5 \
  --config bfile=/path/to/data output_dir=my_run
```

Heavy rules already request the right partition.

## 2. Measure accuracy

Animal-level K-fold cross-validation. Each fold's held-out animals are reduced
to a low-density (LD) marker panel, imputed with the other folds as the
high-density reference, and scored only on the markers they were not given.

```bash
# Simulated LD panel: mask one bfile down to 10k markers
snakemake --snakefile Snakefile_accuracy --use-conda --cores 8 \
  --config bfile=/path/to/hd_data accuracy_output_dir=acc_run \
           cv_n_folds=5 cv_target_n_snps=10000

# Two real arrays, same animals (e.g. V3 -> V4)
snakemake --snakefile Snakefile_accuracy --use-conda --cores 8 \
  --config accuracy_mode=cross_array accuracy_output_dir=acc_v3_v4 \
           cross_array_ld_bfile=/path/to/Ssa70kv3 \
           cross_array_hd_bfile=/path/to/Ssa70kv4 cv_n_folds=5
```

- **Quick check first:** add `cv_folds_to_run=1` to run a single fold.
- **Engines:** `cv.imputers` defaults to Beagle and FImpute; add AlphaImpute2
  with `cv_imputers="beagle fimpute alphaimpute2"`.
- **A real LD array instead of random markers:** build the list with
  `python scripts/shared_marker_list.py --bim ld.bim hd.bim --out panel.txt`
  and pass `cv_target_snp_list=panel.txt`.
- **`cross_array` checks identity first.** Animals whose two typings disagree
  (different fish under one ID) are dropped before scoring and listed in
  `setup/identity_fail.ids`.
- **Inputs must be QC'd upstream**, per chip, before pairing. The pipeline does
  not filter genotypes.

**Outputs** (under `accuracy_output_dir`):

| File | What |
|---|---|
| `cv_imputer_summary.tsv` | mean / SD of r² and concordance per imputer |
| `cv_summary.tsv` | the same, per fold |
| `reliable_markers.txt` | markers whose **worst** fold reaches r² ≥ 0.90 |
| `snp_reliability.tsv` | per-marker r² across folds |
| `{imputer}/fold{N}/metrics_by_{snp,maf_bin,individual}.tsv` | detail |
| `setup/identity.tsv` | `cross_array` only: the identity check per animal |

To keep only markers that hold up in several runs, intersect them:

```bash
python scripts/aggregate_snp_reliability.py \
  --root acc_v3_v4 acc_50k_masked_v4 --imputers beagle fimpute --folds 1 2 3 4 5 \
  --bim /path/to/Ssa70kv4.bim --out v4_reliability.tsv --reliable-out v4_reliable.txt
```

## 3. Build a reference panel

Phases one array into a Beagle reference, optionally capping animals per
full-sib family (needs a pedigree) and holding animals out before phasing.

```bash
snakemake --snakefile Snakefile_refpanel --use-conda --executor slurm --jobs 30 \
  --config refpanel_bfile=/path/to/qc_Ssa70kv1 refpanel_name=Ssa70kv1 \
           refpanel_pedigree=/path/to/pedigree.tsv refpanel_max_per_family=4 \
           output_dir=refpanel_output
```

Writes `<variant>/phased/chr{N}.vcf.gz`, `<variant>/bref3/chr{N}.bref3`, and a
`manifest.tsv` recording why each animal was kept or dropped. Use it with:

```bash
snakemake --use-conda --cores 8 \
  --config bfile=/path/to/target \
           reference_vcf=refpanel_output/full/phased/chr{chrom}.vcf.gz
```

## Tests

```bash
pytest                 # needs pytest + snakemake; tools that are missing are skipped
pytest -m "not slow"   # skip tests that need bgzip/tabix/plink2/snakemake
```

Fixtures are generated at test time; no binary data is stored in git.
