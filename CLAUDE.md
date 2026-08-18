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

Reads `config.yaml` + `config_accuracy.yaml`. Set `accuracy_mode` to one of:

| Mode | What it does |
|------|--------------|
| `mask_and_impute` | Hold out one validation set, mask it to LD density, impute, compare to truth |
| `cross_array` | Same animals on two arrays, scored K-fold; identity-gated. See below |
| `kfold_mask_and_impute` | Animal-level K-fold CV against a shared LD panel; can benchmark Beagle vs AlphaImpute2 vs FImpute in one run |

### `cross_array` is a K-fold, not a held-out panel

The two CV modes share one rule family (`acc_cv_*`, guarded by
`if _acc_mode in ("kfold_mask_and_impute", "cross_array")`). They differ only in
which fileset each per-fold output is cut from:

| Per-fold output | `kfold_mask_and_impute` | `cross_array` |
|---|---|---|
| `to_impute/masked` | the bfile, thinned to the panel | the **LD array**, that fold's animals |
| `reference/panel` | the bfile, minus the fold | the **HD array**, minus the fold |
| `truth/hd` | the bfile, fold only, panel removed | the **HD array**, fold only, panel removed |
| `to_impute/combined` | HD masked for the fold | HD masked, **the fold's LD calls spliced in** |

An external HD panel is not an option for the pairs this exists to measure.
Practically every Ssa70kv4 fish is also V3-typed, so exactly **one** V4 fish
sits outside the cohort — a held-out fold is the only construction in which a
V3 → V4 test runs at all, and using it everywhere keeps the truth pairs
comparable.

`cross_array` derives its LD panel rather than sampling one: it is exactly
LD ∩ HD (`acc_cv_shared_markers` → `scripts/shared_marker_list.py`). That
forces LD ⊆ HD, so no LD-only probe is handed to the imputer as observed and
then never scored.

**The identity gate is mandatory and runs before fold assignment.** Pairing an
animal's LD and HD genotypes on `individ_id` assumes both typings came from one
animal. When an ID covers two physical samples the pair is two unrelated fish,
and the step scores an imputation that never had a chance. Concordance between
LD and HD at their shared markers separates the two cases cleanly — the same
animal lands at 0.95–1.0, two different salmon near 0.5, nothing in between. On
the 2026-08 benchmark this caught 974 of 3,881 animals (25%) in one step, which
had scored mean R² 0.53 purely because a quarter of its truth was a different
fish.

`acc_cv_pair_animals` → `acc_cv_identity_export` → `acc_cv_identity_metrics` →
`acc_cv_identity_gate` writes `setup/identity.tsv`, `setup/identity_fail.ids`
and `setup/animals.txt`; only `animals.txt` reaches `acc_cv_setup`. The run
aborts outright below `cross_array.min_identity_pass_rate` — half a mispaired
cohort is a data problem, not a result.

**Animals pair on IID, not on (FID, IID).** Per-array PLINK exports carry
`FID=0` and are produced in separate runs, so an (FID, IID) match returns
nothing and the run fails with no hint why. The FID is taken from the HD side.

**Allele orientation is forced onto the HD `.bim`** wherever the LD fileset is
exported (`--ref-allele force <hd>.bim 6 2`) or spliced
(`flip_needed()` in `scripts/mask_validation_genotypes.py`). Two exports of the
same marker can disagree about which allele is A1, and **allelic R² cannot see
the difference** — squaring the correlation hides the sign. It surfaces only as
collapsed concordance, long after the run. Incompatible allele pairs raise
rather than being guessed at.

**The fold's LD genotypes must come from the LD array.** `to_impute/combined`
(what FImpute and AlphaImpute2 eat) is built with
`mask_validation_genotypes.py --replace-from <ld_bfile>`. Reusing the HD array's
calls at the same positions would delete the probe differences and
between-array genotyping error the run exists to measure, quietly turning a
cross-array benchmark into a masking benchmark.

### Choosing markers for a panel

`acc_cv_snp_reliability` (`scripts/aggregate_snp_reliability.py`) fans in every
`metrics_by_snp.tsv` and writes `snp_reliability.tsv` + `reliable_markers.txt`.
Two rules make its answer usable:

- a marker is filtered on its **minimum** R² across folds, not its mean — one
  collapsed fold out of five is exactly what a mean hides;
- passing `--root` more than once **intersects runs**, and a marker absent from
  one run is not eligible however well it did in the others.

Markers are keyed by `CHROM:POS:REF:ALT` in the metrics and mapped back to IDs
through the HD `.bim` by (chrom, pos), because a `--extract` list needs IDs.

**All three modes score imputed markers only.** The truth bfile is built with
`--exclude` on the LD panel, because those markers are handed to the imputer as
observed genotypes and come back unchanged — scoring them would inflate accuracy
by roughly the panel's share of all markers. When adding a new accuracy mode,
make sure its truth excludes whatever the imputer was given as input.

### QC is an input contract, not a pipeline step

Dataset preparation is upstream and deliberately outside this repo. The bfiles
handed to any accuracy mode are expected to be already filtered **per chip,
before merging or pairing**: genotype status `OK`, `--geno 0.02`, `--maf 0.025`,
and `--hwe <p> <k> midp keep-fewhet` (the `k` term scales the threshold with
sample size; `keep-fewhet` drops only excess heterozygosity, because too *few*
hets is what family structure and the Wahlund effect produce and those markers
are real). Two ordering invariants come with it: run QC **before** the pedigree
is written into `.fam`, since plink2's `--hwe` considers founders only; and
apply sample-missingness exclusions to **both** arrays, or the pairing breaks.

Note what `--maf 0.025` costs when it reaches the scored markers and not just
the input panel: the markers it removes are the array-specific probes, which
are systematically low-MAF, and they are exactly the ones being imputed. On the
2026-08 step 1 it cut the scoring set from 19,499 markers to 9,470 and empties
the lowest bins of the stratified report. Accuracy will read higher than a
full-spectrum run. That is the current agreed setting; it is a reporting
decision, not a defect.

Both CV modes are configured under the `cv:` key (`n_folds`, `target_n_snps`,
`target_snp_list`, `random_seed`, `imputers`, `reference_max_animals`,
`reliable_r2_threshold`) plus `fimpute_params:`; `cross_array` adds
`cross_array:` (`ld_bfile`, `hd_bfile`, `identity_threshold`,
`min_identity_pass_rate`, `min_shared_markers`). They write `cv_summary.tsv`
(one row per imputer/fold/metric), `cv_imputer_summary.tsv` (mean/SD by
imputer) and `snp_reliability.tsv` + `reliable_markers.txt`, whereas
`mask_and_impute` writes `summary.tsv` + `metrics_by_{maf_bin,snp,individual}.tsv`.
`rule accuracy_all` switches its target list on the mode.

`cv.reference_max_animals` caps the per-fold reference panel, sampled at
`cv.random_seed`. Sweeping panel size is re-running with different values and a
different `accuracy_output_dir` — there is no matrix machinery, on purpose.

Snakemake rejects dotted keys in `--config`, so CLI overrides use flat aliases
(`cv_n_folds=10`, `cv_target_n_snps=10000`, `cv_imputers="beagle alphaimpute2 fimpute"`,
`cross_array_ld_bfile=…`, `cross_array_hd_bfile=…`).

**A flat alias must be resolved before the YAML block, not after.**
`config_accuracy.yaml` defines every `cv.*` key, so `_nested_config` reading the
block first meant a CLI override was accepted, ignored, and the run silently
used the YAML value — `cv_n_folds=2` produced ten folds. The same ordering
applies to `_acc_nested` in `Snakefile_accuracy`, which needs the same lookup
before the include runs.

**`mask_and_impute` has no FImpute path** and now raises rather than running
Beagle under FImpute's name. FImpute is reachable through the two CV modes,
which pick the engine per fold from `cv.imputers`.

## Reference panel construction (separate Snakefile)

```bash
snakemake --snakefile Snakefile_refpanel --use-conda --cores 8 \
  --config refpanel_bfile=/path/to/qc_Ssa70kv1 refpanel_name=Ssa70kv1
```

Reads `config.yaml` + `config_refpanel.yaml`. Builds one array's panel per run;
the array is a parameter, so the same rules produce the v1, v3 and v4 panels.

```
select → keep → per-chrom VCF → bcftools norm -d snps → PHASE → bref3
```

**This is where `reference_vcf` comes from.** `config.yaml` documents it as
"must be phased" and trusts the caller — nothing in the main pipeline ever
produced one. Beagle enforces it and aborts a chromosome with `unphased or
missing genotype for reference sample <id> at marker <...>` the moment it is
handed a plain PLINK export. `refpanel_phase` runs Beagle with `gt=` and no
`ref=`, which both phases the panel and fills its missing calls; both are
required, so a phase-only tool would not be enough.

`bref3.jar` is auto-downloaded into `bin/` like `beagle.jar`; the main
Snakefile's `onstart` never fetched it because nothing there built a panel.

### Panel composition

`refpanel.max_per_family` caps animals per sire × dam pair. A breeding cohort is
not a diverse sample: full sibs share long haplotype blocks, so the tenth sib
adds almost nothing a panel can use while costing the same phasing time, and it
drags the panel's allele frequencies toward whichever families were bred most.
On step 1's 1,457 animals a cap of 4 removes 474 (33%), a cap of 2 removes 788.

Capping needs a pedigree and `scripts/select_refpanel.py` refuses without one —
an array export straight from the database has `0` in both parent columns, so
there are no families to cap. Resolve them from GPA first.

**Animals with unknown parents are each their own family.** Collapsing them into
one "unknown" group would cap the entire founder set to N animals and discard
the broadest haplotype diversity the panel has.

### Variants, and why exclusions happen before phasing

`refpanel.variants` maps a name to a list of ID files to hold out. `full: []` is
the production panel; a benchmark variant names the animals it must not contain.
One selection pass, one manifest, several panels.

Exclusions are applied at selection, **not** by subsetting an already-phased
panel. The cheap route leaks: the excluded animals' genotypes inform the phase
estimates of the relatives that remain, so an accuracy benchmark built that way
is scored against a panel its own test animals helped build. Each variant is
phased separately; at roughly 25–30 node-minutes genome-wide that is affordable.

Every run writes `manifest.tsv` (one row per animal with the reason it was kept
or dropped) and `panel_report.tsv` (source bfile, pedigree, cap, seed,
exclusions, counts). A bref3 blob carries no provenance of its own.

## Architecture

### Mode flags (set at parse time in `Snakefile`)

Flags controlling which rules and files are included:
- `_use_ref` — `reference_vcf` is set → includes `rules/intersect_and_conform.smk` and routes Beagle input through harmonization
- `_use_bref3` — `bref3_jar` is set → adds a bref3 conversion step before Beagle for faster reference loading
- `_use_alphaimpute2` — `imputer: "alphaimpute2"` → includes `rules/alphaimpute2.smk` instead of Beagle rules
- `_use_fimpute` — `imputer: "fimpute"` → includes `rules/fimpute.smk` instead of Beagle rules

`imputer` is validated at parse time; an unknown value raises rather than
silently falling back to Beagle.

### FImpute mode

`rules/fimpute.smk` runs per chromosome and concatenates, like Beagle. FImpute3
is a licensed binary that conda does **not** install — it is located via
`fimpute_params.executable` and the rule that runs it deliberately has no
`conda:` directive.

Three behaviours share one mechanism, FImpute's chip model:
- **Sporadic missing calls** are filled whatever the setup (code `5` on input).
- **Phasing** comes back through the resolved heterozygote codes `3`/`4`, which
  become `0|1` / `1|0`. Code `1` is a genuinely unresolved het and stays `0/1`.
- **A reference panel** (`fimpute_params.reference_bfile`) makes the run
  two-chip: reference = chip 1 (`ref_chip=1`), target = chip 2, and the SNP info
  table carries one index column per chip with `0` where a chip lacks a marker.
  Each animal's genotype string covers only its own chip's markers, in that
  chip's position order — getting that order wrong shifts every call silently.

Note that `save_genotype` and phasing are mutually exclusive: with it set,
FImpute returns every het as code `1` and the phase is gone. It is therefore
only written to the control file when `phase_output: false`.

#### What FImpute demands of its inputs

FImpute validates hard and aborts the entire run rather than skipping a bad
record, so every one of these has to hold before it will start. All are
enforced in `scripts/fimpute_io.py`; do not weaken them.

| Requirement | What happens otherwise | Where enforced |
|---|---|---|
| Sex stated, not `U` | Infers sex itself, then rejects its own inference: `X appeard as both sire and dam` | `assign_sex_from_role()` |
| No parent in both the sire and dam column | Same error | pedigree build, upstream |
| No ancestry cycles | `Pedigree loop is detected!` | pedigree build, upstream |
| One marker per physical position, across **all** panels jointly | `SNPs with the same physical position found` | `drop_duplicate_positions()` |

Two more traps that fail silently rather than loudly:

- **`plink2 --export A` does not necessarily count a1.** On our exports it
  counts a2 for every marker. The `.raw` column suffix (`<snp>_<allele>`) names
  which allele was counted — read it, never assume. Getting this wrong mirrors
  every genotype, and **allelic R² cannot see it** because squaring the
  correlation hides the sign; only concordance collapses.
- **FImpute drops markers of its own accord** (a target-only marker is logged in
  `excluded_snp_list.txt` as `Not On HD`) and rewrites the SNP map in its output
  folder. Read back *that* map, not the one you supplied, or the genotype string
  and the map disagree by however many it removed.

Ungenotyped parents belong in the pedigree with their own rows. Mapping them to
`0` severs full-sib links whenever the parents were not typed in the same run,
which is the normal case when the reference panel is a different generation.

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
- `rules/refpanel.smk` — phased + bref3 reference panel construction (only used
  via `Snakefile_refpanel`)
- `rules/accuracy.smk` — imputation accuracy evaluation (only used via `Snakefile_accuracy`).
  Holds two parallel rule families: `acc_*` for `mask_and_impute`, and `acc_cv_*`
  shared by `kfold_mask_and_impute` and `cross_array`.
- `scripts/alphaimpute2_to_vcf.py` — converts AlphaImpute2 output to VCF
- `scripts/compute_accuracy_metrics.py` — concordance/r² metrics for accuracy evaluation
- `scripts/make_accuracy_cv_setup.py` — deterministic CV fold assignment + shared LD SNP panel
- `scripts/fimpute_io.py` — PLINK raw ↔ FImpute input/output conversion
- `scripts/aggregate_cv_metrics.py` — aggregates per-fold CV summaries into `cv_summary.tsv`
- `scripts/select_refpanel.py` — panel membership: family capping, exclusions,
  and the manifest that explains both
- `scripts/aggregate_snp_reliability.py` — per-marker R² across folds and runs →
  `snp_reliability.tsv` + `reliable_markers.txt`
- `scripts/pair_animals.py` — matches animals across two arrays on IID
- `scripts/shared_marker_list.py` — intersects two or more `.bim` files; also
  generates the `cv.target_snp_list` for a masked run, so the masked and
  cross-array tests thin to the same markers
- `scripts/identity_gate.py` — turns LD-vs-HD concordance into the animal list
  every cross-array rule is built from
- `scripts/mask_validation_genotypes.py` — byte-level `.bed` masking, plus
  `--replace-from` to splice a second array's calls in (plink2's `--pmerge` is
  unfinished, hence the hand-rolled rewrite)

### Conda environments

- `envs/workflow_env.yaml` — main env: openjdk, bcftools, htslib, pandas, pytest
- `envs/alphaimpute2_env.yaml` — AlphaImpute2 env (Python 3.10 required; 3.12/3.14 incompatible)
- `envs/accuracy_env.yaml` — accuracy evaluation env

### SLURM resource groups and partitions

Rules in the same group (`intersect`, `conform`, `beagle`) are submitted together. The `--group-components` flag controls how many chromosomes per job.

This cluster's partitions are split by node size, so every heavy rule must
declare `slurm_partition` in `resources:` — without it Snakemake submits to
the default `r7i-ondemand-large` (15 GiB / 2 CPU) and sbatch rejects the job
("CPU count per node can not be satisfied").

| Rule | `mem_mb` | `slurm_partition` |
|------|---------:|-------------------|
| `run_beagle` | 70000 | `r7i-ondemand-4xlarge` (124 GiB / 16 CPU) |
| `concat_chromosomes` | 64000 | `r7i-ondemand-4xlarge` |
| `bcftools_isec` / `conform_gt` | 32000 | `r7i-ondemand-2xlarge` (62 GiB / 8 CPU) |
| `merge_imputed_with_target_only` / `convert_ref_to_bref3` | 16000 | `r7i-ondemand-2xlarge` |
| `make_per_chrom_vcf` / `normalize_vcf` / `vcf_to_plink` | (none) | default `r7i-ondemand-large` |
| `acc_cv_run_beagle` | 70000 | `r7i-ondemand-4xlarge` |
| `acc_cv_concat_*` | 64000 | `r7i-ondemand-4xlarge` |
| `acc_cv_alphaimpute2_to_vcf` / `acc_cv_run_fimpute` | 32000 | `r7i-ondemand-2xlarge` |
| `acc_cv_identity_metrics` | 32000 | `r7i-ondemand-2xlarge` |
| `acc_cv_beagle_phase_ref` | 70000 | `r7i-ondemand-4xlarge` |
| `refpanel_phase` | 70000 | `r7i-ondemand-4xlarge` |
| `refpanel_bref3` | 16000 | `r7i-ondemand-2xlarge` |
| `acc_cv_run_alphaimpute2` | 16000 | `r7i-ondemand-2xlarge` |
| `run_alphaimpute2` (main pipeline) | 16000 | `r7i-ondemand-2xlarge` |
| `alphaimpute2_to_vcf` (main pipeline) | 32000 | `r7i-ondemand-2xlarge` |
| `acc_run_beagle` | 70000 | `r7i-ondemand-4xlarge` |
| `acc_concat_imputed` / `acc_concat_alphaimpute2` | 64000 | `r7i-ondemand-4xlarge` |
| `acc_alphaimpute2_to_vcf` | 32000 | `r7i-ondemand-2xlarge` |
| `acc_run_alphaimpute2` | 16000 | `r7i-ondemand-2xlarge` |

All accuracy modes (`mask_and_impute`, `cross_array`, `kfold_mask_and_impute`)
now declare partitions on every heavy rule and run under the SLURM executor.
`acc_cv_identity_metrics` compares the whole paired cohort at every shared
marker in one go, so it is the largest single matrix in a cross-array run —
larger than any per-fold comparison.

Partition names track the cluster's current node generation — as of
2026-08-13 that is `r7i-*` (the `r6i-*` names used earlier no longer exist;
only `r6i-ondemand-24xlarge` survives). Re-check with `sinfo` if sbatch starts
reporting `invalid partition specified`.

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
