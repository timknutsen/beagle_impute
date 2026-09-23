# Recipe: test imputation Ssa70kv3 → Ssa70kv4 on real data

For an agent with read access to genoDB, the AquaGen filesystem and the SLURM
cluster. Follow the phases in order; each ends in a checkpoint. **If a checkpoint
fails, stop and report — do not tune thresholds or edit pipeline code to get past
it.**

## Goal

Measure how accurately the V4-only markers can be imputed from real v3 typings
of the same fish, using `Snakefile_accuracy` in `cross_array` mode, with Beagle
and FImpute. Deliverables:

1. `accuracy_v3_to_v4/cv_imputer_summary.tsv` — r² and concordance per imputer
2. `accuracy_v3_to_v4/reliable_markers.txt` — V4 markers that impute reliably
3. A short report (template at the end) including every number the checkpoints ask for

## Ground rules

- **Read-only on genoDB.** Queries only; never write, update or delete.
- **Do not edit files in this repo.** Configure everything with `--config`.
  If the pipeline itself looks wrong, stop and report it.
- **Work in one project directory**, e.g.
  `/mnt/efshome/aquagen/projects/imputation_v3_v4/` (call it `$PROJ`), with
  `data/`, `qc/`, `logs/` and the accuracy output inside it. Record every query
  and command you run in `$PROJ/README_run.md`.
- **Read first:** `README.md`, then the `cross_array` and "QC is an input
  contract" sections of `CLAUDE.md`. They explain why each rule below exists.
- Known facts to check your numbers against (from `ONBOARDING.md`, 2026-08):
  1,629 V4-typed fish, 1,621 of them also v3-typed; v3 lacks 9,583 V4 markers
  (13.9% of V4). Large deviations are a finding, not something to fix.

## Phase 1 — Extract genotypes from genoDB

Discover the schema before querying; do not guess table or column names. Find
where array type, genotype status, batch/delivery, typing ID and individual ID
live, and how genotypes are exported to PLINK. Prefer the existing production
extraction (`ggf2.py`) over writing a new exporter; `aqg_db2plink` is fine for
metadata.

1. **Cohort:** fish with at least one `OK` typing on **both** Ssa70kv3 and
   Ssa70kv4. Also note the fish typed on only one of them (counts only).
2. **One typing per fish per array.** Where a fish has several `OK` typings on
   the same array, pick one deterministically (e.g. highest call rate, then
   latest date) and log the rule and the fish affected.
3. **IDs:** the PLINK IID in *both* exports must be the **individual (fish) ID**,
   not the typing/sample ID — the pipeline pairs the arrays on IID. FID may be 0.
4. **Coordinates:** both exports on the same assembly/map (SIMONv31 or whatever
   the DB uses for both). Record which.
5. Export `data/raw_v3.{bed,bim,fam}` and `data/raw_v4.{bed,bim,fam}`.

**Checkpoint 1**

- Count fish in each export and in the overlap (expect ~1,621 on both arrays).
- Count markers per array and the overlap:
  `python scripts/shared_marker_list.py --bim data/raw_v3.bim data/raw_v4.bim --out qc/shared_raw.txt`
  — its log reports shared markers, `position_mismatch` and `allele_mismatch`.
  Expect roughly 86% of V4 to be shared. A few hundred position mismatches is a
  known annotation issue (they are dropped automatically); thousands means the
  two exports use different maps — stop.
- No duplicate IIDs within either `.fam`.

## Phase 2 — QC, per chip, before pairing

Run on each array **separately**, with plink2 and `--dog`:

```bash
plink2 --dog --bfile data/raw_v3 \
  --geno 0.02 --mind 0.10 --maf 0.025 \
  --hwe 1e-6 0.001 midp keep-fewhet \
  --make-bed --out qc/qc_v3
# same for v4 -> qc/qc_v4
```

Then:

1. **Apply sample exclusions to both arrays:** a fish removed by `--mind` on
   either array is removed from both (`--remove`), so the pairing survives.
2. **Only then write the pedigree** into the `.fam` of both filesets (sire and
   dam columns, from GPA). plink2's `--hwe` uses founders only, so QC must run
   before parents are filled in. Keep ungenotyped parents as IDs in the sire/dam
   columns — do not replace them with 0. Sex may stay 0; the FImpute step infers
   it from parent roles.
3. Name the results `qc/v3_final` and `qc/v4_final`.

**Checkpoint 2** — report markers and fish removed by each filter on each array,
the final counts, and the fraction of cohort fish with both parents known.
Stop if QC removes more than ~20% of either array's markers or fish.

## Phase 3 — Smoke test on one chromosome, one fold

Cut both arrays to chromosome 29 and run a single fold. This takes minutes and
catches setup problems before the full run.

```bash
for a in v3 v4; do
  plink2 --dog --bfile qc/${a}_final --chr 29 --make-bed --out qc/${a}_chr29
done

cd /mnt/efshome/aquagen/code/timknu/workflows/beagle_impute   # the repo; git pull master first
snakemake --snakefile Snakefile_accuracy --use-conda \
  --executor slurm --jobs 20 --cores 32 \
  --config accuracy_mode=cross_array \
           accuracy_output_dir=$PROJ/smoke_chr29 \
           cross_array_ld_bfile=$PROJ/qc/v3_chr29 \
           cross_array_hd_bfile=$PROJ/qc/v4_chr29 \
           cv_n_folds=5 cv_folds_to_run=1 \
           cv_imputers="beagle fimpute"
```

**Checkpoint 3**

- `setup/identity.tsv`: the pass rate and the concordance distribution. Expect
  it bimodal — genuine pairs at 0.95–1.0, mismatched ones near 0.5, nothing in
  between. Report how many failed (listed in `setup/identity_fail.ids`). The run
  aborts by itself below 50% pass; above that, >5% failing is still worth
  reporting as a data problem.
- `cv_imputer_summary.tsv` for both imputers. Sanity checks:
  - `n_variants_evaluated` ≈ the chr29 V4 markers **not** on v3 (after QC).
    Near zero means the truth or the panel is wrong — stop.
  - Concordance well above 0.9. **r² reasonable but concordance ≈ 0.5–0.6 means
    an allele flip** — stop and report; do not continue to the full run.
- Logs for any failed job are in `$PROJ/smoke_chr29/logs/`.

## Phase 4 — Full run

Same command on the full `qc/v3_final` / `qc/v4_final`, all five folds, output
to `$PROJ/accuracy_v3_to_v4`, and without `cv_folds_to_run`. Use
`--jobs 35 --cores 48`. Beagle's reference phasing is roughly half the compute;
expect roughly 30–60 minutes of wall time at 35 slots for ~1,600 fish (the 50K → V1 estimate in `ONBOARDING.md`). If
jobs die, rerun the same command — Snakemake resumes; add `--rerun-incomplete`
if it asks.

**Checkpoint 4** — all five folds present for both imputers in `cv_summary.tsv`,
and the fold-to-fold SD of mean r² small (roughly < 0.02). One outlying fold
means something about that fold's animals — report which.

## Phase 5 — Report

Write `$PROJ/REPORT.md` with:

| Item | Value |
|---|---|
| Fish on v3 / on V4 / on both (raw) | |
| After QC / after identity check | |
| Identity: passed, failed, concordance range of each group | |
| Markers: V4 total / shared with v3 (LD panel) / scored (V4-only) | |
| Position / allele mismatches dropped | |
| Beagle: mean r², median r², mean concordance (± SD across folds) | |
| FImpute: the same | |
| Markers in `reliable_markers.txt` (worst fold r² ≥ 0.90), and % of the V4-only markers | |
| Accuracy by MAF bin (from `metrics_by_maf_bin.tsv`) — does it collapse below MAF 0.05? | |
| Wall time and failed/re-run jobs | |

Close with three sentences: which imputer wins and by how much, how much of V4
is reliably reachable from v3, and anything anomalous found along the way.
Include the exact commands used and the location of every output.

## Stop and ask, don't improvise, if

- the overlap is far from ~1,621 fish or ~86% of V4 markers;
- more than ~5% of pairs fail the identity check;
- concordance is low while r² looks fine (allele flip);
- `n_variants_evaluated` is near zero;
- FImpute aborts (the log names the offending record — see the FImpute table in
  `CLAUDE.md`), or any step needs a code change to proceed.
