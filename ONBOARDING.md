# Onboarding — beagle_impute

`AGENTS.md` (symlink to `CLAUDE.md`) holds durable architecture and conventions.
This file is the live state: what works, what does not, what is verified, and
where the bodies are buried.

**Handover date:** 2026-08-18
**Branch:** `feature/cross-array-kfold`, six commits ahead of `master` after
this handoff commit; pushed with a draft PR open.
**Tests:** 80 pass. Beagle and FImpute both run end to end with and without a
reference; all accuracy/refpanel DAG tests are clean.

---

## The goal this work serves

**V4 imputation**: build a complete reference panel of V4-genotyped salmon,
which requires knowing *which V4 markers impute consistently across datasets*.
That is measured on three truth pairs:

| Test | Mode | Status |
|---|---|---|
| 50K → V1 | `cross_array` | inputs verified, setup stage run, **never run end to end** |
| V3 → V4 | `cross_array` | not started |
| 50K-masked → V4 | `kfold_mask_and_impute` with `cv.target_snp_list` | not started |

Only the last two score V4 markers, so the reliable-marker list is their
intersection: pass both output directories to
`scripts/aggregate_snp_reliability.py --root`.

---

## What changed on this branch

Six commits including this handoff; read their messages, they carry the
reasoning. The material changes are:

1. **`cross_array` rebuilt** as an identity-gated K-fold sharing the `acc_cv_*`
   family with `kfold_mask_and_impute`. It had never produced a valid number.
2. **ONBOARDING correction** — the 0.955 target was wrong (see "Do not expect").
3. **Beagle reference phasing + cohort gating** — two defects that only appear
   when an imputer actually runs.
4. **`Snakefile_refpanel`** — phased, bref3-encoded reference panels. This is
   where `reference_vcf` is supposed to come from; nothing produced one before.
5. **Main-pipeline reference modes hardened** — Beagle and FImpute settings are
   isolated, explicit empty CLI values are safe, two-chip FImpute emits target
   animals only, and all four paths have DAG and live-runtime coverage.

---

## Verified on real data — do not re-derive

Measured against `/mnt/efshome/aquagen/projects/stepwise_imputation_50k_70k/steps/step1`
(1,457 fish, LD 47,913 ⊆ HD 57,383, 9,470 markers to score).

**Identity gate reproduces the hand-rolled `verify_identity.sh` exactly**: the
same 10 animals, and the distribution is bimodal with nothing in the gap —
10 below 0.6, **zero** between 0.6 and 0.9, 1,447 above 0.9.

**The two arrays genotype the same fish essentially identically.** For the 1,447
genuine pairs: per-animal concordance mean 0.9979 (worst 0.9802); per-marker
mean concordance 0.9979, allelic r² 0.9941; only **8** markers of 47,447 below
0.95, and **zero** below 0.50 — so no allele flips or systematic mismapping.
The between-array noise floor is ~0.2% of calls, which caps achievable
concordance near 0.998.

**Beagle and FImpute both complete.** Fold 1 chr29: 35 s wall for phase +
impute on 4 cores, scoring r² 0.974 / concordance 0.987 against a 1,157-animal
panel. FImpute chr29 in 2 s.

**Panel construction works.** chr29 of a 6,331-animal Ssa70kv1 panel builds in
1 m 50 s on 4 cores, comes back fully phased with zero missing calls, and Beagle
loads the bref3 as `ref=` reporting 6,331 reference samples.

**The main pipeline's four reference combinations run end to end.** On the
bundled three-chromosome salmon fixture, using the real Beagle, conform-gt,
PLINK2 and licensed FImpute3 binaries:

| Engine | Reference | Final samples | Final variants | Result |
|---|---:|---:|---:|---|
| Beagle | none | 100 | 650 | VCF + PLINK complete |
| Beagle | phased VCF (200 fish) | 100 | 3,592 | VCF + PLINK complete |
| FImpute | none | 100 | 650 | VCF + PLINK complete |
| FImpute | 50-fish HD PLINK panel | 50 LD targets | 650 | VCF + PLINK complete; zero reference fish leaked |

The FImpute runs deliberately retained a non-empty stale `reference_vcf`; no
Beagle harmonisation rule entered their DAGs. Snakemake's `reference_vcf=` and
`bref3_jar=` CLI null values are also covered because a live run exposed that
edge case.

### Array facts — constraints, not bugs

| | |
|---|---|
| 50K ∩ v1 | 48,623 shared; 19,499 v1-only to impute (28.6% of v1) |
| v1 ∩ v4 | **0 animals**; v1 lacks 16,728 V4 markers (24.2%) |
| v3 ∩ v4 | 1,621 of 1,629 animals; v3 lacks 9,583 V4 markers (13.9%) |
| Year classes | v1 2014–2019 · v3 2018–2024 · v4 2022–2024 |
| Panel pool for v1 | 6,331 typed + 1,252 ancestors, only 1 overlapping ≈ 7,582 |

**A v1 panel cannot serve the V4 goal.** No shared animals, no shared year
classes, and a quarter of V4 is not on the array at any accuracy. v3 is the only
usable bridge, and 13.9% of V4 is beyond even that.

**453 markers sit at different coordinates on the 50K and v1 exports**, 420 of
them on a *different chromosome*, in contiguous blocks. It is in the source
exports (481 raw between 50K and v1; 452 between v1 and v4), i.e. the database's
per-array annotation. `shared_marker_list.py` drops them from the panel.

---

## What does not work / is not done

1. **No truth pair has been run end to end.** Everything above is one
   chromosome, one fold, plus the setup stage. Nothing has produced a
   `cv_summary.tsv`.

2. **`cross_array` cannot use an external reference panel.** It always builds
   the panel from held-out folds. That is correct and necessary for V3 → V4,
   where no external V4 panel can exist — but it blocks three separate things:
   comparing step 1 like-for-like against the 2026-08 numbers; validating a
   permanent panel against a held-out 50K dataset; and measuring whether family
   capping helps. **This is the highest-value next change**: an optional
   `cross_array.reference_bfile` that replaces the fold-derived panel. Roughly
   an hour.

3. **The 453 position-mismatch markers are scored but should not be.** They are
   correctly kept out of the imputer's input, but they land in the truth, so
   4.6% of the scored set is judged on an annotation disagreement rather than on
   imputation. Conservative, so it does not invalidate anything — but exclude
   them from the truth as well.

4. **`qc_Ssa70kv1` does not exist.** Only `raw_Ssa70kv1` is on disk;
   `qc_arrays.sh` was written but its output is not there. Building a panel from
   `raw_` skips `--geno` / `--hwe` / `--maf` entirely.

5. **The local Ssa70kv1 pedigree is nearly complete, not limited to 1,457 as
   previously written.** It overlaps 6,300/6,331 cohort fish and 6,299 have
   both parents. Across cohort + ancestor candidates it overlaps 6,413/7,582;
   audit or explicitly retain the remaining 1,169 ancestors/founders before a
   production capped panel. `select_refpanel.py` still correctly refuses a cap
   when no pedigree is provided.

6. **AlphaImpute2's place is still unresolved.** ~12 h per step against
   FImpute's minutes, losing on both axes. Still wired into every mode.

7. **`mask_and_impute` has no FImpute path.** It now raises instead of silently
   running Beagle under FImpute's name, but the path does not exist. Use a CV mode.

---

## Deferred Ssa70kv1 reference build — inventory and full plan

This is deliberately **recorded, not launched**. Counts below are the live
database snapshot queried on 2026-08-18; re-query immediately before building
because genodb is mutable.

### Full genodb inventory

| Genotype state | Records | Fish-level interpretation |
|---|---:|---|
| `OK` | 332,304 | 316,998 unique fish with at least one OK Ssa70kv1 typing |
| `LOWCR` | 43,898 | 40,252 fish have LOWCR and no OK typing |
| `DQC_FAILED` | 13,844 | 2,131 fish have DQC_FAILED and no OK typing |
| all statuses | — | 359,381 unique fish; 3,332 fish have multiple v1 records |

The OK records span 159 batches. Operational year classes inferred from the
animal-ID prefix are: 1997=2, 1998=15, 1999=49, 2000=59, 2001=291, 2005=684,
2008=221, 2011=260, 2012=122, 2013=213, 2014=418, 2015=11,598,
2016=52,775, 2017=60,404, 2018=94,845, 2019=89,807, 2020=5,146,
`3322`=87 and invalid/missing=3. These buckets total 316,999 while the database
reports 316,998 distinct IDs under its collation, so preserve a one-ID
alias/case-normalisation discrepancy rather than pretending the prefix is an
authoritative year field.

Breeddb matches 232,928 of the 316,998 OK fish; 84,070 do not match. The
available metadata is useful but does **not** expose a reliable nucleus flag:

- `family.category_id` is null for every matched fish, so it cannot identify
  nucleus membership.
- `line` is complete among matches: Prodkjerne tilvekst 106,007; Prodkjerne
  robust 44,407; Mowi 47,230; NFA 34,870; Akvaforsk 414.
- Treat the two Prodkjerne lines as an operational, not canonical, nucleus
  proxy: 150,414 fish, with fertilisation years 2013=152, 2014=391,
  2017=69,822, 2018=80,042 and 7 missing.
- Country among matches is Norway 179,522; Chile 52,532; missing 874.

### Local candidate files already examined

- Raw v1 metadata lists 6,396 animals, while `raw_Ssa70kv1` contains 6,331;
  65 were lost when one batch failed SNP-order consistency.
- The ancestor PLINK set contains 1,252 fish; only one overlaps the 6,331, so
  the union is 7,582 candidate animals.
- The existing GPA pedigree overlaps 6,300/6,331 cohort fish (6,299 with both
  parents) and 6,413/7,582 of the cohort-plus-ancestor union.
- Before marker/sample QC, family caps retain 4,254 (cap 2) or 5,379 (cap 4)
  from the 6,331 cohort; on the 7,582 union they retain 5,462 or 6,609.

### Build plan when work resumes

1. Re-query inventory and freeze a provenance table: typing ID, individual ID,
   status, batch, inferred/real year, line, country, and pedigree availability.
   Use the `aqg_db2plink` metadata path for overview, but use the production
   `ggf2.py` extraction path for the genotype files if its batch/retest checks
   remain the trusted standard.
2. Resolve duplicate typings deterministically and investigate the failed
   65-animal batch. Do not silently union records with inconsistent SNP order.
3. Build `qc_Ssa70kv1` upstream of this repo: genotype status `OK`, unique IID,
   SIMONv31 coordinates, delivery/all-missing checks, then per-chip
   `--geno 0.02 --maf 0.025 --hwe <p> <k> midp keep-fewhet`. Apply QC before
   writing pedigree into `.fam` and apply sample exclusions consistently.
4. Complete/audit GPA rows for the 7,582 candidates. Preserve ungenotyped
   parents as pedigree rows and treat each genuinely unknown-parent animal as
   its own family.
5. Define panel variants at selection time: `full` (production), `cap4`
   (default candidate for unrelated targets), `cap2`, and any benchmark
   exclusions. Exclude validation animals **before** phasing to prevent leakage.
6. Smoke-test chromosome 29, then run `Snakefile_refpanel` on SLURM for every
   retained variant. Require phased, complete VCF/bref3 plus `manifest.tsv` and
   `panel_report.tsv`; never promote a raw PLINK export as Beagle `ref=`.
7. Add `cross_array.reference_bfile`, then validate panel variants against an
   external held-out 50K truth set with both Beagle and FImpute. Compare panel
   size, allelic R², concordance, evaluated-marker count, MAF bins, and worst
   fold. Choose cap/line/year composition from measured performance, not from
   the Prodkjerne label alone.

The full database pool is far too large and family-skewed to phase blindly.
The present 7,582-fish local union is a tractable first candidate, but the final
reference should be chosen only after the upstream QC and external-reference
validation above.

---

## Gotchas that have each cost a run

- **`temp()` outputs nothing declares as input are swept immediately.** This bit
  three times: identity `.tbi`, the panel `.bim`/`.fam`, the phased ref `.tbi`.
  If a rule reads a file, name it in `input:` even when the tool only needs it
  incidentally.
- **A flat `--config` alias must beat the YAML block.** `config_accuracy.yaml`
  defines every `cv.*` key, so resolving the block first made `cv_n_folds=2`
  accepted, ignored, and silently run as ten folds. Both `_nested_config`
  (`rules/accuracy.smk`) and `_acc_nested` (`Snakefile_accuracy`) now check the
  alias first. Keep that order.
- **Beagle's `ref=` must be phased *and* complete.** A plain PLINK export aborts
  the chromosome. `acc_cv_beagle_phase_ref` / `refpanel_phase` run Beagle with
  `gt=` and no `ref=`, which phases and fills missing calls in one pass.
- **`plink2 --export A` counts a2 on our exports, not a1.** Read the `.raw`
  column suffix. Getting it wrong mirrors every genotype.
- **Allelic R² cannot see an allele flip** — squaring hides the sign. Always
  read concordance alongside it, and check `n_variants_evaluated` is what you
  expect. A flip once gave r² 0.932 with concordance 0.409.
- **Dry-runs prove almost nothing.** Every defect fixed in commits 3 and 4
  dry-ran clean. Run one chromosome before trusting a DAG.

---

## Decisions already made — do not relitigate

- **Repo scope is the accuracy modes and panel construction.** Dataset prep
  stays upstream; the ~17 scripts in `stepwise_imputation_50k_70k/` were
  deliberately not ported.
- **`--maf 0.025` applies everywhere, including the scored markers.** Chosen
  knowing it roughly halves the scoring set and reads optimistic.
- **Panel comparison axis is size** (`cv.reference_max_animals`), not composition.
- **Reference panels: one selection pass, one manifest, several variants.**
  Exclusions are applied *before* phasing, never by subsetting a phased panel.
- **Family capping is the diversity criterion**, and the machinery is per array.

**One caveat on capping, because the justification flips by use case.** For
imputing a *new, unrelated* 50K dataset, redundant sibs contribute haplotypes
the first few already gave, so capping is right. For imputing the *next
generation of the same families*, sibs are the signal and capping removes it —
in fold 1, 85% of targets had a half sib in the panel and chr29 scored r² 0.974
with 1,157 animals, against 0.955 for August's 5,790 non-relatives. Do not apply
one rule to both cases, and settle it by measurement once item 2 above exists.

**Do not expect 50K → V1 to reproduce mean R² 0.955.** That came from an
external 5,790-animal panel. A K=5 held-out fold gives ~1,158 — a fifth of it —
and accuracy scales with panel size, so a lower number is the design working.
What *is* comparable is the identity gate, which already reproduces exactly.

---

## Suggested order for the next session

1. Add `cross_array.reference_bfile` (item 2). It unblocks the most.
2. Run 50K → V1 end to end and read R² and concordance together.
3. Exclude the 453 position-mismatch markers from the truth (item 3).
4. Execute the deferred Ssa70kv1 QC/pedigree/panel plan above, then measure
   capped vs uncapped against a held-out 50K set.
5. Run V3 → V4 and 50K-masked → V4; intersect their marker reliability.
6. Decide whether AlphaImpute2 stays.

```bash
# 50K -> V1, the run that is ready to go
snakemake --snakefile Snakefile_accuracy --use-conda --executor slurm --jobs 35 \
  --config accuracy_mode=cross_array \
           cross_array_ld_bfile=<...>/steps/step1/ld \
           cross_array_hd_bfile=<...>/steps/step1/hd \
           accuracy_output_dir=accuracy_50k_to_v1 \
           beagle_jar=bin/beagle.jar \
           cv_n_folds=5 cv_imputers="beagle fimpute"
```

1,671 jobs; the gate serialises the first ~4 minutes, then 290 Beagle
chromosome-jobs fan out. Expect roughly 30–45 minutes wall at 35 slots, half of
it panel phasing.
