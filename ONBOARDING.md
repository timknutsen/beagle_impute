# Onboarding — beagle_impute

`AGENTS.md` (symlink to `CLAUDE.md`) holds durable architecture and conventions.
This file is the live state: what works, what does not, what is verified, and
where the bodies are buried.

**Handover date:** 2026-08-18
**Branch:** `feature/cross-array-kfold`, 4 commits ahead of `master`, **not
pushed, no PR open**
**Tests:** 74 pass. All eight DAGs dry-run clean.

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

Four commits; read their messages, they carry the reasoning.

1. **`cross_array` rebuilt** as an identity-gated K-fold sharing the `acc_cv_*`
   family with `kfold_mask_and_impute`. It had never produced a valid number.
2. **ONBOARDING correction** — the 0.955 target was wrong (see "Do not expect").
3. **Beagle reference phasing + cohort gating** — two defects that only appear
   when an imputer actually runs.
4. **`Snakefile_refpanel`** — phased, bref3-encoded reference panels. This is
   where `reference_vcf` is supposed to come from; nothing produced one before.

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

5. **The pedigree covers 1,457 animals, not the 6,331 on the array.** Family
   capping in `Snakefile_refpanel` needs GPA resolved for all of them, and
   `select_refpanel.py` refuses to cap without a pedigree rather than silently
   building an uncapped panel.

6. **AlphaImpute2's place is still unresolved.** ~12 h per step against
   FImpute's minutes, losing on both axes. Still wired into every mode.

7. **`mask_and_impute` has no FImpute path.** It now raises instead of silently
   running Beagle under FImpute's name, but the path does not exist. Use a CV mode.

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
4. Produce `qc_Ssa70kv1` and the full-array pedigree, then build the v1 panel
   and measure capped vs uncapped against a held-out 50K set.
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
