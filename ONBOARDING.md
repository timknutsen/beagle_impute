# Onboarding — beagle_impute

Read `AGENTS.md` (symlink to `CLAUDE.md`) for durable architecture and
conventions. This file is the live state at handover: what works, what does
not, and where the bodies are buried.

**Handover date:** 2026-08-18
**Branch:** `feature/cross-array-kfold`

---

## What works

All three imputers run end to end from a PLINK bfile, a pedigree in `.fam`
columns 3–4, and a reference panel:

| Imputer | Reference panel | Pedigree | Parallelism |
|---|---|---|---|
| Beagle 5.5 | `reference_vcf`, must be phased | ignores it | per chromosome |
| FImpute3 | `fimpute_params.reference_bfile` → chip 1 | uses it | per chromosome |
| AlphaImpute2 | merge into the input bfile yourself | uses it | whole genome, one thread |

`imputer:` selects the engine and is validated at parse time. 69 tests pass; all
DAGs dry-run clean, including `cross_array` across all three imputers.

The FImpute path was **not** usable against real two-chip data until the
2026-08-16 work. Six separate defects each aborted the whole run; the invariants
are documented in `CLAUDE.md` under "What FImpute demands of its inputs". Do not
weaken those checks — every one of them corresponds to a run that failed.

The fixtures now express the failures that matter. `synth_two_chip` in
`tests/conftest.py` carries a low- and a high-density fileset over shared
animals, chip-only markers on both sides, a duplicate physical position, three
markers with A1/A2 reversed between the chips, and two animals whose two typings
are deliberately different fish. `synth_raw_counting_a2` is a `.raw` whose
counted allele is a2, which is what our real exports produce — the old fixtures
encoded the *same wrong assumption* the production code had, which is how 51
green tests coexisted with a module that could not complete a single real run.

`cross_array` is now a real mode. It shares the `acc_cv_*` family with
`kfold_mask_and_impute`, cutting the LD input from one array and the reference
panel and truth from the other, fold by fold. A held-out fold rather than an
external panel is not a convenience: practically every Ssa70kv4 fish is also
V3-typed, so only one V4 fish exists outside the cohort and no external V4 panel
can be built. `CLAUDE.md` has the construction.

The identity gate runs before fold assignment and is not optional. On the
two-chip fixture it separates the planted mispairs at 0.29–0.38 concordance from
everyone else at 0.95–1.0, and a run aborts below
`cross_array.min_identity_pass_rate`.

`snp_reliability.tsv` / `reliable_markers.txt` are the per-marker output the V4
panel is chosen from — filtered on each marker's worst fold, and intersected
across runs when `--root` is passed more than once.

## What does not work

**AlphaImpute2's place in the benchmark is unresolved.** ~12 h per step against
FImpute's 7 minutes, losing on both accuracy and cost. It is still wired into
every mode; nobody has decided to drop it.

**`mask_and_impute` still has no FImpute path.** It now raises instead of
silently running Beagle under FImpute's name, but the path itself does not
exist. Use a CV mode.

**The dataset-preparation scripts still live outside the repo** with hardcoded
paths — see below. That was deliberately left out of this change: prep is an
upstream process, not part of the imputation pipeline.

## The work that lives outside this repo

`/mnt/efshome/aquagen/projects/stepwise_imputation_50k_70k/` holds ~17 scripts
that do everything between "a database" and "a bfile the pipeline can eat":

cohort and ancestor resolution from genodb + GPA, pedigree construction with the
four FImpute invariants, reference-panel assembly, cross-array dataset matching,
per-array export with duplicate-cluster handling, and the identity gate
(`verify_identity.sh`).

**Every one of them hardcodes that project path.** The logic is general; the form
is not. This is the single biggest obstacle to the repo being a reusable tool,
and it is where all six FImpute bugs and the identity problem were actually
found. Moving them into `scripts/` with the base path as an argument is the
highest-value refactor available.

## Live state of the salmon benchmark

Stepwise 50K → 70Kv1 → v2 → v3 → v4, three imputers, accuracy and speed.
Working dir as above; see its `README.md` for the full write-up.

Step 1 (50K → 70Kv1), 1,457 fish, 9,466 markers scored, 5,790-animal panel:

| Imputer | mean R² | concordance | node-minutes |
|---|---:|---:|---:|
| Beagle 5.5 | 0.955 | 0.981 | 38 |
| FImpute3 | 0.924 | 0.970 | 7 |
| AlphaImpute2 | pending | pending | ~716 |

Beagle leads on accuracy; FImpute costs ~6× less compute for ~2 points of R².
AlphaImpute2 was still running at handover (job 415, ~12 h, on the ancestor
panel) and is not competitive on either axis.

**Two results worth carrying forward:**

*Ancestors in the panel help Beagle and hurt FImpute.* Adding 948 genotyped
ancestors moved Beagle +0.0013 and FImpute −0.0064. Pedigree depth made no
difference at all (0.9246 one generation vs 0.9244 two), so the effect is the
panel's composition, not the pedigree. Hypothesis, untested: Beagle's HMM gains
from more haplotype diversity while FImpute's long-block matching is diluted by it.

*Step 2 is not a result, it is a data problem.* It scored mean R² 0.53 because
**974 of its 3,881 animals (25%) have LD and HD genotypes from different
physical fish**. The distribution is sharply bimodal at 0.4–0.6 versus 0.9–1.0.
Steps 1 and 3 have 10 and 72 such animals. Any conclusion from step 2 before
those animals are removed is worthless.

Steps 2–3 have ancestor panels built and Beagle/FImpute scored (step 2 invalid
per above). **Step 4 cannot be run in this design at all**: only one Ssa70kv4
fish exists outside the cohort, because each array version was used for one
period and an animal's parents were therefore typed on the *previous* array. It
needs a held-out-fold reference or it stays out.

## The three truth pairs this was built for

Target: a complete reference panel of V4-genotyped salmon, which needs to know
which V4 markers impute consistently across datasets.

| Test | Mode | Notes |
|---|---|---|
| 50K → V1 | `cross_array` | method benchmark; scores V1 markers. Reproduces the 2026-08 numbers if it is right |
| V3 → V4 | `cross_array` | scores V4 markers on a real array transition |
| 50K-masked → V4 | `kfold_mask_and_impute` | `cv.target_snp_list` = the 50K ∩ V4 markers from `scripts/shared_marker_list.py` |

Only the last two score V4 markers, so the reliable-marker list is their
intersection — pass both output directories to
`scripts/aggregate_snp_reliability.py --root`.

## Suggested direction for the next session

1. Run all three truth pairs on real data. Nothing here has been through Beagle
   or FImpute yet: the rules, the identity gate, the splice and the allele
   orientation are verified on a synthetic two-chip fixture and on plink2, but
   no imputer has run. Check 50K → V1 against the 2026-08 numbers before
   trusting anything else.
2. Sweep `cv.reference_max_animals` to find where accuracy plateaus.
3. Move the dataset-preparation scripts into `scripts/` and parameterise the
   base path — still the biggest obstacle to this being a reusable tool.
4. Decide whether AlphaImpute2 stays in the benchmark.
