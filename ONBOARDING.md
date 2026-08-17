# Onboarding — beagle_impute

Read `AGENTS.md` (symlink to `CLAUDE.md`) for durable architecture and
conventions. This file is the live state at handover: what works, what does
not, and where the bodies are buried.

**Handover date:** 2026-08-16
**Branch:** `fix/fimpute-real-data-invariants` (1 commit ahead of `master`,
not pushed)

---

## What works

All three imputers run end to end from a PLINK bfile, a pedigree in `.fam`
columns 3–4, and a reference panel:

| Imputer | Reference panel | Pedigree | Parallelism |
|---|---|---|---|
| Beagle 5.5 | `reference_vcf`, must be phased | ignores it | per chromosome |
| FImpute3 | `fimpute_params.reference_bfile` → chip 1 | uses it | per chromosome |
| AlphaImpute2 | merge into the input bfile yourself | uses it | whole genome, one thread |

`imputer:` selects the engine and is validated at parse time. 51 tests pass;
all three DAGs dry-run clean.

The FImpute path was **not** usable against real two-chip data until the commit
on this branch. Six separate defects each aborted the whole run. They are listed
in the commit message and the invariants are now documented in `CLAUDE.md`
under "What FImpute demands of its inputs". Do not weaken those checks — every
one of them corresponds to a run that failed.

## What does not work

**`cross_array` accuracy mode is broken.** It hands Beagle the LD fileset with
no reference panel, so the HD-only markers are absent from the input and there
is nothing to impute. It has never produced a valid number. `kfold_mask_and_impute`
already has the missing piece (`acc_cv_make_reference_bfile`); porting it is the
smallest fix.

**No identity gate.** Nothing verifies that an animal's LD and HD genotypes came
from the same physical sample. On the salmon benchmark this silently ruined a
whole step — see below.

**The test fixtures cannot catch these bugs.** `tests/conftest.py` builds 12
animals × 80 SNPs with no duplicate positions, no two-chip setup, and — the
important part — it encodes the *same wrong assumption* the production code had
about which allele `plink2 --export A` counts. Fixture and code agreed with each
other and both disagreed with reality, which is why 51 green tests coexisted with
a module that could not complete a single real run. One fixture with a duplicate
position, a two-chip layout, and a `.raw` whose counted allele is a2 would have
caught three of the six defects.

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

## Suggested direction for the next session

1. Fix `cross_array` — it is the mode this benchmark actually needs.
2. Add the identity gate as a required step of any cross-array mode, not an
   optional check.
3. Move the dataset-preparation scripts into `scripts/` and parameterise the
   base path.
4. Rebuild the fixtures around the failure modes above.
5. Decide whether AlphaImpute2 stays in the benchmark; at ~12 h per step against
   FImpute's 7 minutes it costs two days of wall time for a result that loses on
   both axes.
