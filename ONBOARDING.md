# Onboarding — beagle_impute

You are taking over an in-progress session on this repo. Read `AGENTS.md`
(symlinked to `CLAUDE.md`) for the durable architecture + conventions. This
file is the *live* state at handover — what's running, what's broken, what
the next move is.

## Live state at handover (2026-05-22, updated)

- Working dir: `/mnt/efshome/aquagen/code/timknu/workflows/beagle_impute`
- Branch: `master`
- Backup tag: `backup/pre-rebase-2026-05-22` — created before today's
  rebase, safe to delete once you trust the master tip.
- Last commit: `82d5868` "Route heavy rules to size-appropriate SLURM partitions"

### Phase 2 run: `runs/RT24_noref/`

Dataset: `/mnt/efshome/aquagen/projects/AG_global_breeding/Rainbowtrout/RT24/plink/RT24_fish`
(25,300 animals × 55,737 SNPs across 32 chromosomes, no reference panel).

Original handover state said the Beagle step was blocked by grouped SLURM memory
requests. The practical fix was to run Beagle with one chromosome per group
component and a real runtime in minutes:

```bash
snakemake --use-conda --cores 48 --executor slurm --jobs 35 \
    --group-components beagle=1 \
    --default-resources runtime=14400 \
    --rerun-incomplete \
    --config bfile=/mnt/efshome/aquagen/projects/AG_global_breeding/Rainbowtrout/RT24/plink/RT24_fish \
             output_dir=runs/RT24_noref
```

Important Snakemake details verified against current docs:

- `mem_mb` is total job memory.
- Group jobs aggregate resources across grouped components.
- `--group-components` targets group names. The old `intersect=5 conform=5`
  examples are ineffective unless those rules have `group:` declarations.
- Runtime values are in minutes for the SLURM executor path; `runtime=240`
  timed out as 4 minutes, while `runtime=14400` requested 240 minutes.

Progress observed during this session:

| Rule | Done |
|------|------|
| `make_per_chrom_vcf` | 32/32 (intermediates were `temp()`, only 11 .vcf.gz left on disk) |
| `normalize_vcf` | 32/32 |
| `run_beagle` | relaunched as one chromosome per grouped job; check current disk/SLURM state before reporting |
| `concat_chromosomes` | check current DAG state |
| `vcf_to_plink` | check current DAG state |

### AlphaImpute2 and accuracy developments

AlphaImpute2 now has a working real-data smoke/accuracy path.

Real test completed under `/tmp/ai2_accuracy_1000`:

- Input: first 1000 RT24 samples, full 55,737 SNPs.
- Validation: 200 animals masked to 10,000 sampled SNPs.
- Reference: 800 animals left full density.
- Accuracy target: 10k to full density.

Key bugs fixed:

- Accuracy PLINK rules now activate `envs/workflow_env.yaml`.
- Accuracy workflow prepends `--dog` like the main Snakefile.
- AlphaImpute2 PLINK export now activates `envs/workflow_env.yaml`.
- AlphaImpute2 mask-and-impute no longer uses PLINK2 non-concatenating
  `--pmerge`/old `--bmerge`. It uses `scripts/mask_validation_genotypes.py`
  to keep reference animals full density and set validation non-LD genotypes
  to missing in-place.
- Accuracy AlphaImpute2 now runs chromosome-wise and concatenates VCFs.
- `alphaimpute2_params.maxthreads` defaults to 1. AlphaImpute2 0.0.3
  segfaulted with `maxthreads=2` on RT24 chr5; `maxthreads=1` completed.
- `scripts/alphaimpute2_to_vcf.py` emits `##contig` lines for bcftools concat.
- AlphaImpute2 genotype code orientation is reversed relative to the earlier
  converter assumption: `0 -> 1/1`, `1 -> 0/1`, `2 -> 0/0`, `9 -> ./.`.
  Without this, allelic r2 remains high but concordance is artificially low.

Evidence:

- Before orientation fix, mean concordance was 0.35575992.
- Diagnostic flipped-dose concordance was 0.99182510.
- Signed r was negative for 97.0% of SNPs, proving the dosage flip.
- See `reports/alphaimpute2_concordance_investigation.md`.

### Things I changed this session that are easy to miss

- `plink_extra_flags` gets `--dog ` prepended at parse time (Snakefile
  line ~33). Every plink2 rule inherits it. **Don't add `--chr-set` or
  `--cow`/`--horse`/etc. anywhere else.** This unblocks any species with
  chromosome codes ≤ 38.
- All rule outputs now route under `${output_dir}` (incl. plink_binary).
- Beagle JAR auto-downloader bumped to `beagle.27Feb25.75f.jar` (5.5).
  Old `beagle.17Dec24.jar` URL 404s.
- `runs/` is in `.gitignore`.
- `Snakefile_accuracy` now inherits the `--dog` convention too.
- AlphaImpute2 accuracy is chromosome-wise; do not revert to one genome-wide
  AlphaImpute2 call unless the segfault is resolved upstream.

### Commands you'll need

```bash
# Resume the run (Snakemake re-uses normalize_vcf outputs already on disk).
snakemake --use-conda --cores 48 --executor slurm --jobs 35 \
    --group-components beagle=1 \
    --default-resources runtime=14400 \
    --rerun-incomplete \
    --config bfile=/mnt/efshome/aquagen/projects/AG_global_breeding/Rainbowtrout/RT24/plink/RT24_fish \
             output_dir=runs/RT24_noref

# Watch
tail -f runs/RT24_noref/_logs/snakemake.log
squeue -u $USER
sacct -X -S "07:00:00" -u $USER -o JobID,JobName%20,State,Elapsed
ls runs/RT24_noref/imputed/   # imputed VCFs appear here

# AlphaImpute2 real accuracy smoke test config used in this session:
snakemake --snakefile Snakefile_accuracy --directory /tmp/ai2_accuracy_1000 \
    --use-conda --conda-prefix /tmp/ai2_smoke/.snakemake/conda \
    --cores 4 --rerun-incomplete
```

### After Phase 2 finishes

Phase 3 is queued in conversation history but not yet started:
imputation **with reference** on `issue_191_trout_parent_control` using the
**Latemat 2023** phased reference VCF, then accuracy evaluation via
`Snakefile_accuracy`. Dataset paths are in
`~/.claude/projects/-mnt-efshome-aquagen-code-timknu-workflows-beagle-impute/memory/reference_data.md`.

### Things I would clean up if I had more time

- Every rule emits "No wall time information given" — add
  `--default-resources runtime=240` to the snakemake invocation, or
  `resources: runtime=...` per rule.
- `make_per_chrom_vcf` / `normalize_vcf` have no `mem_mb` declared; they
  worked on `r6i-ondemand-large` (15 GiB) for 25k animals, but a real
  number would be cleaner.
- The 200+ failed jobs from earlier (`sacct` shows them) are from the
  pre-`--dog` and pre-partition launches. Ignore them; they're not the
  current run.
