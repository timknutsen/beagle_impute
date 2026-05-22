# AlphaImpute2 Concordance Investigation

Date: 2026-05-22

## Question

Why did the 1000-sample real RT24 AlphaImpute2 accuracy run report high allelic r2 but low mean concordance?

Observed summary from `/tmp/ai2_accuracy_1000/accuracy_output/summary.tsv`:

```text
n_samples              200
n_variants_truth       55737
n_variants_evaluated   55737
mean_allelic_r2        0.976178598575565
median_allelic_r2      1.0
pct_r2_ge_0.8          97.9909064024841
pct_r2_ge_0.9          95.27391690078368
mean_concordance       0.35575992
```

## Short Answer

The low concordance is not evidence that AlphaImpute2 imputed badly. It is caused by reversed genotype dosage coding in the AlphaImpute2-to-VCF conversion.

The imputed VCF is effectively using the opposite REF/ALT dosage direction from the truth VCF. Allelic r2 remains high because the metric squares the correlation and therefore hides a negative sign. Exact genotype concordance does not hide that sign, so homozygotes are counted as mismatches.

## Evidence

I compared the current imputed VCF against the truth VCF using the same 200 validation animals and 55,737 evaluated variants.

Files:

- Imputed VCF: `/tmp/ai2_accuracy_1000/accuracy_output/alphaimpute2/all_chromosomes.vcf.gz`
- Truth VCF: `/tmp/ai2_accuracy_1000/accuracy_output/truth/all_chromosomes.vcf.gz`
- Metrics: `/tmp/ai2_accuracy_1000/accuracy_output/summary.tsv`

Counts:

```text
imputed samples:      1000
truth/eval samples:   200
evaluated variants:   55737
```

Diagnostic results:

```text
mean original concordance:       0.3557599159793241
mean flipped-dose concordance:   0.9918251023582954
mean best concordance:           0.991825336886444
signed r mean:                  -0.9865117736403604
signed r median:                -1.0
SNPs with negative signed r:     97.04863914455389%
SNPs where flipped > original:   55736 / 55737
```

This is decisive: if imputed dosages are transformed as `2 - imputed_dosage`, concordance rises from 0.356 to 0.992. Also, 55,736 of 55,737 SNPs have better concordance after flipping.

The MAF-bin pattern also matches a dosage flip:

```text
MAF bin       original concordance   flipped concordance   negative-r SNPs
0.00-0.01     0.0017                 0.9987                17.2%
0.01-0.05     0.0586                 0.9968                96.7%
0.05-0.10     0.1398                 0.9959                99.9%
0.10-0.20     0.2583                 0.9936                100.0%
0.20-0.50     0.4522                 0.9901                100.0%
```

Low-MAF sites look especially bad under the wrong coding because most correct homozygous-major calls are counted as the opposite homozygote.

## Root Cause

The likely source is `scripts/alphaimpute2_to_vcf.py`.

Current assumptions in that script:

```text
0 = 0 copies of A1 -> REF/REF -> VCF 0/0
1 = 1 copy  of A1 -> REF/ALT -> VCF 0/1
2 = 2 copies of A1 -> ALT/ALT -> VCF 1/1
```

Current mapping:

```python
GT_MAP = {"0": "0/0", "1": "0/1", "2": "1/1", "9": "./."}
```

The empirical validation shows this is reversed for the AlphaImpute2 genotype output produced by this workflow. For this run, the mapping should be interpreted in the opposite homozygote direction when writing VCF GT:

```python
GT_MAP = {"0": "1/1", "1": "0/1", "2": "0/0", "9": "./."}
```

Equivalently, current imputed dosage must be evaluated as `2 - dosage` to align to the truth VCF.

## Why r2 Was High

The metrics script computes signed Pearson r and then stores `r ** 2`.

If the true dosage is `x` and the imputed dosage is approximately `2 - x`, the signed correlation is close to `-1`. Squaring gives r2 close to `1`, so the r2 summary looks excellent even though exact GT concordance is low.

This explains the apparently contradictory result:

- High mean allelic r2: the imputed genotypes track the truth very well, but in reversed dosage orientation.
- Low mean concordance: exact GT comparison treats reversed homozygotes as mismatches.

## Fix Applied

The converter and tests were updated after this investigation:

- `scripts/alphaimpute2_to_vcf.py` now maps AlphaImpute2 `0 -> 1/1`, `1 -> 0/1`, `2 -> 0/0`, `9 -> ./.`.
- `tests/test_convert.py` now asserts this orientation so the bug is covered by regression tests.

## Recommendation

Keep `scripts/alphaimpute2_to_vcf.py` and the converter tests aligned so AlphaImpute2 codes are mapped to the same REF/ALT orientation as the PLINK truth VCF.

The applied mapping is:

```python
GT_MAP = {"0": "1/1", "1": "0/1", "2": "0/0", "9": "./."}
```

After changes to this area, rerun at least:

```bash
pytest tests/test_convert.py tests/test_dryrun.py::TestAlphaImpute2Dryrun -q
snakemake --snakefile Snakefile_accuracy --directory /tmp/ai2_accuracy_1000 \
  --use-conda --conda-prefix /tmp/ai2_smoke/.snakemake/conda \
  --cores 4 --rerun-incomplete --rerun-triggers mtime
```

Expected result after the converter fix: mean concordance should be near 0.992 for this run, while mean allelic r2 should remain near 0.976.

## Notes

This investigation does not point to AlphaImpute2 itself as the cause of the low concordance. It points to the VCF conversion layer between AlphaImpute2 output and the accuracy metric script.
