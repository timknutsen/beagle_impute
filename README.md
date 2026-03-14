
# Chromosome-wise Imputation Pipeline

This pipeline converts PLINK files to VCF per chromosome, optionally harmonizes against a reference VCF, and imputes with Beagle 5.5.

## 1. Requirements

- **Conda** and **Snakemake** installed
- **PLINK2** and **Beagle 5.5 JAR** (see paths in `config.yaml`)
- Everything else (bcftools, pandas, Java) is installed automatically via `envs/workflow_env.yaml`

## 2. Configuration

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

You can also override individual values on the command line without editing the file:

```bash
snakemake --use-conda --cores 8 \
          --config bfile=/path/to/mydata \
                   beagle_jar=/path/to/beagle.17Dec24.jar \
                   plink_path=/path/to/plink2
```

## 3. Running the Pipeline

```bash
snakemake --use-conda --cores 8
```

For SLURM cluster execution, see `snakemake_slurm_example.sh` and add `--set-resources` flags for partition names as needed.

## 4. Outputs

- Per-chromosome imputed VCFs: `vcf_output/imputed/chr{chrom}.vcf.gz`
- Single combined VCF + index: `vcf_output/all_chromosomes.vcf.gz` (+ `.tbi`)
- Final PLINK files: `plink_binary/imputed_data.bed/.bim/.fam`
- Logs per step and chromosome: `logs/`

## 5. Notes

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
