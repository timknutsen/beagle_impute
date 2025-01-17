# Chromosome-wise VCF Processing and Beagle Imputation Pipeline

This Snakemake workflow processes PLINK binary files chromosome by chromosome, converts them to VCF format, normalizes them using bcftools, and performs imputation using Beagle.

## Prerequisites

1. Software Requirements:
   - Snakemake (>=8.0)
   - PLINK 2.0
   - bcftools
   - Beagle (v5.4)
   - Java (for Beagle)
   - Python with pandas

2. Required conda environments:

```yaml
# envs/bcftools.yaml
name: bcftools
channels:
  - bioconda
  - conda-forge
dependencies:
  - bcftools
  - htslib

# envs/beagle.yaml
name: beagle
channels:
  - conda-forge
dependencies:
  - openjdk=8
```

## Directory Structure
```
.
├── Snakefile
├── envs/
│   ├── bcftools.yaml
│   └── beagle.yaml
├── logs/
└── vcf_output/
    ├── normalized/
    └── imputed/
```

## Configuration

The workflow configuration is defined in the Snakefile. Key parameters include:

```python
config = {
    "bfile": "path/to/your/plink/files",
    "plink_path": "path/to/plink2",
    "beagle_jar": "path/to/beagle.jar",
    "output_dir": "vcf_output",
    "beagle_params": {
        "window": 80,
        "overlap": 10,
        "ne": 500,
        "nthreads": 9
    }
}
```

## Resource Requirements

The workflow uses three different AWS instance types for different stages:
- PLINK conversion: r6i.xlarge (4 vCPU, 32GB RAM)
- VCF normalization: r6i.2xlarge (8 vCPU, 64GB RAM)
- Beagle imputation: r6i.12xlarge (48 vCPU, 384GB RAM)

## Running the Pipeline

Execute the workflow using:

```bash
snakemake --executor slurm \
          -j135 \
          --use-conda \
          --group-components beagle=5 \
          --cores 48
```

Command explanation:
- `--executor slurm`: Use SLURM for job submission
- `-j35`: Allow up to 35 jobs to run concurrently
- `--use-conda`: Use conda environments for software dependencies
- `--group-components beagle=5`: Group Beagle jobs in batches of 5
- `--cores 48`: Use 48 cores total (matches r6i.12xlarge instance)

## Running with Custom Input Files

To run the workflow with a different input PLINK bed file set:

```bash
# Run with custom bed file input
snakemake --executor slurm \
          -j100 \
          --use-conda \
          --group-components beagle=5 \
          --cores 48 \
          --config bfile=/path/to/your/custom/bedfile
```

Notes:
- Provide the path without the .bed/.bim/.fam extension
- The workflow will automatically look for the corresponding .bim and .fam files
- All other configuration parameters remain as defined in the Snakefile
- Output directories and structure remain the same

## Workflow Steps

1. **make_per_chrom_vcf**: 
   - Converts PLINK binary files to VCF format per chromosome
   - Uses r6i.xlarge instance
   - Runtime: 60 minutes per job

2. **normalize_vcf**: 
   - Normalizes VCF files using bcftools
   - Uses r6i.2xlarge instance
   - Runtime: 120 minutes per job

3. **run_beagle**: 
   - Performs imputation using Beagle
   - Uses r6i.12xlarge instance
   - Groups 5 chromosomes per job
   - 70GB memory and 9 cores per chromosome
   - Runtime: 240 minutes per job

## Output Files

- Raw VCF files: `vcf_output/chr{chrom}.vcf.gz`
- Normalized VCF files: `vcf_output/normalized/chr{chrom}.vcf.gz`
- Imputed VCF files: `vcf_output/imputed/chr{chrom}.vcf.gz`

## Logs

All logs are stored in the `logs/` directory:
- PLINK conversion: `logs/chr{chrom}.log`
- Normalization: `logs/normalize_chr{chrom}.log`
- Beagle imputation: `logs/beagle_chr{chrom}.log`

## Notes

- The workflow automatically excludes chromosome 0
- Local rules (make_per_chrom_vcf and normalize_vcf) run on the head node
- Beagle jobs are grouped to optimize resource usage on the r6i.12xlarge instance