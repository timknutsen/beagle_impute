#!/bin/bash

# Define the input bfile path
BFILE="/mnt/efshome/aquagen/projects/AG_global_breeding/Rainbowtrout/deformity_RT23_2024/deform_main_analysis/genotypes/raw_genos_parents/parents_offspring_combined"

# Define the Snakefile path
SNAKEFILE="/mnt/efshome/aquagen/code/timknu/workflows/beagle_impute/Snakefile"

# Run Snakemake with SLURM executor
snakemake --use-conda \
    --cores 48 \
    --executor slurm \
    --jobs 35 \
    --group-components intersect=5 \
    --group-components conform=5 \
    --group-components beagle=5 \
    --config bfile="${BFILE}" \
    --config output_dir="deformity_RT23_2024" \
    --snakefile "${SNAKEFILE}" \
   --configfile /mnt/efshome/aquagen/code/timknu/workflows/beagle_impute/config.yaml 

