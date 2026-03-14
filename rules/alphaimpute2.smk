# rules/alphaimpute2.smk
#
# AlphaImpute2 pedigree + population imputation pipeline.
# Activated when config["imputer"] == "alphaimpute2".
#
# Rule graph:
#   plink_to_alphaimpute2_fmt  →  run_alphaimpute2  →  alphaimpute2_to_vcf
#
# The final output ({output_dir}/alphaimpute2/all_chromosomes.vcf.gz) is
# consumed by the shared vcf_to_plink rule in the main Snakefile.


# ---------------------------------------------------------------------------
# Step 1: Convert PLINK binary to AlphaImpute2 text format
# ---------------------------------------------------------------------------

rule plink_to_alphaimpute2_fmt:
    """
    Export PLINK binary as AlphaImpute2 genotype + pedigree text files.

    Genotype format:  IID g1 g2 ... gN  (space-separated; 0/1/2=dosage, 9=missing)
    Pedigree format:  IID sire_IID dam_IID  (0 = unknown parent)

    Uses plink2 --export A to produce additive-coded dosages (0/1/2/NA for
    each individual × SNP), then reformats with awk.
    """
    input:
        bed = config["bfile"] + ".bed",
        bim = config["bfile"] + ".bim",
        fam = config["bfile"] + ".fam",
    output:
        genotypes = config["output_dir"] + "/alphaimpute2_input/genotypes.txt",
        pedigree  = config["output_dir"] + "/alphaimpute2_input/pedigree.txt",
    params:
        plink      = config["plink_path"],
        bfile      = config["bfile"],
        raw_prefix = config["output_dir"] + "/alphaimpute2_input/raw_tmp",
        extra      = config.get("plink_extra_flags", ""),
    log:
        "logs/plink_to_alphaimpute2.log",
    shell:
        """
        (
        mkdir -p "$(dirname {output.genotypes})"

        # Export additive-coded genotypes:
        #   0/1/2 = copies of A1 allele; NA = missing
        {params.plink} \
            --bfile {params.bfile} \
            --export A \
            {params.extra} \
            --allow-no-sex \
            --out {params.raw_prefix}

        # Reformat to AlphaImpute2 genotype layout:
        #   skip header, keep IID (col 2), map NA→9, output as integers
        awk 'NR>1 {{
            printf $2;
            for (i=7; i<=NF; i++) {{
                if ($i == "NA") printf " 9";
                else printf " %d", $i;
            }}
            printf "\\n";
        }}' {params.raw_prefix}.raw > {output.genotypes}

        # Pedigree: IID  sire_IID  dam_IID  (columns 2 3 4 of .fam; 0 = unknown)
        awk '{{print $2, $3, $4}}' {input.fam} > {output.pedigree}
        ) &> {log}
        """


# ---------------------------------------------------------------------------
# Step 2: Run AlphaImpute2
# ---------------------------------------------------------------------------

rule run_alphaimpute2:
    """
    Run AlphaImpute2 imputation on the full dataset (all chromosomes at once).

    Output file: {out_prefix}.genotypes — same layout as input but with
    imputed values replacing 9 where possible.
    """
    input:
        genotypes = rules.plink_to_alphaimpute2_fmt.output.genotypes,
        pedigree  = rules.plink_to_alphaimpute2_fmt.output.pedigree,
    output:
        genotypes = config["output_dir"] + "/alphaimpute2_output/imputed.genotypes",
    params:
        out_prefix  = config["output_dir"] + "/alphaimpute2_output/imputed",
        cycles      = config["alphaimpute2_params"].get("cycles", 4),
        threshold   = config["alphaimpute2_params"].get("final_peeling_threshold", 0.1),
        hd_thresh   = config["alphaimpute2_params"].get("hd_threshold", 0.95),
        length      = config["alphaimpute2_params"].get("length", 1.0),
        extra_flags = " ".join(
            (["-ped_only"] if config["alphaimpute2_params"].get("ped_only", False) else [])
            + (["-pop_only"] if config["alphaimpute2_params"].get("pop_only", False) else [])
            + (["-phase_output"] if config["alphaimpute2_params"].get("phase_output", False) else [])
        ),
    threads:
        config["alphaimpute2_params"].get("maxthreads", 1)
    conda:
        "../envs/alphaimpute2_env.yaml"
    log:
        "logs/run_alphaimpute2.log",
    resources:
        mem_mb = 16000,
    shell:
        """
        mkdir -p "$(dirname {output.genotypes})"

        (AlphaImpute2 \
            -genotypes {input.genotypes} \
            -pedigree  {input.pedigree} \
            -out       {params.out_prefix} \
            -maxthreads {threads} \
            -cycles    {params.cycles} \
            -final_peeling_threshold {params.threshold} \
            -hd_threshold {params.hd_thresh} \
            -length    {params.length} \
            {params.extra_flags}
        ) &> {log}
        """


# ---------------------------------------------------------------------------
# Step 3: Convert AlphaImpute2 output to bgzipped VCF
# ---------------------------------------------------------------------------

rule alphaimpute2_to_vcf:
    """
    Convert AlphaImpute2 imputed genotypes to a genome-wide bgzipped VCF.

    Uses the original .bim file for chromosome/position/allele metadata.
    Output is sorted by the natural order of SNPs in the .bim file and
    indexed with tabix for downstream bcftools / PLINK2 compatibility.
    """
    input:
        genotypes = rules.run_alphaimpute2.output.genotypes,
        bim       = config["bfile"] + ".bim",
    output:
        vcf = config["output_dir"] + "/alphaimpute2/all_chromosomes.vcf.gz",
        tbi = config["output_dir"] + "/alphaimpute2/all_chromosomes.vcf.gz.tbi",
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/alphaimpute2_to_vcf.log",
    resources:
        mem_mb = 32000,
    script:
        "../scripts/alphaimpute2_to_vcf.py"
