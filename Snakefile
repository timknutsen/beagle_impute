localrules: make_per_chrom_vcf, normalize_vcf, concat_chromosomes, vcf_to_plink

# Configuration with your specific paths
config = {
    # Paths that should be changed
    "bfile": "tests/data/test_salmon",
    "reference_vcf": "tests/data/test_salmon_ref.PHASED.vcf.gz",  # Optional path to reference panel VCF

    # Default paths and parameters
    "plink_path": "/mnt/efshome/home/timknu/bioinf_tools/plink2.0/plink2",
    "beagle_jar": "/mnt/efshome/home/timknu/bioinf_tools/beagle5.5/beagle.17Dec24.224.jar",
    "output_dir": "vcf_output",
    "beagle_params": {
        "window": 80,
        "overlap": 10,
        "ne": 500,
        "nthreads": 1
    },
    "conform_gt_jar": "/mnt/efshome/home/timknu/bioinf_tools/conform-gt.24May16.cee.jar" # Path to conform-gt jar
}

# Helper function to get chromosomes from bim file
def get_chromosomes(bfile):
    import pandas as pd
    bim = pd.read_csv(f"{bfile}.bim", sep='\t', header=None)
    chroms = [chr for chr in bim[0].unique() if chr != 0]
    return sorted(chroms, key=int)

# Get list of chromosomes
CHROMS = get_chromosomes(config["bfile"])

rule all:
    input:
        expand(
            "{output_dir}/imputed/chr{chrom}.vcf.gz",
            output_dir=config["output_dir"],
            chrom=CHROMS
        ),
        config["output_dir"] + "/all_chromosomes.vcf.gz",
        "plink_binary/imputed_data.bed"

rule make_per_chrom_vcf:
    input:
        bed = config["bfile"] + ".bed"
    output:
        vcf = temp(config["output_dir"] + "/dedup" + "/chr{chrom}.vcf.gz"),
        logs = temp(config["output_dir"] + "/dedup" + "/chr{chrom}.log")
    shadow: "shallow"
    params:
        plink = config["plink_path"],
        bfile = config["bfile"],
        out_prefix = config["output_dir"] + "/dedup/chr{chrom}",
        chr = "{chrom}"
    log:
        "logs/dedup_chr{chrom}.log"
    resources:
        slurm_partition="r6i-ondemand-xlarge",
    shell:
        """        
        ({params.plink} \
            --nonfounders \
            --allow-no-sex \
            --bfile {params.bfile} \
            --chr {params.chr} \
            --export vcf bgz \
            --dog \
            --aec \
            --out {params.out_prefix} \
            --snps-only) &> {log}
        """

rule normalize_vcf:
    input:
        vcf = rules.make_per_chrom_vcf.output.vcf
    output:
        vcf = temp(config["output_dir"] + "/normalized/chr{chrom}.vcf.gz"),
        tbi = temp(config["output_dir"] + "/normalized/chr{chrom}.vcf.gz.tbi")
    conda:
        "envs/workflow_env.yaml"
    log:
        "logs/normalize_chr{chrom}.log"
    resources:
        slurm_partition="r6i-ondemand-2xlarge",
    shell:
        """        
        (bcftools norm -d snps {input.vcf} | \
         bgzip > {output.vcf}) 2> {log}

         tabix -f -p vcf {output.vcf}
        """

if config.get("reference_vcf", "").strip():  # If ref provided, remove unique target variants and flip strand in target to match reference
    include: "rules/intersect_and_conform.smk"
    
rule run_beagle:
    input:
        vcf = branch(
            lambda _: config.get("reference_vcf", "").strip(),
            then=config["output_dir"] + "/harmonized/chr{chrom}.vcf.gz",
            otherwise=rules.normalize_vcf.output.vcf
        )
    output:
        vcf = temp(config["output_dir"] + "/imputed/chr{chrom}.vcf.gz"),
        tbi = temp(config["output_dir"] + "/imputed/chr{chrom}.vcf.gz.tbi"),
        impute_logs = temp(config["output_dir"] + "/imputed/chr{chrom}.log")
    params:
        beagle = config["beagle_jar"],
        window = config["beagle_params"]["window"],
        overlap = config["beagle_params"]["overlap"],
        ne = config["beagle_params"]["ne"],
        outbase = lambda wildcards: f"{config['output_dir']}/imputed/chr{wildcards.chrom}",
        ref_param = f"ref={config['reference_vcf']}" if config.get("reference_vcf", "").strip() else ""
    threads: config["beagle_params"]["nthreads"]
    conda:
        "envs/workflow_env.yaml"
    log:
        "logs/beagle_chr{chrom}.log"
    resources:
        mem_mb = 70000,
        runtime = 240,
        slurm_partition = "r6i-ondemand-12xlarge",
        cpus_per_task = 9
    group:
        "beagle"
    shell:
        """
        (
            java -Xmx{resources.mem_mb}m -jar {params.beagle} \
                gt={input.vcf} \
                {params.ref_param} \
                window={params.window} \
                overlap={params.overlap} \
                out={params.outbase} \
                nthreads={threads} \
                ne={params.ne} \
                chrom={wildcards.chrom}
        ) &> {log}

        tabix -f {output.vcf}
        """

# Add the new concatenation rule
rule concat_chromosomes:
    input:
        vcfs = expand(
            "{output_dir}/imputed/chr{chrom}.vcf.gz",
            output_dir=config["output_dir"],
            chrom=CHROMS
        ),
        tbi = expand(
            "{output_dir}/imputed/chr{chrom}.vcf.gz.tbi",
            output_dir=config["output_dir"],
            chrom=CHROMS
        )

    output:
        vcf = config["output_dir"] + "/all_chromosomes.vcf.gz"
    conda:
        "envs/workflow_env.yaml"
    threads: 4
    resources:
        slurm_partition="r6i-ondemand-2xlarge",
        mem_mb=64000,
        runtime=120
    log:
        "logs/concat_chromosomes.log"
    shell:
        """
        (bcftools concat \
            --output {output.vcf} \
            --output-type z \
            --threads {threads} \
            {input.vcfs}) 2> {log}
        """

# Add new rule for VCF to PLINK conversion
rule vcf_to_plink:
    input:
        vcf = config["output_dir"] + "/all_chromosomes.vcf.gz"
    output:
        bed = "plink_binary/imputed_data.bed",
        bim = "plink_binary/imputed_data.bim",
        fam = "plink_binary/imputed_data.fam"
    params:
        plink = config["plink_path"],
        out_prefix = lambda wildcards, output: output.bed.rsplit('.', 1)[0]
    log:
        "logs/vcf_to_plink.log"
    shell:
        """
        ({params.plink} \
            --vcf {input.vcf} \
            --make-bed \
            --out {params.out_prefix} \
            --dog \
            --allow-no-sex \
            --const-fid) &> {log}
        """
