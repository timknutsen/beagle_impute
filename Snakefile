localrules: make_per_chrom_vcf, normalize_vcf, concat_chromosomes, vcf_to_plink

# Configuration with your specific paths
config = {
    "bfile": "/mnt/efshome/aquagen/projects/AG_global_breeding/Rainbowtrout/deformity_RT23_2024/deform_main_analysis/genotypes/raw_genos_parents/parents_offspring_combined",
    "plink_path": "/mnt/efshome/home/timknu/bioinf_tools/plink2.0/plink2",
    "beagle_jar": "/mnt/efshome/home/timknu/bioinf_tools/beagle5.5/beagle.17Dec24.224.jar",
    "tabix_path": "/mnt/efshome/applications/miniconda3/envs/bioinf_tools/bin/tabix",
    "output_dir": "vcf_output",
    "beagle_params": {
        "window": 80,
        "overlap": 10,
        "ne": 500,
        "nthreads": 9
    }
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
        directory("plink_binary") + "/imputed_data.bed"

rule make_per_chrom_vcf:
    input:
        bed = config["bfile"] + ".bed",
        bim = config["bfile"] + ".bim",
        fam = config["bfile"] + ".fam"
    output:
        vcf = temp(directory(config["output_dir"]) + "/chr{chrom}.vcf.gz"),
        vcf_logs = temp(directory(config["output_dir"]) + "/chr{chrom}.log")
    params:
        bfile = config["bfile"],
        plink = config["plink_path"]
    log:
        directory("logs") + "/chr{chrom}.log"
    resources:
        slurm_partition="r6i-ondemand-xlarge",
        mem_mb=32000,
        runtime=60
    shell:
        """
        mkdir -p {config[output_dir]} logs
        
        ({params.plink} \
            --nonfounders \
            --allow-no-sex \
            --bfile {params.bfile} \
            --chr {wildcards.chrom} \
            --export vcf bgz \
            --dog \
            --aec \
            --out {config[output_dir]}/chr{wildcards.chrom} \
            --snps-only just-acgt) &> {log}
        """

rule normalize_vcf:
    input:
        vcf = rules.make_per_chrom_vcf.output.vcf
    output:
        vcf = temp(directory(config["output_dir"]) + "/normalized/chr{chrom}.vcf.gz")
    conda:
        "envs/bcftools.yaml"
    log:
        directory("logs") + "/normalize_chr{chrom}.log"
    resources:
        slurm_partition="r6i-ondemand-2xlarge",
        runtime=120

    shell:
        """
        mkdir -p {config[output_dir]}/normalized logs
        
        (bcftools norm -d snps {input.vcf} | \
         bgzip > {output.vcf}) 2> {log}
        """

rule run_beagle:
    input:
        vcf = rules.normalize_vcf.output.vcf
    output:
        vcf = temp(directory(config["output_dir"]) + "/imputed/chr{chrom}.vcf.gz"),
        tbi = temp(directory(config["output_dir"]) + "/imputed/chr{chrom}.vcf.gz.tbi"),
        impute_logs = temp(directory(config["output_dir"]) + "/imputed/chr{chrom}.log")
    params:
        beagle = config["beagle_jar"],
        window = config["beagle_params"]["window"],
        overlap = config["beagle_params"]["overlap"],
        ne = config["beagle_params"]["ne"],
        outbase = lambda wildcards: f"{config['output_dir']}/imputed/chr{wildcards.chrom}",
        tabix = config["tabix_path"]
    threads: config["beagle_params"]["nthreads"]
    conda:
        "envs/beagle.yaml"
    log:
        directory("logs") + "/beagle_chr{chrom}.log"
    resources:
        mem_mb=70000,  # 70GB per job
        runtime=240,
        slurm_partition="r6i-ondemand-12xlarge",
        cpus_per_task=9
    group:
        "beagle"
    shell:
        """
        mkdir -p {config[output_dir]}/imputed logs

        (java -Xmx{resources.mem_mb}m -jar {params.beagle} \
            gt={input.vcf} \
            window={params.window} \
            overlap={params.overlap} \
            out={params.outbase} \
            nthreads={threads} \
            ne={params.ne} \
            chrom={wildcards.chrom} && \
        {params.tabix} -f {output.vcf}) &> {log}
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
        "envs/bcftools.yaml"
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
        bed = directory("plink_binary") + "/imputed_data.bed",
        bim = directory("plink_binary") + "/imputed_data.bim",
        fam = directory("plink_binary") + "/imputed_data.fam"
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
            --double-id) &> {log}
        """