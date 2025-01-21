localrules: make_per_chrom_vcf, normalize_vcf, concat_chromosomes, vcf_to_plink

# Configuration with your specific paths
config = {
    #"bfile": "/mnt/efshome/aquagen/projects/AG_global_breeding/Rainbowtrout/deformity_RT23_2024/deform_main_analysis/genotypes/raw_genos_parents/parents_offspring_combined",
    "bfile": "tests/data/test_salmon",
    "plink_path": "/mnt/efshome/home/timknu/bioinf_tools/plink2.0/plink2",
    "beagle_jar": "/mnt/efshome/home/timknu/bioinf_tools/beagle5.5/beagle.17Dec24.224.jar",
    "tabix_path": "/mnt/efshome/applications/miniconda3/envs/bioinf_tools/bin/tabix",
    "output_dir": "vcf_output",
    "beagle_params": {
        "window": 80,
        "overlap": 10,
        "ne": 500,
        "nthreads": 1
    }
}

# Add to existing config
config.update({
    "reference_vcf": "tests/data/test_salmon_ref.PHASED.vcf.gz",  # Path to reference panel VCF
    "conform_gt_jar": "/mnt/efshome/home/timknu/bioinf_tools/conform-gt.24May16.cee.jar", # Path to conform-gt jar
})

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
        bed = config["bfile"] + ".bed"
    output:
        vcf = temp(directory(config["output_dir"] + "/dedup") + "/chr{chrom}.vcf.gz"),
        logs = temp(directory(config["output_dir"] + "/dedup") + "/chr{chrom}.log")
    shadow: "shallow"
    params:
        plink = config["plink_path"],
        bfile = config["bfile"],
        out_prefix = config["output_dir"] + "/dedup/chr{chrom}",
        chr = "{chrom}"
    log:
        directory("logs") + "/dedup_chr{chrom}.log"
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
            --snps-only just-acgt) &> {log}
        """

rule normalize_vcf:
    input:
        vcf = rules.make_per_chrom_vcf.output.vcf
    output:
        vcf = temp(directory(config["output_dir"]) + "/normalized/chr{chrom}.vcf.gz"),
        tbi = temp(directory(config["output_dir"]) + "/normalized/chr{chrom}.vcf.gz.tbi")
    conda:
        "envs/bcftools.yaml"
    log:
        directory("logs") + "/normalize_chr{chrom}.log"
    resources:
        slurm_partition="r6i-ondemand-2xlarge",
    shell:
        """        
        (bcftools norm -d snps {input.vcf} | \
         bgzip > {output.vcf}) 2> {log}

         tabix -f -p vcf {output.vcf}
        """

if config.get("reference_vcf", "").strip():  # Will be False for empty string or whitespace
    rule bcftools_isec:
        input:
            target = rules.normalize_vcf.output.vcf,
            tbi = rules.normalize_vcf.output.tbi,
            reference = config["reference_vcf"]
        output:
            vcf = directory(config["output_dir"]) + "/intersect/chr{chrom}_temp/0002.vcf.gz",
            tbi = directory(config["output_dir"]) + "/intersect/chr{chrom}_temp/0002.vcf.gz.tbi"
        params:
            out_prefix = config["output_dir"] + "/intersect/chr{chrom}_temp",
            chr = "{chrom}",
            tabix = config["tabix_path"]
        conda:
            "envs/bcftools.yaml"
        log:
            directory("logs") + "/bcftools_isec_chr{chrom}.log"
        resources:
            slurm_partition = "r6i-ondemand-xlarge",
            mem_mb = 32000,
            runtime = 60
        shell:
            """
            mkdir -p {params.out_prefix}
            
            (bcftools isec \
                -p {params.out_prefix} \
                -O z \
                {input.target} \
                {input.reference}) 2> {log}
            """

if config.get("reference_vcf", "").strip():
    rule conform_gt:
        input:
            vcf = rules.bcftools_isec.output.vcf,
            ref = config["reference_vcf"]
        output:
            vcf = temp(directory(config["output_dir"]) + "/harmonized/chr{chrom}.vcf.gz")
        params:
            conform_jar = config.get("conform_gt_jar", ""),
            tabix = config["tabix_path"],
            outbase = lambda wildcards: f"{config['output_dir']}/harmonized/chr{wildcards.chrom}"
        conda:
            "envs/beagle.yaml"
        log:
            directory("logs") + "/conform_gt_chr{chrom}.log"
        resources:
            slurm_partition = "r6i-ondemand-xlarge",
            mem_mb = 32000,
            runtime = 60
        shell:
            """
            java -Xmx{resources.mem_mb}m -jar {params.conform_jar} \
                ref={input.ref} \
                gt={input.vcf} \
                chrom={wildcards.chrom} \
                match=POS \
                out={params.outbase} 2> {log}

            if [ -f "{params.outbase}.vcf.gz" ]; then
                {params.tabix} -f -p vcf {output.vcf}
            else
                echo "No conform-gt .vcf.gz output found for chr{wildcards.chrom}" >> {log}
                exit 1
            fi
            """


rule run_beagle:
    input:
        vcf = branch(
            lambda _: config.get("reference_vcf", "").strip(),
            then=config["output_dir"] + "/harmonized/chr{chrom}.vcf.gz",
            otherwise=rules.normalize_vcf.output.vcf
        )
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
        tabix = config["tabix_path"],
        ref_param = f"ref={config['reference_vcf']}" if config.get("reference_vcf", "").strip() else ""
    threads: config["beagle_params"]["nthreads"]
    conda:
        "envs/beagle.yaml"
    log:
        directory("logs") + "/beagle_chr{chrom}.log"
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

        {params.tabix} -f {output.vcf}
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