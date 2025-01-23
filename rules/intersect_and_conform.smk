# file: rules/intersect_and_conform.smk
rule bcftools_isec:
    input:
        target = rules.normalize_vcf.output.vcf,
        tbi = rules.normalize_vcf.output.tbi,
        reference = config.get("reference_vcf", "")
    output:
        vcf = config["output_dir"] + "/intersect/chr{chrom}_temp/0002.vcf.gz",
        tbi = config["output_dir"] + "/intersect/chr{chrom}_temp/0002.vcf.gz.tbi"
    params:
        out_prefix = config["output_dir"] + "/intersect/chr{chrom}_temp"
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/bcftools_isec_chr{chrom}.log"
    resources:
        slurm_partition = "r6i-ondemand-xlarge",
        mem_mb = 32000
    shell:
        """       
        (bcftools isec \
            -p {params.out_prefix} \
            -O z \
            {input.target} \
            {input.reference}) 2> {log}
        """

rule conform_gt:
    input:
        vcf = rules.bcftools_isec.output.vcf,
        ref = config.get("reference_vcf", "")
    output:
        vcf = temp(config["output_dir"] + "/harmonized/chr{chrom}.vcf.gz")
    params:
        conform_jar = config.get("conform_gt_jar", ""),
        outbase = config['output_dir'] + "/harmonized/chr{chrom}"
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/conform_gt_chr{chrom}.log"
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
            tabix -f -p vcf {output.vcf}
        else
            echo "No conform-gt .vcf.gz output found for chr{wildcards.chrom}" >> {log}
            exit 1
        fi
        """
