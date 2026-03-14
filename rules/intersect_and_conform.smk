# file: rules/intersect_and_conform.smk
rule bcftools_isec:
    input:
        target = rules.normalize_vcf.output.vcf,
        tbi = rules.normalize_vcf.output.tbi,
        reference = config.get("reference_vcf", "")
    output:
        vcf         = config["output_dir"] + "/intersect/chr{chrom}_temp/0002.vcf.gz",
        tbi         = config["output_dir"] + "/intersect/chr{chrom}_temp/0002.vcf.gz.tbi",
        target_only = config["output_dir"] + "/intersect/chr{chrom}_temp/0000.vcf.gz",
        target_only_tbi = config["output_dir"] + "/intersect/chr{chrom}_temp/0000.vcf.gz.tbi"
    params:
        out_prefix = config["output_dir"] + "/intersect/chr{chrom}_temp"
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/bcftools_isec_chr{chrom}.log"
    resources:
        mem_mb = 32000
    shell:
        """
        (bcftools isec \
            -p {params.out_prefix} \
            -O z \
            {input.target} \
            {input.reference}) 2> {log}

        tabix -f {output.target_only}
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

# Merge Beagle-imputed markers (ref intersection) with target-only markers
# (markers in target not in reference, kept as-is from the chip data)
rule merge_with_target_only:
    input:
        imputed     = config["output_dir"] + "/imputed_ref/chr{chrom}.vcf.gz",
        imputed_tbi = config["output_dir"] + "/imputed_ref/chr{chrom}.vcf.gz.tbi",
        target_only = rules.bcftools_isec.output.target_only,
        target_only_tbi = rules.bcftools_isec.output.target_only_tbi
    output:
        vcf = temp(config["output_dir"] + "/imputed/chr{chrom}.vcf.gz"),
        tbi = temp(config["output_dir"] + "/imputed/chr{chrom}.vcf.gz.tbi")
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/merge_imputed_chr{chrom}.log"
    shell:
        """
        (bcftools concat --allow-overlaps \
            {input.imputed} \
            {input.target_only} | \
        bcftools sort | \
        bgzip > {output.vcf}) 2> {log}

        tabix -f {output.vcf}
        """

if config.get("bref3_jar", "").strip():
    rule convert_ref_to_bref3:
        input:
            ref = config.get("reference_vcf", "")
        output:
            bref3 = config["output_dir"] + "/ref/reference.bref3"
        params:
            bref3_jar = config["bref3_jar"]
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/convert_ref_to_bref3.log"
        shell:
            "java -jar {params.bref3_jar} {input.ref} > {output.bref3} 2> {log}"
