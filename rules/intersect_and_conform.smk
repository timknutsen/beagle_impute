# file: rules/intersect_and_conform.smk

rule bcftools_isec:
    input:
        target    = rules.normalize_vcf.output.vcf,
        tbi       = rules.normalize_vcf.output.tbi,
        reference = config.get("reference_vcf", "")
    output:
        vcf             = config["output_dir"] + "/intersect/chr{chrom}_temp/0002.vcf.gz",
        tbi             = config["output_dir"] + "/intersect/chr{chrom}_temp/0002.vcf.gz.tbi",
        target_only     = config["output_dir"] + "/intersect/chr{chrom}_temp/0000.vcf.gz",
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

        tabix -f {output.vcf}
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
        outbase     = config["output_dir"] + "/harmonized/chr{chrom}"
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/conform_gt_chr{chrom}.log"
    resources:
        mem_mb  = 32000,
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

# Convert one chromosome of the reference VCF to bref3 binary format.
# Beagle reads bref3 3-43x faster than VCF depending on reference panel size.
# Only included when bref3_jar is configured (see _use_bref3 flag in Snakefile).
rule convert_ref_to_bref3:
    input:
        reference = config.get("reference_vcf", "")
    output:
        bref3 = config["output_dir"] + "/bref3/chr{chrom}_ref.bref3"
    params:
        bref3_jar = config.get("bref3_jar", "")
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/bref3_chr{chrom}.log"
    resources:
        mem_mb = 16000
    shell:
        """
        (bcftools view -r {wildcards.chrom} -O z {input.reference} \
            | java -jar {params.bref3_jar} /dev/stdin \
            > {output.bref3}) 2> {log}
        """
