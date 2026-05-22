configfile: "config.yaml"

import pandas as pd

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def get_chromosomes(bfile):
    bim = pd.read_csv(f"{bfile}.bim", sep='\t', header=None)
    chroms = [chr for chr in bim[0].unique() if chr != 0]
    return sorted(chroms, key=int)

def get_chroms():
    return get_chromosomes(config["bfile"])

# ---------------------------------------------------------------------------
# Parse-time feature flags
# ---------------------------------------------------------------------------

_use_ref   = bool(config.get("reference_vcf", "").strip())
_use_bref3 = _use_ref and bool(config.get("bref3_jar", "").strip())

# When a reference is used, Beagle writes to imputed_ref/ (markers at ref
# positions only). The merge rule below then adds back target-only markers
# and writes the final result to imputed/ — the same path concat_chromosomes
# expects in both modes, so no downstream rule changes are needed.
_beagle_subdir = "imputed_ref" if _use_ref else "imputed"

_use_alphaimpute2 = config.get("imputer", "beagle") == "alphaimpute2"

# Always run plink2 with --dog so non-human chromosome codes (1–38) are
# accepted. Salmon (29), trout (32), and most livestock fit inside this range,
# and human data still works because plink2's only check is that codes ≤ 38.
# This is prepended to any user-supplied plink_extra_flags so a single source
# of truth flows to every rule via config["plink_extra_flags"].
config["plink_extra_flags"] = ("--dog " + config.get("plink_extra_flags", "")).strip()

# ---------------------------------------------------------------------------
# Auto-download helpers
# If the configured jar path starts with "bin/", treat it as auto-managed and
# add it as an explicit input dependency so Snakemake downloads it first.
# Users with existing JARs can point config to any absolute path and these
# download rules are never triggered.
# ---------------------------------------------------------------------------

_BEAGLE_URL    = "https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar"
_CONFORM_URL   = "https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar"

_beagle_jar  = config.get("beagle_jar",    "bin/beagle.jar")
_conform_jar = config.get("conform_gt_jar", "bin/conform-gt.jar")


def _auto_download(path, url):
    """Download *url* to *path* if path starts with 'bin/' and does not exist."""
    import os, subprocess as sp
    if path.startswith("bin/") and not os.path.exists(path):
        os.makedirs("bin", exist_ok=True)
        print(f"Downloading {url} → {path}")
        sp.run(["wget", "-q", "-O", path, url], check=True)


onstart:
    if not _use_alphaimpute2:
        _auto_download(_beagle_jar, _BEAGLE_URL)
    if _use_ref:
        _auto_download(_conform_jar, _CONFORM_URL)

# ---------------------------------------------------------------------------
# Shared helpers: resolve final imputed VCF path (used by vcf_to_plink)
# ---------------------------------------------------------------------------

# AlphaImpute2 produces a single genome-wide VCF; Beagle produces per-chrom
# VCFs that are concatenated. Both land in the same variable so vcf_to_plink
# needs no conditional logic.
_final_vcf = (
    config["output_dir"] + "/alphaimpute2/all_chromosomes.vcf.gz"
    if _use_alphaimpute2
    else config["output_dir"] + "/all_chromosomes.vcf.gz"
)

# rule all targets differ between imputers
_rule_all_inputs = (
    [
        _final_vcf,
        _final_vcf + ".tbi",
        config["output_dir"] + "/plink_binary/imputed_data.bed",
    ]
    if _use_alphaimpute2
    else
    expand(
        "{output_dir}/imputed/chr{chrom}.vcf.gz",
        output_dir=config["output_dir"],
        chrom=get_chroms(),
    )
    + [
        config["output_dir"] + "/all_chromosomes.vcf.gz",
        config["output_dir"] + "/all_chromosomes.vcf.gz.tbi",
        config["output_dir"] + "/plink_binary/imputed_data.bed",
    ]
)

# ---------------------------------------------------------------------------
# rule all
# ---------------------------------------------------------------------------

rule all:
    input:
        _rule_all_inputs

# ---------------------------------------------------------------------------
# Per-chromosome VCF extraction from PLINK binary
# ---------------------------------------------------------------------------

rule make_per_chrom_vcf:
    input:
        bed = config["bfile"] + ".bed"
    output:
        vcf = temp(config["output_dir"] + "/dedup/chr{chrom}.vcf.gz")
    params:
        plink       = config["plink_path"],
        bfile       = config["bfile"],
        out_prefix  = config["output_dir"] + "/dedup/chr{chrom}",
        chr         = "{chrom}",
        extra_flags = config.get("plink_extra_flags", "")
    conda:
        "envs/workflow_env.yaml"
    log:
        "logs/dedup_chr{chrom}.log"
    threads: 1
    shell:
        """
        ({params.plink} \
            --nonfounders \
            --allow-no-sex \
            --bfile {params.bfile} \
            --chr {params.chr} \
            --export vcf bgz \
            {params.extra_flags} \
            --out {params.out_prefix} \
            --snps-only) &> {log}
        """

# ---------------------------------------------------------------------------
# Deduplication and normalisation
# ---------------------------------------------------------------------------

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
    threads: 1
    shell:
        """
        (bcftools norm -d snps {input.vcf} | \
         bgzip > {output.vcf}) 2> {log}

        tabix -f -p vcf {output.vcf}
        """

# ---------------------------------------------------------------------------
# Beagle imputation
# ---------------------------------------------------------------------------

rule run_beagle:
    input:
        vcf = branch(
            lambda _: _use_ref,
            then      = config["output_dir"] + "/harmonized/chr{chrom}.vcf.gz",
            otherwise = rules.normalize_vcf.output.vcf
        ),
        # bref3 is a per-chromosome binary reference; only present when bref3_jar is set
        bref3 = branch(
            lambda _: _use_bref3,
            then      = config["output_dir"] + "/bref3/chr{chrom}_ref.bref3",
            otherwise = []
        )
    output:
        vcf = temp(config["output_dir"] + "/" + _beagle_subdir + "/chr{chrom}.vcf.gz"),
        tbi = temp(config["output_dir"] + "/" + _beagle_subdir + "/chr{chrom}.vcf.gz.tbi")
    params:
        beagle    = config["beagle_jar"],
        window    = config["beagle_params"]["window"],
        overlap   = config["beagle_params"]["overlap"],
        ne        = config["beagle_params"]["ne"],
        outbase   = lambda wildcards: f"{config['output_dir']}/{_beagle_subdir}/chr{wildcards.chrom}",
        ref_param = (
            lambda wildcards, input:
                f"ref={input.bref3}" if _use_bref3
                else f"ref={config['reference_vcf']}" if _use_ref
                else ""
        )
    threads: config["beagle_params"]["nthreads"]
    conda:
        "envs/workflow_env.yaml"
    log:
        "logs/beagle_chr{chrom}.log"
    resources:
        mem_mb = 70000,
        slurm_partition = "r6i-ondemand-4xlarge"
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

# ---------------------------------------------------------------------------
# Merge imputed markers with target-only markers (ref mode only)
#
# Beagle with ref= outputs only markers present in the reference panel.
# This rule merges those back with chip markers that were not in the
# reference (0000.vcf.gz from bcftools_isec), so the final imputed/
# VCF contains all original chip markers plus any imputed from the ref.
# ---------------------------------------------------------------------------

if _use_ref:
    rule merge_imputed_with_target_only:
        input:
            imputed         = config["output_dir"] + "/" + _beagle_subdir + "/chr{chrom}.vcf.gz",
            imputed_tbi     = config["output_dir"] + "/" + _beagle_subdir + "/chr{chrom}.vcf.gz.tbi",
            target_only     = config["output_dir"] + "/intersect/chr{chrom}_temp/0000.vcf.gz",
            target_only_tbi = config["output_dir"] + "/intersect/chr{chrom}_temp/0000.vcf.gz.tbi"
        output:
            vcf = temp(config["output_dir"] + "/imputed/chr{chrom}.vcf.gz"),
            tbi = temp(config["output_dir"] + "/imputed/chr{chrom}.vcf.gz.tbi")
        conda:
            "envs/workflow_env.yaml"
        log:
            "logs/merge_imputed_chr{chrom}.log"
        resources:
            mem_mb = 16000,
            slurm_partition = "r6i-ondemand-2xlarge"
        shell:
            """
            (bcftools concat \
                --allow-overlaps \
                -O u \
                {input.imputed} \
                {input.target_only} \
            | bcftools sort \
                -O z \
                -m {resources.mem_mb}M \
                -T $(dirname {output.vcf})/sort_tmp_chr{wildcards.chrom} \
                -o {output.vcf}
            ) 2> {log}

            tabix -f {output.vcf}
            """

# ---------------------------------------------------------------------------
# Concatenate all per-chromosome results
# ---------------------------------------------------------------------------

rule concat_chromosomes:
    input:
        vcfs = expand(
            "{output_dir}/imputed/chr{chrom}.vcf.gz",
            output_dir=config["output_dir"],
            chrom=get_chroms()
        ),
        tbi = expand(
            "{output_dir}/imputed/chr{chrom}.vcf.gz.tbi",
            output_dir=config["output_dir"],
            chrom=get_chroms()
        )
    output:
        vcf = config["output_dir"] + "/all_chromosomes.vcf.gz",
        tbi = config["output_dir"] + "/all_chromosomes.vcf.gz.tbi"
    conda:
        "envs/workflow_env.yaml"
    threads: 4
    resources:
        mem_mb = 64000,
        slurm_partition = "r6i-ondemand-4xlarge"
    log:
        "logs/concat_chromosomes.log"
    shell:
        """
        (bcftools concat \
            --output {output.vcf} \
            --output-type z \
            --threads {threads} \
            {input.vcfs}) 2> {log}

        tabix -f {output.vcf}
        """

# ---------------------------------------------------------------------------
# Final VCF to PLINK binary
# ---------------------------------------------------------------------------

rule vcf_to_plink:
    input:
        vcf = _final_vcf
    output:
        bed = config["output_dir"] + "/plink_binary/imputed_data.bed",
        bim = config["output_dir"] + "/plink_binary/imputed_data.bim",
        fam = config["output_dir"] + "/plink_binary/imputed_data.fam"
    params:
        plink       = config["plink_path"],
        out_prefix  = lambda wildcards, output: output.bed.rsplit('.', 1)[0],
        extra_flags = config.get("plink_extra_flags", "")
    conda:
        "envs/workflow_env.yaml"
    log:
        "logs/vcf_to_plink.log"
    shell:
        """
        ({params.plink} \
            --vcf {input.vcf} \
            --make-bed \
            --out {params.out_prefix} \
            {params.extra_flags} \
            --allow-no-sex \
            --const-fid) &> {log}
        """

# ---------------------------------------------------------------------------
# Conditional includes (placed after all rules so rules.xxx references resolve)
# ---------------------------------------------------------------------------

if _use_ref:
    include: "rules/intersect_and_conform.smk"

if _use_alphaimpute2:
    include: "rules/alphaimpute2.smk"
