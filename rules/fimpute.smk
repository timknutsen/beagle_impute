# rules/fimpute.smk
#
# FImpute3 pedigree + population imputation pipeline.
# Activated when config["imputer"] == "fimpute".
#
# Rule graph (per chromosome, then concatenated):
#   fimpute_chrom_bfile → fimpute_export_raw → fimpute_prepare_inputs
#     → run_fimpute → fimpute_to_vcf → concat_fimpute
#
# The final output ({output_dir}/fimpute/all_chromosomes.vcf.gz) is consumed by
# the shared vcf_to_plink rule in the main Snakefile.
#
# FImpute is a licensed binary that is NOT installed by conda — it must already
# exist on the machine and is located via config["fimpute_params"]["executable"].
# Everything else here runs inside envs/workflow_env.yaml.

_fim_cfg          = config.get("fimpute_params", {})
_fimpute_exe      = _fim_cfg.get(
    "executable", "/mnt/efshome/applications/FImpute3/2026/FImpute3"
)
_fimpute_nthreads = int(_fim_cfg.get("nthreads", 1))
_fim_dir          = config["output_dir"] + "/fimpute_work"

# Optional denser (or simply different) panel used as chip 1. Empty means the
# target is imputed on its own, which still fills sporadic gaps and phases.
_fim_ref_bfile = str(_fim_cfg.get("reference_bfile", "")).strip()
_fim_use_ref   = bool(_fim_ref_bfile)
_fim_keep_og   = bool(_fim_cfg.get("keep_observed", True))
_fim_phased    = bool(_fim_cfg.get("phase_output", True))


# ---------------------------------------------------------------------------
# Step 1: One PLINK binary per chromosome
#
# FImpute treats each run as a single linkage group, so the genome is split
# first and every chromosome is imputed independently.
# ---------------------------------------------------------------------------

rule fimpute_chrom_bfile:
    input:
        bed = config["bfile"] + ".bed",
    output:
        bed = temp(_fim_dir + "/chrom/chr{chrom}.bed"),
        bim = temp(_fim_dir + "/chrom/chr{chrom}.bim"),
        fam = temp(_fim_dir + "/chrom/chr{chrom}.fam"),
    params:
        plink = config["plink_path"],
        bfile = config["bfile"],
        out   = _fim_dir + "/chrom/chr{chrom}",
        extra = config.get("plink_extra_flags", ""),
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/fimpute_chrom_chr{chrom}.log",
    shell:
        """
        ({params.plink} \
            --bfile {params.bfile} \
            --chr {wildcards.chrom} \
            --make-bed \
            --out {params.out} \
            --nonfounders --allow-no-sex {params.extra}) &> {log}
        """


# ---------------------------------------------------------------------------
# Step 2: Additive-coded export (0/1/2 copies of A1, NA = missing)
# ---------------------------------------------------------------------------

rule fimpute_export_raw:
    input:
        bed = rules.fimpute_chrom_bfile.output.bed,
    output:
        raw = temp(_fim_dir + "/input/chr{chrom}.raw"),
    params:
        plink = config["plink_path"],
        bfile = _fim_dir + "/chrom/chr{chrom}",
        out   = _fim_dir + "/input/chr{chrom}",
        extra = config.get("plink_extra_flags", ""),
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/fimpute_export_raw_chr{chrom}.log",
    shell:
        """
        ({params.plink} \
            --bfile {params.bfile} \
            --export A \
            {params.extra} --allow-no-sex \
            --out {params.out}) &> {log}
        """


# ---------------------------------------------------------------------------
# Step 3: FImpute's own text inputs
#
# Animals are renumbered 1..N on the way in because FImpute constrains ID
# length; chr{chrom}.id_map.tsv carries the original IIDs back out again.
# ---------------------------------------------------------------------------

if _fim_use_ref:

    rule fimpute_ref_chrom_bfile:
        """Reference panel restricted to one chromosome (chip 1)."""
        input:
            bed = _fim_ref_bfile + ".bed",
        output:
            bed = temp(_fim_dir + "/ref/chr{chrom}.bed"),
            bim = temp(_fim_dir + "/ref/chr{chrom}.bim"),
            fam = temp(_fim_dir + "/ref/chr{chrom}.fam"),
        params:
            plink = config["plink_path"],
            bfile = _fim_ref_bfile,
            out   = _fim_dir + "/ref/chr{chrom}",
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/fimpute_ref_chrom_chr{chrom}.log",
        shell:
            """
            ({params.plink} \
                --bfile {params.bfile} \
                --chr {wildcards.chrom} \
                --make-bed \
                --out {params.out} \
                --nonfounders --allow-no-sex {params.extra}) &> {log}
            """

    rule fimpute_ref_export_raw:
        input:
            bed = rules.fimpute_ref_chrom_bfile.output.bed,
        output:
            raw = temp(_fim_dir + "/ref/chr{chrom}.raw"),
        params:
            plink = config["plink_path"],
            bfile = _fim_dir + "/ref/chr{chrom}",
            out   = _fim_dir + "/ref/chr{chrom}",
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/fimpute_ref_export_raw_chr{chrom}.log",
        shell:
            """
            ({params.plink} \
                --bfile {params.bfile} \
                --export A \
                {params.extra} --allow-no-sex \
                --out {params.out}) &> {log}
            """


rule fimpute_prepare_inputs:
    input:
        raw = rules.fimpute_export_raw.output.raw,
        bim = rules.fimpute_chrom_bfile.output.bim,
        fam = rules.fimpute_chrom_bfile.output.fam,
        ref = (
            [
                _fim_dir + "/ref/chr{chrom}.raw",
                _fim_dir + "/ref/chr{chrom}.bim",
                _fim_dir + "/ref/chr{chrom}.fam",
            ]
            if _fim_use_ref
            else []
        ),
    output:
        genos = _fim_dir + "/input/chr{chrom}.genos",
        snps  = _fim_dir + "/input/chr{chrom}.snps",
        ped   = _fim_dir + "/input/chr{chrom}.ped",
        ctrl  = _fim_dir + "/input/chr{chrom}.ctrl",
        idmap = _fim_dir + "/input/chr{chrom}.id_map.tsv",
    params:
        out_dir  = _fim_dir + "/input",
        nthreads = _fimpute_nthreads,
        ref_args = (
            lambda wc, input: (
                f"--ref-raw {input.ref[0]} --ref-bim {input.ref[1]} --ref-fam {input.ref[2]}"
                if _fim_use_ref
                else ""
            )
        ),
        keep_og  = "" if _fim_keep_og else "--no-keep-observed",
        phase    = "" if _fim_phased else "--no-phase",
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/fimpute_prepare_chr{chrom}.log",
    shell:
        """
        python {workflow.basedir}/scripts/fimpute_io.py prepare-inputs \
            --raw {input.raw} \
            --bim {input.bim} \
            --fam {input.fam} \
            --out-dir {params.out_dir} \
            --chrom {wildcards.chrom} \
            --nthreads {params.nthreads} \
            {params.ref_args} {params.keep_og} {params.phase} \
            &> {log}
        """


# ---------------------------------------------------------------------------
# Step 4: Run FImpute
#
# -o overwrites the output folder, which makes the rule re-runnable.
# ---------------------------------------------------------------------------

rule run_fimpute:
    input:
        ctrl = rules.fimpute_prepare_inputs.output.ctrl,
    output:
        imp    = _fim_dir + "/input/chr{chrom}/genotypes_imp.txt",
        report = _fim_dir + "/input/chr{chrom}/report.txt",
        # FImpute rewrites the SNP map for the markers it actually kept. It
        # drops some of its own accord -- a marker carried only by the target
        # chip is logged in excluded_snp_list.txt as "Not On HD" -- so the
        # input .snps overstates the genotype string by however many it
        # removed, and the VCF writer must follow this file instead.
        snpinfo = _fim_dir + "/input/chr{chrom}/snp_info.txt",
    params:
        exe = _fimpute_exe,
    threads:
        _fimpute_nthreads
    log:
        "logs/fimpute_run_chr{chrom}.log",
    resources:
        mem_mb = 32000,
        slurm_partition = "r7i-ondemand-2xlarge",
    shell:
        """
        ({params.exe} {input.ctrl} -o) &> {log}
        """


# ---------------------------------------------------------------------------
# Step 5: Back to VCF, using the .bim for position/allele metadata
# ---------------------------------------------------------------------------

rule fimpute_to_vcf:
    input:
        imp   = rules.run_fimpute.output.imp,
        bim   = rules.fimpute_chrom_bfile.output.bim,
        # FImpute's own post-run map, not the one we handed it -- see the note
        # on run_fimpute.output.snpinfo.
        snps  = rules.run_fimpute.output.snpinfo,
        idmap = rules.fimpute_prepare_inputs.output.idmap,
        refbim = (
            _fim_dir + "/ref/chr{chrom}.bim" if _fim_use_ref else []
        ),
    output:
        vcf = temp(_fim_dir + "/imputed/chr{chrom}.vcf.gz"),
        tbi = temp(_fim_dir + "/imputed/chr{chrom}.vcf.gz.tbi"),
    params:
        # Marker order follows the SNP map, which is the union across chips.
        # input.refbim is a bare string when there is one reference .bim and
        # empty when there is none, so normalise before joining -- list() on a
        # string would splice it into single characters.
        bims     = lambda wc, input: " ".join(
            [input.bim]
            + ([input.refbim] if isinstance(input.refbim, str) else list(input.refbim))
        ),
        unphased = "" if _fim_phased else "--unphased",
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/fimpute_to_vcf_chr{chrom}.log",
    shell:
        """
        (python {workflow.basedir}/scripts/fimpute_io.py to-vcf \
            --imputed {input.imp} \
            --bim {params.bims} \
            --snp-info {input.snps} \
            --id-map {input.idmap} \
            {params.unphased} \
            --out-vcf {output.vcf}.tmp

        bgzip -c {output.vcf}.tmp > {output.vcf}
        rm -f {output.vcf}.tmp
        tabix -f -p vcf {output.vcf}) &> {log}
        """


# ---------------------------------------------------------------------------
# Step 6: Concatenate the per-chromosome results
# ---------------------------------------------------------------------------

rule concat_fimpute:
    input:
        vcfs = expand(
            _fim_dir + "/imputed/chr{chrom}.vcf.gz", chrom=get_chroms()
        ),
        tbis = expand(
            _fim_dir + "/imputed/chr{chrom}.vcf.gz.tbi", chrom=get_chroms()
        ),
    output:
        vcf = config["output_dir"] + "/fimpute/all_chromosomes.vcf.gz",
        tbi = config["output_dir"] + "/fimpute/all_chromosomes.vcf.gz.tbi",
    conda:
        "../envs/workflow_env.yaml"
    threads: 4
    log:
        "logs/fimpute_concat.log",
    resources:
        mem_mb = 64000,
        slurm_partition = "r7i-ondemand-4xlarge",
    shell:
        """
        (bcftools concat \
            --output {output.vcf} \
            --output-type z \
            --threads {threads} \
            {input.vcfs}) 2> {log}

        tabix -f {output.vcf}
        """
