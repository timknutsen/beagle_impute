# rules/accuracy.smk — imputation accuracy evaluation rules
#
# Included exclusively by Snakefile_accuracy; never by the main Snakefile.
#
# Structure:
#   1. Mode-conditional data preparation
#      mask_and_impute  → acc_sample_validation_ids, acc_sample_snp_list,
#                         acc_make_truth_bfile, acc_make_masked_ld_bfile,
#                         (+ acc_make_reference_panel, acc_merge_for_alphaimpute2
#                            when imputer=alphaimpute2)
#      cross_array      → acc_find_common_animals, acc_make_ld_subset,
#                         acc_make_hd_truth
#
#   2. Shared imputation rules
#      Beagle path      → acc_make_per_chrom_vcf, acc_normalize_vcf,
#                         acc_run_beagle, acc_concat_imputed
#      AlphaImpute2     → acc_plink_to_alphaimpute2_fmt, acc_run_alphaimpute2,
#                         acc_alphaimpute2_to_vcf
#
#   3. Truth conversion → acc_truth_to_vcf
#   4. Metrics          → acc_compute_metrics
#
# After data preparation each branch defines two path variables:
#   _impute_bfile   PLINK prefix fed into the imputation rules
#   _truth_bfile    PLINK prefix used to produce the truth VCF

_acc_mode         = config.get("accuracy_mode", "mask_and_impute")
_acc_out          = config.get("accuracy_output_dir", "accuracy_output")
_use_alphaimpute2 = config.get("imputer", "beagle") == "alphaimpute2"


def _nested_config(section, key, default=None):
    """Support both YAML nested config and Snakemake --config section.key=value."""
    if isinstance(config.get(section), dict) and key in config[section]:
        return config[section][key]
    if f"{section}_{key}" in config:
        return config[f"{section}_{key}"]
    return config.get(f"{section}.{key}", default)


def _as_list(value, default):
    if value is None:
        value = default
    if isinstance(value, str):
        return [item for item in value.replace(",", " ").split() if item]
    return list(value)


_cv_n_folds = int(_nested_config("cv", "n_folds", 10))
_cv_folds = list(range(1, _cv_n_folds + 1))
_cv_imputers = _as_list(_nested_config("cv", "imputers", ["beagle", "alphaimpute2", "fimpute"]), [])
_unsupported_cv_imputers = sorted(set(_cv_imputers) - {"beagle", "alphaimpute2", "fimpute"})
if _unsupported_cv_imputers:
    raise ValueError(
        "Unsupported cv imputers: "
        + ", ".join(_unsupported_cv_imputers)
        + ". Use beagle, alphaimpute2, and/or fimpute."
    )
_cv_target_n_snps = int(_nested_config("cv", "target_n_snps", 10000))
_cv_seed = int(_nested_config("cv", "random_seed", 42))
_cv_target_snp_list = _nested_config("cv", "target_snp_list", "")
_fimpute_executable = _nested_config(
    "fimpute_params",
    "executable",
    "/mnt/efshome/applications/FImpute3/2026/FImpute3",
)
_fimpute_nthreads = int(_nested_config("fimpute_params", "nthreads", 1))


def _cv_prefix(imputer, fold, stem):
    return f"{_acc_out}/{imputer}/fold{fold}/{stem}"


def _cv_imputed_vcf(wc):
    if wc.imputer == "beagle":
        return _cv_prefix(wc.imputer, wc.fold, "imputed/all_chromosomes.vcf.gz")
    if wc.imputer == "alphaimpute2":
        return _cv_prefix(wc.imputer, wc.fold, "alphaimpute2/all_chromosomes.vcf.gz")
    if wc.imputer == "fimpute":
        return _cv_prefix(wc.imputer, wc.fold, "fimpute/all_chromosomes.vcf.gz")
    raise ValueError(f"Unsupported CV imputer: {wc.imputer}")


def _cv_truth_vcf(wc):
    return _cv_prefix(wc.imputer, wc.fold, "truth/all_chromosomes.vcf.gz")


# ---------------------------------------------------------------------------
# Mode: kfold_mask_and_impute
# ---------------------------------------------------------------------------

if _acc_mode == "kfold_mask_and_impute":

    rule acc_cv_setup:
        """Create the shared fold assignment and LD marker panel."""
        input:
            fam = config["bfile"] + ".fam",
            bim = config["bfile"] + ".bim",
        output:
            folds = _acc_out + "/setup/folds.tsv",
            snps  = _acc_out + "/setup/ld_snp_list.txt",
        params:
            n_folds = _cv_n_folds,
            n_snps  = _cv_target_n_snps,
            seed    = _cv_seed,
            snp_list = _cv_target_snp_list,
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/setup.log",
        shell:
            """
            python {workflow.basedir}/scripts/make_accuracy_cv_setup.py \
                --fam {input.fam} \
                --bim {input.bim} \
                --n-folds {params.n_folds} \
                --target-n-snps {params.n_snps} \
                --seed {params.seed} \
                --target-snp-list "{params.snp_list}" \
                --folds-out {output.folds} \
                --snps-out {output.snps} \
                &> {log}
            """

    rule acc_cv_fold_ids:
        input:
            folds = _acc_out + "/setup/folds.tsv",
        output:
            ids = _acc_out + "/setup/fold{fold}_ids.txt",
        run:
            import pandas as pd

            folds = pd.read_csv(input.folds, sep="\t")
            fold = int(wildcards.fold)
            folds.loc[folds["fold"] == fold, ["fid", "iid"]].to_csv(
                output.ids, sep="\t", header=False, index=False
            )

    rule acc_cv_make_truth_bfile:
        input:
            bed = config["bfile"] + ".bed",
            ids = _acc_out + "/setup/fold{fold}_ids.txt",
        output:
            bed = _acc_out + "/{imputer}/fold{fold}/truth/hd.bed",
            bim = _acc_out + "/{imputer}/fold{fold}/truth/hd.bim",
            fam = _acc_out + "/{imputer}/fold{fold}/truth/hd.fam",
        params:
            plink = config["plink_path"],
            bfile = config["bfile"],
            out   = lambda wc: _cv_prefix(wc.imputer, wc.fold, "truth/hd"),
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/{imputer}_fold{fold}_truth.log",
        shell:
            """
            ({params.plink} \
                --bfile {params.bfile} \
                --keep {input.ids} \
                --make-bed \
                --out {params.out} \
                --nonfounders --allow-no-sex {params.extra}) &> {log}
            """

    rule acc_cv_make_masked_ld_bfile:
        input:
            bed  = config["bfile"] + ".bed",
            ids  = _acc_out + "/setup/fold{fold}_ids.txt",
            snps = _acc_out + "/setup/ld_snp_list.txt",
        output:
            bed = _acc_out + "/{imputer}/fold{fold}/to_impute/masked.bed",
            bim = _acc_out + "/{imputer}/fold{fold}/to_impute/masked.bim",
            fam = _acc_out + "/{imputer}/fold{fold}/to_impute/masked.fam",
        params:
            plink = config["plink_path"],
            bfile = config["bfile"],
            out   = lambda wc: _cv_prefix(wc.imputer, wc.fold, "to_impute/masked"),
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/{imputer}_fold{fold}_masked_ld.log",
        shell:
            """
            ({params.plink} \
                --bfile {params.bfile} \
                --keep {input.ids} \
                --extract {input.snps} \
                --make-bed \
                --out {params.out} \
                --nonfounders --allow-no-sex {params.extra}) &> {log}
            """

    rule acc_cv_make_reference_bfile:
        input:
            bed = config["bfile"] + ".bed",
            ids = _acc_out + "/setup/fold{fold}_ids.txt",
        output:
            bed = _acc_out + "/{imputer}/fold{fold}/reference/panel.bed",
            bim = _acc_out + "/{imputer}/fold{fold}/reference/panel.bim",
            fam = _acc_out + "/{imputer}/fold{fold}/reference/panel.fam",
        params:
            plink = config["plink_path"],
            bfile = config["bfile"],
            out   = lambda wc: _cv_prefix(wc.imputer, wc.fold, "reference/panel"),
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/{imputer}_fold{fold}_reference.log",
        shell:
            """
            ({params.plink} \
                --bfile {params.bfile} \
                --remove {input.ids} \
                --make-bed \
                --out {params.out} \
                --nonfounders --allow-no-sex {params.extra}) &> {log}
            """

    rule acc_cv_make_combined_masked_bfile:
        input:
            bed  = config["bfile"] + ".bed",
            bim  = config["bfile"] + ".bim",
            fam  = config["bfile"] + ".fam",
            ids  = _acc_out + "/setup/fold{fold}_ids.txt",
            snps = _acc_out + "/setup/ld_snp_list.txt",
        output:
            bed = _acc_out + "/{imputer}/fold{fold}/to_impute/combined.bed",
            bim = _acc_out + "/{imputer}/fold{fold}/to_impute/combined.bim",
            fam = _acc_out + "/{imputer}/fold{fold}/to_impute/combined.fam",
        params:
            bfile = config["bfile"],
            out   = lambda wc: _cv_prefix(wc.imputer, wc.fold, "to_impute/combined"),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/{imputer}_fold{fold}_combined_masked.log",
        shell:
            """
            python {workflow.basedir}/scripts/mask_validation_genotypes.py \
                --bfile {params.bfile} \
                --validation-ids {input.ids} \
                --ld-snps {input.snps} \
                --out {params.out} \
                &> {log}
            """

    rule acc_cv_truth_to_vcf:
        input:
            bed = _acc_out + "/{imputer}/fold{fold}/truth/hd.bed",
        output:
            vcf = _acc_out + "/{imputer}/fold{fold}/truth/all_chromosomes.vcf.gz",
            tbi = _acc_out + "/{imputer}/fold{fold}/truth/all_chromosomes.vcf.gz.tbi",
        params:
            plink = config["plink_path"],
            bfile = lambda wc: _cv_prefix(wc.imputer, wc.fold, "truth/hd"),
            out   = lambda wc: _cv_prefix(wc.imputer, wc.fold, "truth/all_chromosomes"),
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/{imputer}_fold{fold}_truth_to_vcf.log",
        shell:
            """
            ({params.plink} \
                --bfile {params.bfile} \
                --export vcf bgz \
                --nonfounders --allow-no-sex \
                {params.extra} \
                --out {params.out} \
                --snps-only) &> {log}
            tabix -f -p vcf {output.vcf}
            """

    # Beagle CV path: held-out LD animals as gt=, remaining full-density animals as ref=.

    rule acc_cv_beagle_gt_vcf:
        input:
            bed = _acc_out + "/beagle/fold{fold}/to_impute/masked.bed",
        output:
            vcf = temp(_acc_out + "/beagle/fold{fold}/dedup/gt_chr{chrom}.vcf.gz"),
        params:
            plink = config["plink_path"],
            bfile = lambda wc: _cv_prefix("beagle", wc.fold, "to_impute/masked"),
            out   = lambda wc: _cv_prefix("beagle", wc.fold, f"dedup/gt_chr{wc.chrom}"),
            chr   = "{chrom}",
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/beagle_fold{fold}_gt_chr{chrom}.log",
        shell:
            """
            ({params.plink} --bfile {params.bfile} --chr {params.chr} \
                --export vcf bgz {params.extra} --out {params.out} \
                --nonfounders --allow-no-sex --snps-only) &> {log}
            """

    rule acc_cv_beagle_ref_vcf:
        input:
            bed = _acc_out + "/beagle/fold{fold}/reference/panel.bed",
        output:
            vcf = temp(_acc_out + "/beagle/fold{fold}/dedup/ref_chr{chrom}.vcf.gz"),
        params:
            plink = config["plink_path"],
            bfile = lambda wc: _cv_prefix("beagle", wc.fold, "reference/panel"),
            out   = lambda wc: _cv_prefix("beagle", wc.fold, f"dedup/ref_chr{wc.chrom}"),
            chr   = "{chrom}",
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/beagle_fold{fold}_ref_chr{chrom}.log",
        shell:
            """
            ({params.plink} --bfile {params.bfile} --chr {params.chr} \
                --export vcf bgz {params.extra} --out {params.out} \
                --nonfounders --allow-no-sex --snps-only) &> {log}
            """

    rule acc_cv_beagle_normalize_gt:
        input:
            vcf = _acc_out + "/beagle/fold{fold}/dedup/gt_chr{chrom}.vcf.gz",
        output:
            vcf = temp(_acc_out + "/beagle/fold{fold}/normalized/gt_chr{chrom}.vcf.gz"),
            tbi = temp(_acc_out + "/beagle/fold{fold}/normalized/gt_chr{chrom}.vcf.gz.tbi"),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/beagle_fold{fold}_normalize_gt_chr{chrom}.log",
        shell:
            """
            (bcftools norm -d snps {input.vcf} | bgzip > {output.vcf}) 2> {log}
            tabix -f -p vcf {output.vcf}
            """

    rule acc_cv_beagle_normalize_ref:
        input:
            vcf = _acc_out + "/beagle/fold{fold}/dedup/ref_chr{chrom}.vcf.gz",
        output:
            vcf = temp(_acc_out + "/beagle/fold{fold}/normalized/ref_chr{chrom}.vcf.gz"),
            tbi = temp(_acc_out + "/beagle/fold{fold}/normalized/ref_chr{chrom}.vcf.gz.tbi"),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/beagle_fold{fold}_normalize_ref_chr{chrom}.log",
        shell:
            """
            (bcftools norm -d snps {input.vcf} | bgzip > {output.vcf}) 2> {log}
            tabix -f -p vcf {output.vcf}
            """

    rule acc_cv_run_beagle:
        input:
            gt  = _acc_out + "/beagle/fold{fold}/normalized/gt_chr{chrom}.vcf.gz",
            ref = _acc_out + "/beagle/fold{fold}/normalized/ref_chr{chrom}.vcf.gz",
        output:
            vcf = temp(_acc_out + "/beagle/fold{fold}/imputed/chr{chrom}.vcf.gz"),
            tbi = temp(_acc_out + "/beagle/fold{fold}/imputed/chr{chrom}.vcf.gz.tbi"),
        params:
            beagle  = config["beagle_jar"],
            window  = config["beagle_params"]["window"],
            overlap = config["beagle_params"]["overlap"],
            ne      = config["beagle_params"]["ne"],
            outbase = lambda wc: _cv_prefix("beagle", wc.fold, f"imputed/chr{wc.chrom}"),
        threads:
            config["beagle_params"]["nthreads"]
        conda:
            "../envs/workflow_env.yaml"
        resources:
            mem_mb = 70000,
            slurm_partition = "r6i-ondemand-4xlarge",
        log:
            "logs/accuracy_cv/beagle_fold{fold}_chr{chrom}.log",
        shell:
            """
            (java -Xmx{resources.mem_mb}m -jar {params.beagle} \
                gt={input.gt} ref={input.ref} \
                window={params.window} overlap={params.overlap} \
                out={params.outbase} nthreads={threads} ne={params.ne} \
                chrom={wildcards.chrom}) &> {log}
            tabix -f {output.vcf}
            """

    rule acc_cv_concat_beagle:
        input:
            vcfs = lambda wc: expand(
                _acc_out + "/beagle/fold{fold}/imputed/chr{chrom}.vcf.gz",
                fold=wc.fold,
                chrom=acc_get_chromosomes(),
            ),
        output:
            vcf = _acc_out + "/beagle/fold{fold}/imputed/all_chromosomes.vcf.gz",
            tbi = _acc_out + "/beagle/fold{fold}/imputed/all_chromosomes.vcf.gz.tbi",
        conda:
            "../envs/workflow_env.yaml"
        threads: 4
        resources:
            mem_mb = 64000,
            slurm_partition = "r6i-ondemand-4xlarge",
        log:
            "logs/accuracy_cv/beagle_fold{fold}_concat.log",
        shell:
            """
            (bcftools concat --output {output.vcf} --output-type z \
                --threads {threads} {input.vcfs}) 2> {log}
            tabix -f {output.vcf}
            """

    # AlphaImpute2 CV path: combined full-density reference + LD held-out animals.

    rule acc_cv_make_alphaimpute2_chrom_bfile:
        input:
            bed = _acc_out + "/alphaimpute2/fold{fold}/to_impute/combined.bed",
            bim = _acc_out + "/alphaimpute2/fold{fold}/to_impute/combined.bim",
            fam = _acc_out + "/alphaimpute2/fold{fold}/to_impute/combined.fam",
        output:
            bed = temp(_acc_out + "/alphaimpute2/fold{fold}/chrom/chr{chrom}.bed"),
            bim = temp(_acc_out + "/alphaimpute2/fold{fold}/chrom/chr{chrom}.bim"),
            fam = temp(_acc_out + "/alphaimpute2/fold{fold}/chrom/chr{chrom}.fam"),
        params:
            plink = config["plink_path"],
            bfile = lambda wc: _cv_prefix("alphaimpute2", wc.fold, "to_impute/combined"),
            out   = lambda wc: _cv_prefix("alphaimpute2", wc.fold, f"chrom/chr{wc.chrom}"),
            chr   = "{chrom}",
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/alphaimpute2_fold{fold}_chrom_chr{chrom}.log",
        shell:
            """
            ({params.plink} --bfile {params.bfile} --chr {params.chr} \
                --make-bed --out {params.out} \
                --nonfounders --allow-no-sex {params.extra}) &> {log}
            """

    rule acc_cv_plink_to_alphaimpute2_fmt:
        input:
            bed = _acc_out + "/alphaimpute2/fold{fold}/chrom/chr{chrom}.bed",
            bim = _acc_out + "/alphaimpute2/fold{fold}/chrom/chr{chrom}.bim",
            fam = _acc_out + "/alphaimpute2/fold{fold}/chrom/chr{chrom}.fam",
        output:
            genotypes = _acc_out + "/alphaimpute2/fold{fold}/alphaimpute2_input/chr{chrom}.genotypes.txt",
            pedigree  = _acc_out + "/alphaimpute2/fold{fold}/alphaimpute2_input/chr{chrom}.pedigree.txt",
        params:
            plink      = config["plink_path"],
            bfile      = lambda wc: _cv_prefix("alphaimpute2", wc.fold, f"chrom/chr{wc.chrom}"),
            raw_prefix = lambda wc: _cv_prefix("alphaimpute2", wc.fold, f"alphaimpute2_input/chr{wc.chrom}.raw_tmp"),
            extra      = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/alphaimpute2_fold{fold}_to_fmt_chr{chrom}.log",
        shell:
            """
            (
            mkdir -p "$(dirname {output.genotypes})"
            {params.plink} --bfile {params.bfile} --export A \
                {params.extra} --allow-no-sex --out {params.raw_prefix}
            awk 'NR>1 {{
                printf $2;
                for (i=7; i<=NF; i++) {{
                    if ($i == "NA") printf " 9";
                    else printf " %d", $i;
                }}
                printf "\\n";
            }}' {params.raw_prefix}.raw > {output.genotypes}
            awk '{{print $2, $3, $4}}' {input.fam} > {output.pedigree}
            ) &> {log}
            """

    rule acc_cv_run_alphaimpute2:
        input:
            genotypes = _acc_out + "/alphaimpute2/fold{fold}/alphaimpute2_input/chr{chrom}.genotypes.txt",
            pedigree  = _acc_out + "/alphaimpute2/fold{fold}/alphaimpute2_input/chr{chrom}.pedigree.txt",
        output:
            genotypes = temp(_acc_out + "/alphaimpute2/fold{fold}/alphaimpute2_output/chr{chrom}.genotypes"),
        params:
            out_prefix  = lambda wc: _cv_prefix("alphaimpute2", wc.fold, f"alphaimpute2_output/chr{wc.chrom}"),
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
        resources:
            mem_mb = 16000,
            slurm_partition = "r6i-ondemand-2xlarge",
        log:
            "logs/accuracy_cv/alphaimpute2_fold{fold}_chr{chrom}.log",
        shell:
            """
            mkdir -p "$(dirname {output.genotypes})"
            (AlphaImpute2 \
                -genotypes {input.genotypes} \
                -pedigree {input.pedigree} \
                -out {params.out_prefix} \
                -maxthreads {threads} \
                -cycles {params.cycles} \
                -final_peeling_threshold {params.threshold} \
                -hd_threshold {params.hd_thresh} \
                -length {params.length} \
                {params.extra_flags}) &> {log}
            """

    rule acc_cv_alphaimpute2_to_vcf:
        input:
            genotypes = _acc_out + "/alphaimpute2/fold{fold}/alphaimpute2_output/chr{chrom}.genotypes",
            bim       = _acc_out + "/alphaimpute2/fold{fold}/chrom/chr{chrom}.bim",
        output:
            vcf = temp(_acc_out + "/alphaimpute2/fold{fold}/alphaimpute2/chr{chrom}.vcf.gz"),
            tbi = temp(_acc_out + "/alphaimpute2/fold{fold}/alphaimpute2/chr{chrom}.vcf.gz.tbi"),
        conda:
            "../envs/workflow_env.yaml"
        resources:
            mem_mb = 32000,
            slurm_partition = "r6i-ondemand-2xlarge",
        log:
            "logs/accuracy_cv/alphaimpute2_fold{fold}_to_vcf_chr{chrom}.log",
        script:
            "../scripts/alphaimpute2_to_vcf.py"

    rule acc_cv_concat_alphaimpute2:
        input:
            vcfs = lambda wc: expand(
                _acc_out + "/alphaimpute2/fold{fold}/alphaimpute2/chr{chrom}.vcf.gz",
                fold=wc.fold,
                chrom=acc_get_chromosomes(),
            ),
        output:
            vcf = _acc_out + "/alphaimpute2/fold{fold}/alphaimpute2/all_chromosomes.vcf.gz",
            tbi = _acc_out + "/alphaimpute2/fold{fold}/alphaimpute2/all_chromosomes.vcf.gz.tbi",
        conda:
            "../envs/workflow_env.yaml"
        threads: 4
        resources:
            mem_mb = 64000,
            slurm_partition = "r6i-ondemand-4xlarge",
        log:
            "logs/accuracy_cv/alphaimpute2_fold{fold}_concat.log",
        shell:
            """
            (bcftools concat --output {output.vcf} --output-type z \
                --threads {threads} {input.vcfs}) 2> {log}
            tabix -f {output.vcf}
            """

    # FImpute CV path: combined full-density reference + LD held-out animals.

    rule acc_cv_make_fimpute_chrom_bfile:
        input:
            bed = _acc_out + "/fimpute/fold{fold}/to_impute/combined.bed",
            bim = _acc_out + "/fimpute/fold{fold}/to_impute/combined.bim",
            fam = _acc_out + "/fimpute/fold{fold}/to_impute/combined.fam",
        output:
            bed = temp(_acc_out + "/fimpute/fold{fold}/chrom/chr{chrom}.bed"),
            bim = temp(_acc_out + "/fimpute/fold{fold}/chrom/chr{chrom}.bim"),
            fam = temp(_acc_out + "/fimpute/fold{fold}/chrom/chr{chrom}.fam"),
        params:
            plink = config["plink_path"],
            bfile = lambda wc: _cv_prefix("fimpute", wc.fold, "to_impute/combined"),
            out   = lambda wc: _cv_prefix("fimpute", wc.fold, f"chrom/chr{wc.chrom}"),
            chr   = "{chrom}",
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/fimpute_fold{fold}_chrom_chr{chrom}.log",
        shell:
            """
            ({params.plink} --bfile {params.bfile} --chr {params.chr} \
                --make-bed --out {params.out} \
                --nonfounders --allow-no-sex {params.extra}) &> {log}
            """

    rule acc_cv_export_fimpute_raw:
        input:
            bed = _acc_out + "/fimpute/fold{fold}/chrom/chr{chrom}.bed",
        output:
            raw = _acc_out + "/fimpute/fold{fold}/fimpute_input/chr{chrom}.raw",
        params:
            plink = config["plink_path"],
            bfile = lambda wc: _cv_prefix("fimpute", wc.fold, f"chrom/chr{wc.chrom}"),
            out   = lambda wc: _cv_prefix("fimpute", wc.fold, f"fimpute_input/chr{wc.chrom}"),
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/fimpute_fold{fold}_export_raw_chr{chrom}.log",
        shell:
            """
            ({params.plink} --bfile {params.bfile} --export A \
                {params.extra} --allow-no-sex --out {params.out}) &> {log}
            """

    rule acc_cv_prepare_fimpute_inputs:
        input:
            raw = _acc_out + "/fimpute/fold{fold}/fimpute_input/chr{chrom}.raw",
            bim = _acc_out + "/fimpute/fold{fold}/chrom/chr{chrom}.bim",
            fam = _acc_out + "/fimpute/fold{fold}/chrom/chr{chrom}.fam",
        output:
            genos = _acc_out + "/fimpute/fold{fold}/fimpute_input/chr{chrom}.genos",
            snps  = _acc_out + "/fimpute/fold{fold}/fimpute_input/chr{chrom}.snps",
            ped   = _acc_out + "/fimpute/fold{fold}/fimpute_input/chr{chrom}.ped",
            ctrl  = _acc_out + "/fimpute/fold{fold}/fimpute_input/chr{chrom}.ctrl",
            idmap = _acc_out + "/fimpute/fold{fold}/fimpute_input/chr{chrom}.id_map.tsv",
        params:
            out_dir  = lambda wc: _cv_prefix("fimpute", wc.fold, "fimpute_input"),
            nthreads = _fimpute_nthreads,
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/fimpute_fold{fold}_prepare_chr{chrom}.log",
        shell:
            """
            python {workflow.basedir}/scripts/fimpute_io.py prepare-inputs \
                --raw {input.raw} \
                --bim {input.bim} \
                --fam {input.fam} \
                --out-dir {params.out_dir} \
                --chrom {wildcards.chrom} \
                --nthreads {params.nthreads} \
                &> {log}
            """

    rule acc_cv_run_fimpute:
        input:
            ctrl = _acc_out + "/fimpute/fold{fold}/fimpute_input/chr{chrom}.ctrl",
        output:
            imp = _acc_out + "/fimpute/fold{fold}/fimpute_input/chr{chrom}/genotypes_imp.txt",
            report = _acc_out + "/fimpute/fold{fold}/fimpute_input/chr{chrom}/report.txt",
        params:
            exe = _fimpute_executable,
        threads:
            _fimpute_nthreads
        resources:
            mem_mb = 32000,
            slurm_partition = "r6i-ondemand-2xlarge",
        log:
            "logs/accuracy_cv/fimpute_fold{fold}_run_chr{chrom}.log",
        shell:
            """
            ({params.exe} {input.ctrl} -o) &> {log}
            """

    rule acc_cv_fimpute_to_vcf:
        input:
            imp   = _acc_out + "/fimpute/fold{fold}/fimpute_input/chr{chrom}/genotypes_imp.txt",
            bim   = _acc_out + "/fimpute/fold{fold}/chrom/chr{chrom}.bim",
            idmap = _acc_out + "/fimpute/fold{fold}/fimpute_input/chr{chrom}.id_map.tsv",
        output:
            vcf = temp(_acc_out + "/fimpute/fold{fold}/fimpute/chr{chrom}.vcf.gz"),
            tbi = temp(_acc_out + "/fimpute/fold{fold}/fimpute/chr{chrom}.vcf.gz.tbi"),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy_cv/fimpute_fold{fold}_to_vcf_chr{chrom}.log",
        shell:
            """
            python {workflow.basedir}/scripts/fimpute_io.py to-vcf \
                --imputed {input.imp} \
                --bim {input.bim} \
                --id-map {input.idmap} \
                --out-vcf {output.vcf}.tmp \
                &> {log}
            bgzip -c {output.vcf}.tmp > {output.vcf}
            tabix -f -p vcf {output.vcf}
            """

    rule acc_cv_concat_fimpute:
        input:
            vcfs = lambda wc: expand(
                _acc_out + "/fimpute/fold{fold}/fimpute/chr{chrom}.vcf.gz",
                fold=wc.fold,
                chrom=acc_get_chromosomes(),
            ),
        output:
            vcf = _acc_out + "/fimpute/fold{fold}/fimpute/all_chromosomes.vcf.gz",
            tbi = _acc_out + "/fimpute/fold{fold}/fimpute/all_chromosomes.vcf.gz.tbi",
        conda:
            "../envs/workflow_env.yaml"
        threads: 4
        resources:
            mem_mb = 64000,
            slurm_partition = "r6i-ondemand-4xlarge",
        log:
            "logs/accuracy_cv/fimpute_fold{fold}_concat.log",
        shell:
            """
            (bcftools concat --output {output.vcf} --output-type z \
                --threads {threads} {input.vcfs}) 2> {log}
            tabix -f {output.vcf}
            """

    rule acc_cv_compute_metrics:
        input:
            imputed = _cv_imputed_vcf,
            truth   = _cv_truth_vcf,
        output:
            by_snp   = _acc_out + "/{imputer}/fold{fold}/metrics_by_snp.tsv",
            by_maf   = _acc_out + "/{imputer}/fold{fold}/metrics_by_maf_bin.tsv",
            by_indiv = _acc_out + "/{imputer}/fold{fold}/metrics_by_individual.tsv",
            summary  = _acc_out + "/{imputer}/fold{fold}/summary.tsv",
        params:
            out_dir  = lambda wc: f"{_acc_out}/{wc.imputer}/fold{wc.fold}",
            maf_bins = " ".join(
                str(b) for b in config.get("maf_bins", [0.01, 0.05, 0.1, 0.2, 0.5])
            ),
        conda:
            "../envs/accuracy_env.yaml"
        log:
            "logs/accuracy_cv/{imputer}_fold{fold}_metrics.log",
        shell:
            """
            python {workflow.basedir}/scripts/compute_accuracy_metrics.py \
                --imputed {input.imputed} \
                --truth {input.truth} \
                --out-dir {params.out_dir} \
                --maf-bins {params.maf_bins} \
                &> {log}
            """

    rule acc_cv_aggregate:
        input:
            summaries = expand(
                _acc_out + "/{imputer}/fold{fold}/summary.tsv",
                imputer=_cv_imputers,
                fold=_cv_folds,
            ),
        output:
            summary = _acc_out + "/cv_summary.tsv",
            imputer_summary = _acc_out + "/cv_imputer_summary.tsv",
        params:
            imputers = " ".join(_cv_imputers),
            folds    = " ".join(str(fold) for fold in _cv_folds),
            root     = _acc_out,
        conda:
            "../envs/accuracy_env.yaml"
        log:
            "logs/accuracy_cv/aggregate.log",
        shell:
            """
            python {workflow.basedir}/scripts/aggregate_cv_metrics.py \
                --root {params.root} \
                --imputers {params.imputers} \
                --folds {params.folds} \
                --summary-out {output.summary} \
                --imputer-summary-out {output.imputer_summary} \
                &> {log}
            """

    _impute_bfile = _acc_out + "/unused/to_impute"
    _truth_bfile = _acc_out + "/unused/truth"
    _acc_final_vcf = _acc_out + "/unused/imputed.vcf.gz"

# ---------------------------------------------------------------------------
# Mode: mask_and_impute
# ---------------------------------------------------------------------------

elif _acc_mode == "mask_and_impute":

    rule acc_sample_validation_ids:
        """
        Determine which animals to withhold as the validation set.
        Either copies a user-supplied FID/IID file or randomly samples
        validation_fraction of all animals (reproducible via random_state=42).
        """
        input:
            fam = config["bfile"] + ".fam",
        output:
            ids = _acc_out + "/setup/validation_ids.txt",
        params:
            fraction    = lambda _: config.get("mask", {}).get("validation_fraction", 0.2),
            animal_list = lambda _: config.get("mask", {}).get("validation_animals", ""),
        run:
            import pandas as pd
            import shutil

            fam = pd.read_csv(
                input.fam, sep=r"\s+", header=None,
                names=["fid", "iid", "pid", "mid", "sex", "phen"],
            )
            if params.animal_list:
                shutil.copy(params.animal_list, output.ids)
            else:
                n = max(1, int(len(fam) * params.fraction))
                fam.sample(n, random_state=42)[["fid", "iid"]].to_csv(
                    output.ids, sep="\t", header=False, index=False
                )

    rule acc_sample_snp_list:
        """
        Determine which SNPs form the low-density panel.
        Either copies a user-supplied SNP ID file or randomly samples
        target_n_snps from the bfile .bim (reproducible via random_state=42).
        """
        input:
            bim = config["bfile"] + ".bim",
        output:
            snps = _acc_out + "/setup/ld_snp_list.txt",
        params:
            n_snps   = lambda _: config.get("mask", {}).get("target_n_snps", 50000),
            snp_list = lambda _: config.get("mask", {}).get("target_snp_list", ""),
        run:
            import pandas as pd
            import shutil

            if params.snp_list:
                shutil.copy(params.snp_list, output.snps)
            else:
                bim = pd.read_csv(input.bim, sep="\t", header=None)
                n = min(params.n_snps, len(bim))
                bim[1].sample(n, random_state=42).to_csv(
                    output.snps, header=False, index=False
                )

    rule acc_make_truth_bfile:
        """Validation animals at full density — serves as the truth for comparison."""
        input:
            bed = config["bfile"] + ".bed",
            ids = _acc_out + "/setup/validation_ids.txt",
        output:
            bed = _acc_out + "/truth/hd.bed",
            bim = _acc_out + "/truth/hd.bim",
            fam = _acc_out + "/truth/hd.fam",
        params:
            plink = config["plink_path"],
            bfile = config["bfile"],
            out   = _acc_out + "/truth/hd",
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy/make_truth_bfile.log",
        shell:
            """
            ({params.plink} \
                --bfile {params.bfile} \
                --keep  {input.ids} \
                --make-bed \
                --out   {params.out} \
                --nonfounders --allow-no-sex {params.extra}) &> {log}
            """

    rule acc_make_masked_ld_bfile:
        """
        Validation animals restricted to the low-density SNP panel.
        This is the imputation input; the missing high-density genotypes
        are what the imputer must recover.
        """
        input:
            bed  = config["bfile"] + ".bed",
            ids  = _acc_out + "/setup/validation_ids.txt",
            snps = _acc_out + "/setup/ld_snp_list.txt",
        output:
            bed = _acc_out + "/to_impute/masked.bed",
            bim = _acc_out + "/to_impute/masked.bim",
            fam = _acc_out + "/to_impute/masked.fam",
        params:
            plink = config["plink_path"],
            bfile = config["bfile"],
            out   = _acc_out + "/to_impute/masked",
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy/make_masked_ld_bfile.log",
        shell:
            """
            ({params.plink} \
                --bfile   {params.bfile} \
                --keep    {input.ids} \
                --extract {input.snps} \
                --make-bed \
                --out     {params.out} \
                --nonfounders --allow-no-sex {params.extra}) &> {log}
            """

    if _use_alphaimpute2:
        # AlphaImpute2 imputes by having both HD (reference) and LD (validation)
        # animals in the same dataset.  For mask-and-impute, build that directly
        # from the source BED by setting validation genotypes outside the LD
        # marker panel to missing.  This avoids PLINK2's unfinished
        # non-concatenating --pmerge path.

        rule acc_make_reference_panel:
            """
            Non-validation animals at full high density.
            Merged with the masked LD animals so AlphaImpute2 can use the
            HD individuals as a reference to impute the LD individuals.
            """
            input:
                bed = config["bfile"] + ".bed",
                ids = _acc_out + "/setup/validation_ids.txt",
            output:
                bed = _acc_out + "/reference/panel.bed",
                bim = _acc_out + "/reference/panel.bim",
                fam = _acc_out + "/reference/panel.fam",
            params:
                plink = config["plink_path"],
                bfile = config["bfile"],
                out   = _acc_out + "/reference/panel",
                extra = config.get("plink_extra_flags", ""),
            conda:
                "../envs/workflow_env.yaml"
            log:
                "logs/accuracy/make_reference_panel.log",
            shell:
                """
                ({params.plink} \
                    --bfile  {params.bfile} \
                    --remove {input.ids} \
                    --make-bed \
                    --out    {params.out} \
                    --nonfounders --allow-no-sex {params.extra}) &> {log}
                """

        rule acc_merge_for_alphaimpute2:
            """
            Create a joint AlphaImpute2 input BED.
            Reference animals remain full-density.  Validation animals keep
            LD-panel genotypes and have all other markers set missing, which
            AlphaImpute2 treats as the imputation target.
            """
            input:
                bed  = config["bfile"] + ".bed",
                bim  = config["bfile"] + ".bim",
                fam  = config["bfile"] + ".fam",
                ids  = _acc_out + "/setup/validation_ids.txt",
                snps = _acc_out + "/setup/ld_snp_list.txt",
            output:
                bed = _acc_out + "/to_impute/combined.bed",
                bim = _acc_out + "/to_impute/combined.bim",
                fam = _acc_out + "/to_impute/combined.fam",
            params:
                bfile = config["bfile"],
                out   = _acc_out + "/to_impute/combined",
            conda:
                "../envs/workflow_env.yaml"
            log:
                "logs/accuracy/merge_for_alphaimpute2.log",
            shell:
                """
                python {workflow.basedir}/scripts/mask_validation_genotypes.py \
                    --bfile {params.bfile} \
                    --validation-ids {input.ids} \
                    --ld-snps {input.snps} \
                    --out {params.out} \
                    &> {log}
                """

        _impute_bfile = _acc_out + "/to_impute/combined"
    else:
        _impute_bfile = _acc_out + "/to_impute/masked"

    _truth_bfile = _acc_out + "/truth/hd"

# ---------------------------------------------------------------------------
# Mode: cross_array
# ---------------------------------------------------------------------------

elif _acc_mode == "cross_array":

    rule acc_find_common_animals:
        """
        Intersect .fam files to identify animals present on both arrays.
        Only common animals are used so the comparison is fair.
        """
        input:
            ld_fam = config["cross_array"]["ld_bfile"] + ".fam",
            hd_fam = config["cross_array"]["hd_bfile"] + ".fam",
        output:
            ids = _acc_out + "/setup/common_ids.txt",
        run:
            import pandas as pd

            ld = pd.read_csv(
                input.ld_fam, sep=r"\s+", header=None,
                names=["fid", "iid", "pid", "mid", "sex", "phen"],
            )
            hd = pd.read_csv(
                input.hd_fam, sep=r"\s+", header=None,
                names=["fid", "iid", "pid", "mid", "sex", "phen"],
            )
            common = pd.merge(ld[["fid", "iid"]], hd[["fid", "iid"]], on=["fid", "iid"])
            n_ld, n_hd = len(ld), len(hd)
            print(
                f"LD animals: {n_ld}  HD animals: {n_hd}  "
                f"Common: {len(common)}"
            )
            if len(common) == 0:
                raise ValueError(
                    "No animals in common between LD and HD bfiles. "
                    "Check that IIDs match across the two datasets."
                )
            common.to_csv(output.ids, sep="\t", header=False, index=False)

    rule acc_make_ld_subset:
        """Low-density array, common animals only — imputation input."""
        input:
            bed = config["cross_array"]["ld_bfile"] + ".bed",
            ids = _acc_out + "/setup/common_ids.txt",
        output:
            bed = _acc_out + "/to_impute/ld.bed",
            bim = _acc_out + "/to_impute/ld.bim",
            fam = _acc_out + "/to_impute/ld.fam",
        params:
            plink = config["plink_path"],
            bfile = config["cross_array"]["ld_bfile"],
            out   = _acc_out + "/to_impute/ld",
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy/make_ld_subset.log",
        shell:
            """
            ({params.plink} \
                --bfile {params.bfile} \
                --keep  {input.ids} \
                --make-bed \
                --out   {params.out} \
                --nonfounders --allow-no-sex {params.extra}) &> {log}
            """

    rule acc_make_hd_truth:
        """High-density array, common animals only — serves as the truth."""
        input:
            bed = config["cross_array"]["hd_bfile"] + ".bed",
            ids = _acc_out + "/setup/common_ids.txt",
        output:
            bed = _acc_out + "/truth/hd.bed",
            bim = _acc_out + "/truth/hd.bim",
            fam = _acc_out + "/truth/hd.fam",
        params:
            plink = config["plink_path"],
            bfile = config["cross_array"]["hd_bfile"],
            out   = _acc_out + "/truth/hd",
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy/make_hd_truth.log",
        shell:
            """
            ({params.plink} \
                --bfile {params.bfile} \
                --keep  {input.ids} \
                --make-bed \
                --out   {params.out} \
                --nonfounders --allow-no-sex {params.extra}) &> {log}
            """

    _impute_bfile = _acc_out + "/to_impute/ld"
    _truth_bfile  = _acc_out + "/truth/hd"

else:
    raise ValueError(
        f"Unknown accuracy_mode: {_acc_mode!r}. "
        "Use 'mask_and_impute', 'cross_array', or 'kfold_mask_and_impute'."
    )

# ---------------------------------------------------------------------------
# Shared: imputation — Beagle path (per-chromosome)
# ---------------------------------------------------------------------------

if not _use_alphaimpute2:

    rule acc_make_per_chrom_vcf:
        input:
            bed = _impute_bfile + ".bed",
        output:
            vcf = temp(_acc_out + "/dedup/chr{chrom}.vcf.gz"),
        params:
            plink = config["plink_path"],
            bfile = _impute_bfile,
            out   = _acc_out + "/dedup/chr{chrom}",
            chr   = "{chrom}",
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy/dedup_chr{chrom}.log",
        shell:
            """
            ({params.plink} \
                --nonfounders --allow-no-sex \
                --bfile {params.bfile} \
                --chr   {params.chr} \
                --export vcf bgz {params.extra} \
                --out    {params.out} \
                --snps-only) &> {log}
            """

    rule acc_normalize_vcf:
        input:
            vcf = _acc_out + "/dedup/chr{chrom}.vcf.gz",
        output:
            vcf = temp(_acc_out + "/normalized/chr{chrom}.vcf.gz"),
            tbi = temp(_acc_out + "/normalized/chr{chrom}.vcf.gz.tbi"),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy/normalize_chr{chrom}.log",
        shell:
            """
            (bcftools norm -d snps {input.vcf} | bgzip > {output.vcf}) 2> {log}
            tabix -f -p vcf {output.vcf}
            """

    rule acc_run_beagle:
        input:
            vcf = _acc_out + "/normalized/chr{chrom}.vcf.gz",
        output:
            vcf = temp(_acc_out + "/imputed/chr{chrom}.vcf.gz"),
            tbi = temp(_acc_out + "/imputed/chr{chrom}.vcf.gz.tbi"),
        params:
            beagle    = config["beagle_jar"],
            window    = config["beagle_params"]["window"],
            overlap   = config["beagle_params"]["overlap"],
            ne        = config["beagle_params"]["ne"],
            outbase   = lambda wc: _acc_out + f"/imputed/chr{wc.chrom}",
            # Use reference_vcf from the main config if provided.
            # For mask_and_impute without an external reference, Beagle uses
            # population LD within the study sample alone.
            ref_param = (
                ""
                if _acc_mode == "mask_and_impute"
                else (
                    f"ref={config['reference_vcf']}"
                    if config.get("reference_vcf", "").strip()
                    else ""
                )
            ),
        threads:
            config["beagle_params"]["nthreads"]
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy/beagle_chr{chrom}.log",
        resources:
            mem_mb = 70000,
        shell:
            """
            (java -Xmx{resources.mem_mb}m -jar {params.beagle} \
                gt={input.vcf} {params.ref_param} \
                window={params.window} \
                overlap={params.overlap} \
                out={params.outbase} \
                nthreads={threads} \
                ne={params.ne} \
                chrom={wildcards.chrom}) &> {log}
            tabix -f {output.vcf}
            """

    rule acc_concat_imputed:
        input:
            vcfs = lambda _: expand(
                _acc_out + "/imputed/chr{chrom}.vcf.gz",
                chrom=acc_get_chromosomes(),
            ),
        output:
            vcf = _acc_out + "/imputed/all_chromosomes.vcf.gz",
            tbi = _acc_out + "/imputed/all_chromosomes.vcf.gz.tbi",
        conda:
            "../envs/workflow_env.yaml"
        threads: 4
        resources:
            mem_mb = 64000,
        log:
            "logs/accuracy/concat_imputed.log",
        shell:
            """
            (bcftools concat \
                --output      {output.vcf} \
                --output-type z \
                --threads     {threads} \
                {input.vcfs}) 2> {log}
            tabix -f {output.vcf}
            """

    _acc_final_vcf = _acc_out + "/imputed/all_chromosomes.vcf.gz"

# ---------------------------------------------------------------------------
# Shared: imputation — AlphaImpute2 path (genome-wide)
# ---------------------------------------------------------------------------

else:

    rule acc_make_alphaimpute2_chrom_bfile:
        """Split the joint AlphaImpute2 PLINK input by chromosome."""
        input:
            bed = _impute_bfile + ".bed",
            bim = _impute_bfile + ".bim",
            fam = _impute_bfile + ".fam",
        output:
            bed = temp(_acc_out + "/alphaimpute2_chrom/chr{chrom}.bed"),
            bim = temp(_acc_out + "/alphaimpute2_chrom/chr{chrom}.bim"),
            fam = temp(_acc_out + "/alphaimpute2_chrom/chr{chrom}.fam"),
        params:
            plink = config["plink_path"],
            bfile = _impute_bfile,
            out   = _acc_out + "/alphaimpute2_chrom/chr{chrom}",
            chr   = "{chrom}",
            extra = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy/alphaimpute2_chrom_chr{chrom}.log",
        shell:
            """
            ({params.plink} \
                --bfile {params.bfile} \
                --chr {params.chr} \
                --make-bed \
                --out {params.out} \
                --nonfounders --allow-no-sex {params.extra}) &> {log}
            """

    rule acc_plink_to_alphaimpute2_fmt:
        """
        Export one chromosome of the imputation-input PLINK file to AlphaImpute2
        text format.  AlphaImpute2 is run chromosome-wise to avoid treating all
        chromosomes as one long linkage group and to avoid genome-wide segfaults.
        """
        input:
            bed = _acc_out + "/alphaimpute2_chrom/chr{chrom}.bed",
            bim = _acc_out + "/alphaimpute2_chrom/chr{chrom}.bim",
            fam = _acc_out + "/alphaimpute2_chrom/chr{chrom}.fam",
        output:
            genotypes = _acc_out + "/alphaimpute2_input/chr{chrom}.genotypes.txt",
            pedigree  = _acc_out + "/alphaimpute2_input/chr{chrom}.pedigree.txt",
        params:
            plink      = config["plink_path"],
            bfile      = _acc_out + "/alphaimpute2_chrom/chr{chrom}",
            raw_prefix = _acc_out + "/alphaimpute2_input/chr{chrom}.raw_tmp",
            extra      = config.get("plink_extra_flags", ""),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy/plink_to_alphaimpute2_chr{chrom}.log",
        shell:
            """
            (
            mkdir -p "$(dirname {output.genotypes})"

            {params.plink} --bfile {params.bfile} --export A \
                {params.extra} --allow-no-sex --out {params.raw_prefix}

            awk 'NR>1 {{
                printf $2;
                for (i=7; i<=NF; i++) {{
                    if ($i == "NA") printf " 9";
                    else printf " %d", $i;
                }}
                printf "\\n";
            }}' {params.raw_prefix}.raw > {output.genotypes}

            awk '{{print $2, $3, $4}}' {input.fam} > {output.pedigree}
            ) &> {log}
            """

    rule acc_run_alphaimpute2:
        input:
            genotypes = _acc_out + "/alphaimpute2_input/chr{chrom}.genotypes.txt",
            pedigree  = _acc_out + "/alphaimpute2_input/chr{chrom}.pedigree.txt",
        output:
            genotypes = temp(_acc_out + "/alphaimpute2_output/chr{chrom}.genotypes"),
        params:
            out_prefix  = _acc_out + "/alphaimpute2_output/chr{chrom}",
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
            "logs/accuracy/run_alphaimpute2_chr{chrom}.log",
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
                {params.extra_flags}) &> {log}
            """

    rule acc_alphaimpute2_to_vcf:
        """
        Convert per-chromosome AlphaImpute2 output to VCF using the matching
        chromosome .bim for variant metadata.
        """
        input:
            genotypes = _acc_out + "/alphaimpute2_output/chr{chrom}.genotypes",
            bim       = _acc_out + "/alphaimpute2_chrom/chr{chrom}.bim",
        output:
            vcf = temp(_acc_out + "/alphaimpute2/chr{chrom}.vcf.gz"),
            tbi = temp(_acc_out + "/alphaimpute2/chr{chrom}.vcf.gz.tbi"),
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy/alphaimpute2_to_vcf_chr{chrom}.log",
        resources:
            mem_mb = 32000,
        script:
            "../scripts/alphaimpute2_to_vcf.py"

    rule acc_concat_alphaimpute2:
        input:
            vcfs = lambda _: expand(
                _acc_out + "/alphaimpute2/chr{chrom}.vcf.gz",
                chrom=acc_get_chromosomes(),
            ),
        output:
            vcf = _acc_out + "/alphaimpute2/all_chromosomes.vcf.gz",
            tbi = _acc_out + "/alphaimpute2/all_chromosomes.vcf.gz.tbi",
        conda:
            "../envs/workflow_env.yaml"
        threads: 4
        resources:
            mem_mb = 64000,
        log:
            "logs/accuracy/concat_alphaimpute2.log",
        shell:
            """
            (bcftools concat \
                --output      {output.vcf} \
                --output-type z \
                --threads     {threads} \
                {input.vcfs}) 2> {log}
            tabix -f {output.vcf}
            """

    _acc_final_vcf = _acc_out + "/alphaimpute2/all_chromosomes.vcf.gz"

# ---------------------------------------------------------------------------
# Shared: convert truth PLINK to genome-wide VCF
# ---------------------------------------------------------------------------

rule acc_truth_to_vcf:
    """Export truth PLINK file to a genome-wide bgzipped VCF for comparison."""
    input:
        bed = _truth_bfile + ".bed",
    output:
        vcf = _acc_out + "/truth/all_chromosomes.vcf.gz",
        tbi = _acc_out + "/truth/all_chromosomes.vcf.gz.tbi",
    params:
        plink = config["plink_path"],
        bfile = _truth_bfile,
        out   = _acc_out + "/truth/all_chromosomes",
        extra = config.get("plink_extra_flags", ""),
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/accuracy/truth_to_vcf.log",
    shell:
        """
        ({params.plink} \
            --bfile {params.bfile} \
            --export vcf bgz \
            --nonfounders --allow-no-sex \
            {params.extra} \
            --out    {params.out} \
            --snps-only) &> {log}
        tabix -f -p vcf {output.vcf}
        """

# ---------------------------------------------------------------------------
# Shared: compute accuracy metrics
# ---------------------------------------------------------------------------

rule acc_compute_metrics:
    """
    Compare imputed genotypes against truth genotypes and write accuracy tables.

    Outputs:
      metrics_by_snp.tsv         allelic_r2, concordance, MAF, n_evaluated per SNP
      metrics_by_maf_bin.tsv     mean/median R² and concordance stratified by MAF
      metrics_by_individual.tsv  mean concordance per animal
      summary.tsv                overall mean R², % SNPs with R² ≥ 0.8/0.9, concordance
    """
    input:
        imputed = _acc_final_vcf,
        truth   = _acc_out + "/truth/all_chromosomes.vcf.gz",
    output:
        by_snp   = _acc_out + "/metrics_by_snp.tsv",
        by_maf   = _acc_out + "/metrics_by_maf_bin.tsv",
        by_indiv = _acc_out + "/metrics_by_individual.tsv",
        summary  = _acc_out + "/summary.tsv",
    params:
        out_dir  = _acc_out,
        maf_bins = " ".join(
            str(b) for b in config.get("maf_bins", [0.01, 0.05, 0.1, 0.2, 0.5])
        ),
    conda:
        "../envs/accuracy_env.yaml"
    log:
        "logs/accuracy/compute_metrics.log",
    shell:
        """
        python {workflow.basedir}/scripts/compute_accuracy_metrics.py \
            --imputed  {input.imputed} \
            --truth    {input.truth} \
            --out-dir  {params.out_dir} \
            --maf-bins {params.maf_bins} \
            &> {log}
        """
