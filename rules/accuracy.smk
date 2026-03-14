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

# ---------------------------------------------------------------------------
# Mode: mask_and_impute
# ---------------------------------------------------------------------------

if _acc_mode == "mask_and_impute":

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
        # animals in the same dataset.  Create the reference panel and merge.

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
            Merge HD reference panel + masked LD validation animals.
            After the merge, validation animals have missing values at all
            non-LD positions, which AlphaImpute2 treats as the imputation
            target (non-HD animals according to hd_threshold).
            """
            input:
                ref_bed = _acc_out + "/reference/panel.bed",
                ref_bim = _acc_out + "/reference/panel.bim",
                ref_fam = _acc_out + "/reference/panel.fam",
                ld_bed  = _acc_out + "/to_impute/masked.bed",
                ld_bim  = _acc_out + "/to_impute/masked.bim",
                ld_fam  = _acc_out + "/to_impute/masked.fam",
            output:
                bed = _acc_out + "/to_impute/combined.bed",
                bim = _acc_out + "/to_impute/combined.bim",
                fam = _acc_out + "/to_impute/combined.fam",
            params:
                plink    = config["plink_path"],
                ref_pref = _acc_out + "/reference/panel",
                ld_pref  = _acc_out + "/to_impute/masked",
                out      = _acc_out + "/to_impute/combined",
                extra    = config.get("plink_extra_flags", ""),
            log:
                "logs/accuracy/merge_for_alphaimpute2.log",
            shell:
                """
                ({params.plink} \
                    --bfile   {params.ref_pref} \
                    --bmerge  {params.ld_pref} \
                    --make-bed \
                    --out     {params.out} \
                    --nonfounders --allow-no-sex {params.extra}) &> {log}
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
        "Use 'mask_and_impute' or 'cross_array'."
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
                f"ref={config['reference_vcf']}"
                if config.get("reference_vcf", "").strip()
                else ""
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

    rule acc_plink_to_alphaimpute2_fmt:
        """
        Export the imputation-input PLINK file to AlphaImpute2 text format.
        For mask_and_impute this is the merged HD-reference + masked-LD file,
        so AlphaImpute2 sees both HD and LD animals together and can use the
        HD animals as a reference panel during pedigree/population imputation.
        """
        input:
            bed = _impute_bfile + ".bed",
            bim = _impute_bfile + ".bim",
            fam = _impute_bfile + ".fam",
        output:
            genotypes = _acc_out + "/alphaimpute2_input/genotypes.txt",
            pedigree  = _acc_out + "/alphaimpute2_input/pedigree.txt",
        params:
            plink      = config["plink_path"],
            bfile      = _impute_bfile,
            raw_prefix = _acc_out + "/alphaimpute2_input/raw_tmp",
            extra      = config.get("plink_extra_flags", ""),
        log:
            "logs/accuracy/plink_to_alphaimpute2.log",
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
            genotypes = _acc_out + "/alphaimpute2_input/genotypes.txt",
            pedigree  = _acc_out + "/alphaimpute2_input/pedigree.txt",
        output:
            genotypes = _acc_out + "/alphaimpute2_output/imputed.genotypes",
        params:
            out_prefix  = _acc_out + "/alphaimpute2_output/imputed",
            cycles      = config["alphaimpute2_params"].get("cycles", 4),
            threshold   = config["alphaimpute2_params"].get("final_peeling_threshold", 0.1),
            hd_thresh   = config["alphaimpute2_params"].get("hd_threshold", 0.95),
            length      = config["alphaimpute2_params"].get("length", 1.0),
            extra_flags = " ".join(
                (["-ped_only"] if config["alphaimpute2_params"].get("ped_only", False) else [])
                + (["-pop_only"] if config["alphaimpute2_params"].get("pop_only", False) else [])
            ),
        threads:
            config["alphaimpute2_params"].get("maxthreads", 1)
        conda:
            "../envs/alphaimpute2_env.yaml"
        log:
            "logs/accuracy/run_alphaimpute2.log",
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
        Convert AlphaImpute2 output to VCF using the imputation-input .bim
        for variant metadata (contains all SNPs after the merge step).
        """
        input:
            genotypes = _acc_out + "/alphaimpute2_output/imputed.genotypes",
            bim       = _impute_bfile + ".bim",
        output:
            vcf = _acc_out + "/alphaimpute2/all_chromosomes.vcf.gz",
            tbi = _acc_out + "/alphaimpute2/all_chromosomes.vcf.gz.tbi",
        conda:
            "../envs/workflow_env.yaml"
        log:
            "logs/accuracy/alphaimpute2_to_vcf.log",
        resources:
            mem_mb = 32000,
        script:
            "../scripts/alphaimpute2_to_vcf.py"

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
        mv {params.out}.vcf.gz {output.vcf}
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
        python scripts/compute_accuracy_metrics.py \
            --imputed  {input.imputed} \
            --truth    {input.truth} \
            --out-dir  {params.out_dir} \
            --maf-bins {params.maf_bins} \
            &> {log}
        """
