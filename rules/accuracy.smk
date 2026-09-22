# rules/accuracy.smk — imputation accuracy evaluation rules
#
# Included exclusively by Snakefile_accuracy; never by the main Snakefile.
#
# Both modes are animal-level K-fold CV and share one rule family (acc_cv_*).
# They differ only in where each fold's inputs are cut from:
#
#   masked LD input   the array the held-out animals were typed on
#   reference panel   the high-density array, minus the held-out fold
#   truth             the high-density array, held-out fold, LD panel removed
#
# In kfold_mask_and_impute both sources are the one bfile, so the LD input is
# that same fileset thinned to the panel. In cross_array they are two real
# arrays, which is the only construction in which a V3 -> V4 test can run at
# all: practically every V4 fish is also V3-typed, so no external V4 panel
# exists to hold out.
#
# Layout under accuracy_output_dir:
#
#   setup/                  fold assignment, LD panel, identity gate (cross_array)
#   folds/fold{N}/          everything a fold needs that does not depend on the
#                           imputer: truth, reference panel, masked inputs,
#                           per-chromosome splits. Built once, shared by every
#                           engine in cv.imputers.
#   {imputer}/fold{N}/      that engine's imputed VCF and its metrics
#   logs/

_acc_out = config.get("accuracy_output_dir", "accuracy_output")
_log     = _acc_out + "/logs/"


def _as_list(value, default):
    if value is None:
        value = default
    if isinstance(value, (int, float)):
        return [str(value)]
    if isinstance(value, str):
        return [item for item in value.replace(",", " ").split() if item]
    return [str(v) for v in value]


_known_imputers = {"beagle", "alphaimpute2", "fimpute"}
_cv_imputers = _as_list(nested_config("cv", "imputers"), ["beagle", "fimpute"])
_unsupported = sorted(set(_cv_imputers) - _known_imputers)
if _unsupported or not _cv_imputers:
    raise ValueError(
        f"Unsupported cv imputers: {', '.join(_unsupported) or '(none given)'}. "
        "Use beagle, alphaimpute2, and/or fimpute."
    )

_cv_n_folds = int(nested_config("cv", "n_folds", 10))
# A subset of folds to actually run. [1] is a single hold-out of 1/n_folds of
# the animals -- a quick check before committing to the full K-fold.
_cv_folds = [int(f) for f in _as_list(nested_config("cv", "folds_to_run"), [])] \
    or list(range(1, _cv_n_folds + 1))
if not set(_cv_folds) <= set(range(1, _cv_n_folds + 1)):
    raise ValueError(f"cv.folds_to_run {_cv_folds} must lie within 1..{_cv_n_folds}.")

_cv_target_n_snps = int(nested_config("cv", "target_n_snps", 10000))
_cv_seed = int(nested_config("cv", "random_seed", 42))
_cv_target_snp_list = str(nested_config("cv", "target_snp_list", "") or "")
_cv_reference_max_animals = int(nested_config("cv", "reference_max_animals", 0))
_cv_reliable_r2 = float(nested_config("cv", "reliable_r2_threshold", 0.90))

_fimpute_executable = nested_config(
    "fimpute_params", "executable", "/mnt/efshome/applications/FImpute3/2026/FImpute3"
)
_fimpute_nthreads = int(nested_config("fimpute_params", "nthreads", 1))

_cv_identity_threshold = float(nested_config("cross_array", "identity_threshold", 0.90))
_cv_min_pass_rate = float(nested_config("cross_array", "min_identity_pass_rate", 0.5))
_cv_min_shared_markers = int(nested_config("cross_array", "min_shared_markers", 100))

_cv_cross = _acc_mode == "cross_array"
if _cv_cross:
    _cv_ld_bfile = _cross_array_bfile("ld")
    _cv_hd_bfile = _cross_array_bfile("hd")
else:
    _cv_ld_bfile = config["bfile"]
    _cv_hd_bfile = config["bfile"]

# The LD panel is fixed by the two arrays in cross_array -- it is exactly the
# markers both carry -- so acc_cv_setup is handed that list instead of sampling
# one. Everything downstream reads setup/ld_snp_list.txt either way.
_cv_setup_snp_list = (
    _acc_out + "/setup/shared_snps.txt" if _cv_cross else _cv_target_snp_list
)
_cv_setup_fam = (
    _acc_out + "/setup/animals.txt" if _cv_cross else config["bfile"] + ".fam"
)

# The HD fileset restricted to the animals the identity gate passed. Beagle's
# panel and truth are cut with --keep off folds.tsv and are therefore already
# gated, but the combined fileset FImpute and AlphaImpute2 eat is built by
# masking the whole HD file -- so without this it hands them the failed animals
# as reference donors, complete with the pedigree links the gate just declared
# unreliable, and the engines are no longer scored on the same cohort.
_cv_cohort_bfile = _acc_out + "/setup/cohort" if _cv_cross else config["bfile"]

_plink = config["plink_path"]
_plink_extra = config.get("plink_extra_flags", "")


def _fold(stem):
    """Path of an imputer-independent per-fold file."""
    return _acc_out + "/folds/fold{fold}/" + stem


# Fold numbers are digits and imputers are named engines. Without this,
# setup/fold{fold}_ids.txt also matches setup/fold1_reference_ids.txt with
# fold="1_reference", and the two rules collide.
wildcard_constraints:
    fold = r"\d+",
    imputer = "|".join(_cv_imputers),
    chrom = r"[^/]+",


# ---------------------------------------------------------------------------
# cross_array only: pair the two arrays and gate on identity
# ---------------------------------------------------------------------------

if _cv_cross:

    rule acc_cv_pair_animals:
        """Animals carried by both arrays, matched on individual ID."""
        input:
            ld_fam = _cv_ld_bfile + ".fam",
            hd_fam = _cv_hd_bfile + ".fam",
        output:
            pairs = _acc_out + "/setup/paired_ids.txt",
        conda:
            "../envs/workflow_env.yaml"
        log:
            _log + "pair_animals.log",
        shell:
            """
            python {workflow.basedir}/scripts/pair_animals.py \
                --ld-fam {input.ld_fam} \
                --hd-fam {input.hd_fam} \
                --out {output.pairs} \
                &> {log}
            """

    rule acc_cv_shared_markers:
        """
        The markers both arrays carry -- the LD panel, and nothing else.

        Forcing the panel to be LD n HD does two jobs at once. An LD-only
        probe would be handed to the imputer as an observed genotype and
        then never scored, because the truth fileset does not contain it;
        and an HD marker absent from LD is exactly what the run is meant to
        recover.
        """
        input:
            ld_bim = _cv_ld_bfile + ".bim",
            hd_bim = _cv_hd_bfile + ".bim",
        output:
            snps = _acc_out + "/setup/shared_snps.txt",
        params:
            min_shared = _cv_min_shared_markers,
        conda:
            "../envs/workflow_env.yaml"
        log:
            _log + "shared_markers.log",
        shell:
            """
            python {workflow.basedir}/scripts/shared_marker_list.py \
                --bim {input.ld_bim} {input.hd_bim} \
                --out {output.snps} \
                --min-shared {params.min_shared} \
                &> {log}
            """

    rule acc_cv_identity_export:
        """Both arrays cut to the shared markers, ready for a direct comparison."""
        input:
            bed  = lambda wc: (_cv_ld_bfile if wc.side == "ld" else _cv_hd_bfile) + ".bed",
            ids  = _acc_out + "/setup/paired_ids.txt",
            snps = _acc_out + "/setup/shared_snps.txt",
        output:
            vcf = temp(_acc_out + "/setup/identity/{side}.vcf.gz"),
            tbi = temp(_acc_out + "/setup/identity/{side}.vcf.gz.tbi"),
        wildcard_constraints:
            side = "ld|hd",
        params:
            bfile = lambda wc: _cv_ld_bfile if wc.side == "ld" else _cv_hd_bfile,
            out   = lambda wc: f"{_acc_out}/setup/identity/{wc.side}_raw",
        conda:
            "../envs/workflow_env.yaml"
        log:
            _log + "identity_export_{side}.log",
        shell:
            """
            ({_plink} --bfile {params.bfile} --keep {input.ids} \
                --extract {input.snps} --snps-only \
                --export vcf bgz id-paste=iid --out {params.out} \
                {_plink_extra}) &> {log}
            bcftools norm -d snps {params.out}.vcf.gz 2>> {log} | bgzip > {output.vcf}
            tabix -f -p vcf {output.vcf}
            rm -f {params.out}.vcf.gz
            """

    rule acc_cv_identity_metrics:
        """LD-vs-HD concordance per animal, at the markers both arrays measured."""
        input:
            ld = _acc_out + "/setup/identity/ld.vcf.gz",
            hd = _acc_out + "/setup/identity/hd.vcf.gz",
            # compute_accuracy_metrics.py reads through "bcftools view -R",
            # which needs the index; temp() indexes nothing names as input are
            # deleted as soon as the producing job finishes.
            ld_tbi = _acc_out + "/setup/identity/ld.vcf.gz.tbi",
            hd_tbi = _acc_out + "/setup/identity/hd.vcf.gz.tbi",
        output:
            by_indiv = _acc_out + "/setup/identity/metrics_by_individual.tsv",
        params:
            out_dir = _acc_out + "/setup/identity",
        conda:
            "../envs/accuracy_env.yaml"
        # The whole paired cohort at every shared marker: the largest single
        # matrix in the run, bigger than any per-fold comparison.
        resources:
            mem_mb = 32000,
            slurm_partition = "r7i-ondemand-2xlarge",
        log:
            _log + "identity_metrics.log",
        shell:
            """
            python {workflow.basedir}/scripts/compute_accuracy_metrics.py \
                --imputed {input.ld} \
                --truth {input.hd} \
                --out-dir {params.out_dir} \
                &> {log}
            """

    rule acc_cv_identity_gate:
        """
        Drop animals whose two typings are not the same physical fish.

        Mandatory, not a report. On the 2026-08 salmon benchmark 974 of
        3,881 animals in one step had LD and HD genotypes from different
        fish; the step scored mean R2 0.53 and looked like a result.
        """
        input:
            pairs   = _acc_out + "/setup/paired_ids.txt",
            metrics = _acc_out + "/setup/identity/metrics_by_individual.tsv",
        output:
            identity = _acc_out + "/setup/identity.tsv",
            failed   = _acc_out + "/setup/identity_fail.ids",
            animals  = _acc_out + "/setup/animals.txt",
        params:
            threshold     = _cv_identity_threshold,
            min_pass_rate = _cv_min_pass_rate,
        conda:
            "../envs/accuracy_env.yaml"
        log:
            _log + "identity_gate.log",
        shell:
            """
            python {workflow.basedir}/scripts/identity_gate.py \
                --pairs {input.pairs} \
                --metrics {input.metrics} \
                --threshold {params.threshold} \
                --min-pass-rate {params.min_pass_rate} \
                --identity-out {output.identity} \
                --fail-out {output.failed} \
                --keep-out {output.animals} \
                &> {log}
            """

    rule acc_cv_cohort_bfile:
        """The high-density fileset cut down to the animals that passed the gate."""
        input:
            bed     = _cv_hd_bfile + ".bed",
            animals = _acc_out + "/setup/animals.txt",
        output:
            bed = _acc_out + "/setup/cohort.bed",
            bim = _acc_out + "/setup/cohort.bim",
            fam = _acc_out + "/setup/cohort.fam",
        params:
            out = _acc_out + "/setup/cohort",
        conda:
            "../envs/workflow_env.yaml"
        log:
            _log + "cohort_bfile.log",
        shell:
            """
            ({_plink} --bfile {_cv_hd_bfile} --keep {input.animals} \
                --make-bed --out {params.out} {_plink_extra}) &> {log}
            """


# ---------------------------------------------------------------------------
# Fold assignment and the shared LD panel
# ---------------------------------------------------------------------------

rule acc_cv_setup:
    """Create the shared fold assignment and LD marker panel."""
    input:
        fam = _cv_setup_fam,
        bim = _cv_hd_bfile + ".bim",
        snp_list = [_acc_out + "/setup/shared_snps.txt"] if _cv_cross else [],
    output:
        folds = _acc_out + "/setup/folds.tsv",
        snps  = _acc_out + "/setup/ld_snp_list.txt",
    params:
        n_folds  = _cv_n_folds,
        n_snps   = _cv_target_n_snps,
        seed     = _cv_seed,
        snp_list = _cv_setup_snp_list,
    conda:
        "../envs/workflow_env.yaml"
    log:
        _log + "setup.log",
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
        folds.loc[folds["fold"] == int(wildcards.fold), ["fid", "iid"]].to_csv(
            output.ids, sep="\t", header=False, index=False
        )


rule acc_cv_reference_ids:
    """
    Trim the reference panel to a fixed size, when a size is asked for.

    Panel size is the axis worth sweeping: accuracy climbs with it and then
    plateaus, and past the plateau the panel is only costing compute. Every
    fold draws from the same seed, so a size is comparable across folds and
    across truth pairs. cv.reference_max_animals = 0 keeps everyone.
    """
    input:
        folds = _acc_out + "/setup/folds.tsv",
    output:
        ids = _acc_out + "/setup/fold{fold}_reference_ids.txt",
    params:
        max_animals = _cv_reference_max_animals,
        seed = _cv_seed,
    run:
        import pandas as pd

        folds = pd.read_csv(input.folds, sep="\t")
        panel = folds.loc[folds["fold"] != int(wildcards.fold), ["fid", "iid"]]
        if params.max_animals and len(panel) > params.max_animals:
            panel = panel.sample(
                n=params.max_animals, random_state=params.seed
            ).sort_values(["fid", "iid"])
        panel.to_csv(output.ids, sep="\t", header=False, index=False)


# ---------------------------------------------------------------------------
# Per-fold inputs, shared by every imputer
# ---------------------------------------------------------------------------

rule acc_cv_truth_vcf:
    """
    Held-out animals at high density, LD panel removed -- the truth.

    The LD panel markers are handed to the imputer as observed genotypes, so
    they come back unchanged and would score r2 ~ 1.0 by construction.
    Excluding them keeps the metrics on markers the imputer had to recover.
    """
    input:
        bed  = _cv_hd_bfile + ".bed",
        ids  = _acc_out + "/setup/fold{fold}_ids.txt",
        snps = _acc_out + "/setup/ld_snp_list.txt",
    output:
        vcf = _fold("truth.vcf.gz"),
        tbi = _fold("truth.vcf.gz.tbi"),
    params:
        out = _fold("truth"),
    conda:
        "../envs/workflow_env.yaml"
    log:
        _log + "fold{fold}_truth.log",
    shell:
        """
        ({_plink} --bfile {_cv_hd_bfile} --keep {input.ids} \
            --exclude {input.snps} --snps-only \
            --export vcf bgz id-paste=iid --out {params.out} \
            {_plink_extra}) &> {log}
        tabix -f -p vcf {output.vcf}
        """


rule acc_cv_masked_ld_bfile:
    """
    The held-out animals as Beagle will see them: LD density only.

    In cross_array these are the animals' real low-density typings, and the
    allele coding is forced onto the HD .bim. Two exports of the same marker
    can disagree about which allele is A1, and nothing downstream notices: an
    allele flip leaves allelic r2 looking reasonable, because squaring the
    correlation hides the sign, and surfaces only as collapsed concordance.
    plink2 refuses outright if an HD allele is absent from the LD marker,
    which is the loud failure that case deserves.
    """
    input:
        bed  = _cv_ld_bfile + ".bed",
        ids  = _acc_out + "/setup/fold{fold}_ids.txt",
        snps = _acc_out + "/setup/ld_snp_list.txt",
        hd_bim = [_cv_hd_bfile + ".bim"] if _cv_cross else [],
    output:
        bed = temp(_fold("masked.bed")),
        bim = temp(_fold("masked.bim")),
        fam = temp(_fold("masked.fam")),
    params:
        out = _fold("masked"),
        # Column 6 of a .bim is A2, column 2 the marker ID.
        orient = f"--ref-allele force {_cv_hd_bfile}.bim 6 2" if _cv_cross else "",
    conda:
        "../envs/workflow_env.yaml"
    log:
        _log + "fold{fold}_masked_ld.log",
    shell:
        """
        ({_plink} --bfile {_cv_ld_bfile} --keep {input.ids} \
            --extract {input.snps} {params.orient} \
            --make-bed --out {params.out} {_plink_extra}) &> {log}
        """


rule acc_cv_reference_bfile:
    """The other folds at full density -- the panel Beagle draws haplotypes from."""
    input:
        bed = _cv_hd_bfile + ".bed",
        ids = _acc_out + "/setup/fold{fold}_reference_ids.txt",
    output:
        bed = temp(_fold("reference.bed")),
        bim = temp(_fold("reference.bim")),
        fam = temp(_fold("reference.fam")),
    params:
        out = _fold("reference"),
    conda:
        "../envs/workflow_env.yaml"
    log:
        _log + "fold{fold}_reference.log",
    shell:
        """
        ({_plink} --bfile {_cv_hd_bfile} --keep {input.ids} \
            --make-bed --out {params.out} {_plink_extra}) &> {log}
        """


rule acc_cv_combined_bfile:
    """
    One fileset: the held-out fold at LD density, everyone else at HD.

    FImpute and AlphaImpute2 take a single input rather than a separate panel,
    so the panel has to be folded in. plink2's --pmerge is unfinished, hence
    the byte-level rewrite.

    In cross_array the held-out animals' panel genotypes are spliced in from
    the array they were really typed on. Leaving the HD array's own calls
    there would delete the array-transition effect the run exists to measure
    and quietly turn it into a masking test.
    """
    input:
        bed  = _cv_cohort_bfile + ".bed",
        bim  = _cv_cohort_bfile + ".bim",
        fam  = _cv_cohort_bfile + ".fam",
        ids  = _acc_out + "/setup/fold{fold}_ids.txt",
        snps = _acc_out + "/setup/ld_snp_list.txt",
        ld   = [_cv_ld_bfile + ".bed"] if _cv_cross else [],
    output:
        bed = temp(_fold("combined.bed")),
        bim = temp(_fold("combined.bim")),
        fam = temp(_fold("combined.fam")),
    params:
        out     = _fold("combined"),
        replace = f"--replace-from {_cv_ld_bfile}" if _cv_cross else "",
    conda:
        "../envs/workflow_env.yaml"
    log:
        _log + "fold{fold}_combined.log",
    shell:
        """
        python {workflow.basedir}/scripts/mask_validation_genotypes.py \
            --bfile {_cv_cohort_bfile} \
            --validation-ids {input.ids} \
            --ld-snps {input.snps} \
            {params.replace} \
            --out {params.out} \
            &> {log}
        """


rule acc_cv_combined_chrom:
    """One chromosome of the combined fileset, plus its allele-count export."""
    input:
        bed = _fold("combined.bed"),
        bim = _fold("combined.bim"),
        fam = _fold("combined.fam"),
    output:
        bed = temp(_fold("chrom/chr{chrom}.bed")),
        bim = temp(_fold("chrom/chr{chrom}.bim")),
        fam = temp(_fold("chrom/chr{chrom}.fam")),
        raw = temp(_fold("chrom/chr{chrom}.raw")),
    params:
        src = _fold("combined"),
        out = _fold("chrom/chr{chrom}"),
    conda:
        "../envs/workflow_env.yaml"
    log:
        _log + "fold{fold}_chrom_chr{chrom}.log",
    shell:
        """
        ({_plink} --bfile {params.src} --chr {wildcards.chrom} \
            --make-bed --out {params.out} {_plink_extra}
         {_plink} --bfile {params.out} --export A \
            --out {params.out} {_plink_extra}) &> {log}
        """


# ---------------------------------------------------------------------------
# Beagle: held-out LD animals as gt=, the phased remaining folds as ref=
# ---------------------------------------------------------------------------

rule acc_cv_beagle_vcf:
    """One chromosome of the masked targets (gt) or the reference panel (ref)."""
    input:
        bed = lambda wc: _fold("masked.bed" if wc.panel == "gt" else "reference.bed"),
        bim = lambda wc: _fold("masked.bim" if wc.panel == "gt" else "reference.bim"),
        fam = lambda wc: _fold("masked.fam" if wc.panel == "gt" else "reference.fam"),
    output:
        vcf = temp(_fold("vcf/{panel}_chr{chrom}.vcf.gz")),
        tbi = temp(_fold("vcf/{panel}_chr{chrom}.vcf.gz.tbi")),
    wildcard_constraints:
        panel = "gt|ref",
    params:
        src = lambda wc, input: input.bed[:-4],
        raw = _fold("vcf/{panel}_chr{chrom}_raw"),
    conda:
        "../envs/workflow_env.yaml"
    log:
        _log + "fold{fold}_{panel}_vcf_chr{chrom}.log",
    shell:
        """
        ({_plink} --bfile {params.src} --chr {wildcards.chrom} --snps-only \
            --export vcf bgz id-paste=iid --out {params.raw} {_plink_extra}
         bcftools norm -d snps {params.raw}.vcf.gz | bgzip > {output.vcf}
         tabix -f -p vcf {output.vcf}
         rm -f {params.raw}.vcf.gz) &> {log}
        """


rule acc_cv_beagle_phase_ref:
    """
    Phase the fold's reference panel; Beagle will not accept it otherwise.

    `ref=` takes "bref3 or VCF file with phased genotypes". A panel exported
    straight from PLINK is unphased, and Beagle aborts the chromosome with
    "unphased or missing genotype for reference sample". Running Beagle on the
    panel alone, with no ref=, phases it and fills its sporadic missing calls.

    This is roughly half of Beagle's CV compute and cannot be shared: every
    fold has a different panel, and phasing the whole cohort once would let
    the held-out animals' HD genotypes shape the panel's haplotypes.
    """
    input:
        vcf = _fold("vcf/ref_chr{chrom}.vcf.gz"),
        tbi = _fold("vcf/ref_chr{chrom}.vcf.gz.tbi"),
    output:
        vcf = temp(_fold("vcf/ref_phased_chr{chrom}.vcf.gz")),
        tbi = temp(_fold("vcf/ref_phased_chr{chrom}.vcf.gz.tbi")),
    params:
        beagle  = config["beagle_jar"],
        window  = config["beagle_params"]["window"],
        overlap = config["beagle_params"]["overlap"],
        ne      = config["beagle_params"]["ne"],
        outbase = _fold("vcf/ref_phased_chr{chrom}"),
        heap_mb = java_heap_mb(70000),
    threads:
        config["beagle_params"]["nthreads"]
    conda:
        "../envs/workflow_env.yaml"
    resources:
        mem_mb = 70000,
        slurm_partition = "r7i-ondemand-4xlarge",
    log:
        _log + "fold{fold}_phase_ref_chr{chrom}.log",
    shell:
        """
        (java -Xmx{params.heap_mb}m -jar {params.beagle} \
            gt={input.vcf} \
            window={params.window} overlap={params.overlap} \
            out={params.outbase} nthreads={threads} ne={params.ne} \
            chrom={wildcards.chrom}) &> {log}
        tabix -f {output.vcf}
        """


rule acc_cv_run_beagle:
    input:
        gt      = _fold("vcf/gt_chr{chrom}.vcf.gz"),
        gt_tbi  = _fold("vcf/gt_chr{chrom}.vcf.gz.tbi"),
        ref     = _fold("vcf/ref_phased_chr{chrom}.vcf.gz"),
        ref_tbi = _fold("vcf/ref_phased_chr{chrom}.vcf.gz.tbi"),
    output:
        vcf = temp(_acc_out + "/beagle/fold{fold}/imputed/chr{chrom}.vcf.gz"),
        tbi = temp(_acc_out + "/beagle/fold{fold}/imputed/chr{chrom}.vcf.gz.tbi"),
    params:
        beagle  = config["beagle_jar"],
        window  = config["beagle_params"]["window"],
        overlap = config["beagle_params"]["overlap"],
        ne      = config["beagle_params"]["ne"],
        outbase = _acc_out + "/beagle/fold{fold}/imputed/chr{chrom}",
        heap_mb = java_heap_mb(70000),
    threads:
        config["beagle_params"]["nthreads"]
    conda:
        "../envs/workflow_env.yaml"
    resources:
        mem_mb = 70000,
        slurm_partition = "r7i-ondemand-4xlarge",
    log:
        _log + "beagle_fold{fold}_chr{chrom}.log",
    shell:
        """
        (java -Xmx{params.heap_mb}m -jar {params.beagle} \
            gt={input.gt} ref={input.ref} \
            window={params.window} overlap={params.overlap} \
            out={params.outbase} nthreads={threads} ne={params.ne} \
            chrom={wildcards.chrom}) &> {log}
        tabix -f {output.vcf}
        """


# ---------------------------------------------------------------------------
# AlphaImpute2: the combined fileset, one chromosome at a time
# ---------------------------------------------------------------------------

_ai2 = config.get("alphaimpute2_params", {})


rule acc_cv_alphaimpute2_input:
    input:
        raw = _fold("chrom/chr{chrom}.raw"),
        fam = _fold("chrom/chr{chrom}.fam"),
    output:
        genotypes = temp(_acc_out + "/alphaimpute2/fold{fold}/input/chr{chrom}.genotypes.txt"),
        pedigree  = temp(_acc_out + "/alphaimpute2/fold{fold}/input/chr{chrom}.pedigree.txt"),
    log:
        _log + "alphaimpute2_fold{fold}_input_chr{chrom}.log",
    shell:
        """
        (awk 'NR>1 {{
            printf $2;
            for (i=7; i<=NF; i++) {{
                if ($i == "NA") printf " 9";
                else printf " %d", $i;
            }}
            printf "\\n";
        }}' {input.raw} > {output.genotypes}
        awk '{{print $2, $3, $4}}' {input.fam} > {output.pedigree}) &> {log}
        """


rule acc_cv_run_alphaimpute2:
    input:
        genotypes = _acc_out + "/alphaimpute2/fold{fold}/input/chr{chrom}.genotypes.txt",
        pedigree  = _acc_out + "/alphaimpute2/fold{fold}/input/chr{chrom}.pedigree.txt",
    output:
        genotypes = temp(_acc_out + "/alphaimpute2/fold{fold}/output/chr{chrom}.genotypes"),
    params:
        out_prefix  = _acc_out + "/alphaimpute2/fold{fold}/output/chr{chrom}",
        cycles      = _ai2.get("cycles", 4),
        threshold   = _ai2.get("final_peeling_threshold", 0.1),
        hd_thresh   = _ai2.get("hd_threshold", 0.95),
        length      = _ai2.get("length", 1.0),
        extra_flags = " ".join(
            flag for key, flag in (
                ("ped_only", "-ped_only"),
                ("pop_only", "-pop_only"),
                ("phase_output", "-phase_output"),
            ) if _ai2.get(key, False)
        ),
    threads:
        _ai2.get("maxthreads", 1)
    conda:
        "../envs/alphaimpute2_env.yaml"
    resources:
        mem_mb = 16000,
        slurm_partition = "r7i-ondemand-2xlarge",
    log:
        _log + "alphaimpute2_fold{fold}_chr{chrom}.log",
    shell:
        """
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
        genotypes = _acc_out + "/alphaimpute2/fold{fold}/output/chr{chrom}.genotypes",
        bim       = _fold("chrom/chr{chrom}.bim"),
    output:
        vcf = temp(_acc_out + "/alphaimpute2/fold{fold}/imputed/chr{chrom}.vcf.gz"),
        tbi = temp(_acc_out + "/alphaimpute2/fold{fold}/imputed/chr{chrom}.vcf.gz.tbi"),
    conda:
        "../envs/workflow_env.yaml"
    resources:
        mem_mb = 32000,
        slurm_partition = "r7i-ondemand-2xlarge",
    log:
        _log + "alphaimpute2_fold{fold}_to_vcf_chr{chrom}.log",
    script:
        "../scripts/alphaimpute2_to_vcf.py"


# ---------------------------------------------------------------------------
# FImpute: the combined fileset, one chromosome at a time
# ---------------------------------------------------------------------------

rule acc_cv_prepare_fimpute_inputs:
    input:
        raw = _fold("chrom/chr{chrom}.raw"),
        bim = _fold("chrom/chr{chrom}.bim"),
        fam = _fold("chrom/chr{chrom}.fam"),
    output:
        genos = _acc_out + "/fimpute/fold{fold}/input/chr{chrom}.genos",
        snps  = _acc_out + "/fimpute/fold{fold}/input/chr{chrom}.snps",
        ped   = _acc_out + "/fimpute/fold{fold}/input/chr{chrom}.ped",
        ctrl  = _acc_out + "/fimpute/fold{fold}/input/chr{chrom}.ctrl",
        idmap = _acc_out + "/fimpute/fold{fold}/input/chr{chrom}.id_map.tsv",
    params:
        out_dir  = _acc_out + "/fimpute/fold{fold}/input",
        nthreads = _fimpute_nthreads,
    conda:
        "../envs/workflow_env.yaml"
    log:
        _log + "fimpute_fold{fold}_prepare_chr{chrom}.log",
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
    # FImpute3 is a licensed binary that conda does not install, so no conda:.
    input:
        ctrl = _acc_out + "/fimpute/fold{fold}/input/chr{chrom}.ctrl",
    output:
        imp    = _acc_out + "/fimpute/fold{fold}/input/chr{chrom}/genotypes_imp.txt",
        report = _acc_out + "/fimpute/fold{fold}/input/chr{chrom}/report.txt",
        # FImpute drops markers of its own accord and rewrites the SNP map for
        # what it kept. Read back this file, not the .snps handed in, or the
        # genotype string and the map disagree by however many it removed.
        snpinfo = _acc_out + "/fimpute/fold{fold}/input/chr{chrom}/snp_info.txt",
    params:
        exe = _fimpute_executable,
    threads:
        _fimpute_nthreads
    resources:
        mem_mb = 32000,
        slurm_partition = "r7i-ondemand-2xlarge",
    log:
        _log + "fimpute_fold{fold}_run_chr{chrom}.log",
    shell:
        """
        ({params.exe} {input.ctrl} -o) &> {log}
        """


rule acc_cv_fimpute_to_vcf:
    input:
        imp   = _acc_out + "/fimpute/fold{fold}/input/chr{chrom}/genotypes_imp.txt",
        bim   = _fold("chrom/chr{chrom}.bim"),
        snps  = _acc_out + "/fimpute/fold{fold}/input/chr{chrom}/snp_info.txt",
        idmap = _acc_out + "/fimpute/fold{fold}/input/chr{chrom}.id_map.tsv",
    output:
        vcf = temp(_acc_out + "/fimpute/fold{fold}/imputed/chr{chrom}.vcf.gz"),
        tbi = temp(_acc_out + "/fimpute/fold{fold}/imputed/chr{chrom}.vcf.gz.tbi"),
    conda:
        "../envs/workflow_env.yaml"
    log:
        _log + "fimpute_fold{fold}_to_vcf_chr{chrom}.log",
    shell:
        """
        (python {workflow.basedir}/scripts/fimpute_io.py to-vcf \
            --imputed {input.imp} \
            --bim {input.bim} \
            --snp-info {input.snps} \
            --id-map {input.idmap} \
            --out-vcf {output.vcf}.tmp
        bgzip -c {output.vcf}.tmp > {output.vcf}
        rm -f {output.vcf}.tmp
        tabix -f -p vcf {output.vcf}) &> {log}
        """


# ---------------------------------------------------------------------------
# Every imputer: concatenate, score, aggregate
# ---------------------------------------------------------------------------

rule acc_cv_concat:
    input:
        vcfs = lambda wc: expand(
            _acc_out + "/{imputer}/fold{fold}/imputed/chr{chrom}.vcf.gz",
            imputer=wc.imputer, fold=wc.fold, chrom=acc_get_chromosomes(),
        ),
        # bcftools concat needs the indexes: Beagle writes no ##contig lines,
        # and temp() indexes that nothing names as input are deleted before
        # this rule runs. Leaving them out is what stopped every Beagle CV run.
        tbis = lambda wc: expand(
            _acc_out + "/{imputer}/fold{fold}/imputed/chr{chrom}.vcf.gz.tbi",
            imputer=wc.imputer, fold=wc.fold, chrom=acc_get_chromosomes(),
        ),
    output:
        vcf = _acc_out + "/{imputer}/fold{fold}/imputed/all_chromosomes.vcf.gz",
        tbi = _acc_out + "/{imputer}/fold{fold}/imputed/all_chromosomes.vcf.gz.tbi",
    conda:
        "../envs/workflow_env.yaml"
    threads: 4
    resources:
        mem_mb = 64000,
        slurm_partition = "r7i-ondemand-4xlarge",
    log:
        _log + "{imputer}_fold{fold}_concat.log",
    shell:
        """
        (bcftools concat --output {output.vcf} --output-type z \
            --threads {threads} {input.vcfs}) 2> {log}
        tabix -f {output.vcf}
        """


rule acc_cv_compute_metrics:
    input:
        imputed     = _acc_out + "/{imputer}/fold{fold}/imputed/all_chromosomes.vcf.gz",
        imputed_tbi = _acc_out + "/{imputer}/fold{fold}/imputed/all_chromosomes.vcf.gz.tbi",
        truth       = _fold("truth.vcf.gz"),
        truth_tbi   = _fold("truth.vcf.gz.tbi"),
    output:
        by_snp   = _acc_out + "/{imputer}/fold{fold}/metrics_by_snp.tsv",
        by_maf   = _acc_out + "/{imputer}/fold{fold}/metrics_by_maf_bin.tsv",
        by_indiv = _acc_out + "/{imputer}/fold{fold}/metrics_by_individual.tsv",
        summary  = _acc_out + "/{imputer}/fold{fold}/summary.tsv",
    params:
        out_dir  = _acc_out + "/{imputer}/fold{fold}",
        maf_bins = " ".join(
            str(b) for b in config.get("maf_bins", [0.01, 0.05, 0.1, 0.2, 0.5])
        ),
    conda:
        "../envs/accuracy_env.yaml"
    log:
        _log + "{imputer}_fold{fold}_metrics.log",
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
            imputer=_cv_imputers, fold=_cv_folds,
        ),
    output:
        summary = _acc_out + "/cv_summary.tsv",
        imputer_summary = _acc_out + "/cv_imputer_summary.tsv",
    params:
        imputers = " ".join(_cv_imputers),
        folds    = " ".join(str(fold) for fold in _cv_folds),
    conda:
        "../envs/accuracy_env.yaml"
    log:
        _log + "aggregate.log",
    shell:
        """
        python {workflow.basedir}/scripts/aggregate_cv_metrics.py \
            --root {_acc_out} \
            --imputers {params.imputers} \
            --folds {params.folds} \
            --summary-out {output.summary} \
            --imputer-summary-out {output.imputer_summary} \
            &> {log}
        """


rule acc_cv_snp_reliability:
    """
    Per-marker accuracy across folds -- the table the reference panel is picked from.

    Markers are filtered on their worst fold rather than their mean, because a
    marker that collapses in one fold out of five is not one to build a panel on.
    """
    input:
        by_snp = expand(
            _acc_out + "/{imputer}/fold{fold}/metrics_by_snp.tsv",
            imputer=_cv_imputers, fold=_cv_folds,
        ),
        bim = _cv_hd_bfile + ".bim",
    output:
        table    = _acc_out + "/snp_reliability.tsv",
        reliable = _acc_out + "/reliable_markers.txt",
    params:
        imputers  = " ".join(_cv_imputers),
        folds     = " ".join(str(fold) for fold in _cv_folds),
        threshold = _cv_reliable_r2,
    conda:
        "../envs/accuracy_env.yaml"
    log:
        _log + "snp_reliability.log",
    shell:
        """
        python {workflow.basedir}/scripts/aggregate_snp_reliability.py \
            --root {_acc_out} \
            --imputers {params.imputers} \
            --folds {params.folds} \
            --bim {input.bim} \
            --r2-threshold {params.threshold} \
            --out {output.table} \
            --reliable-out {output.reliable} \
            &> {log}
        """
