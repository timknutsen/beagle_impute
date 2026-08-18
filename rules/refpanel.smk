# rules/refpanel.smk — build a phased, bref3-encoded reference panel.
#
# Included exclusively by Snakefile_refpanel.
#
#   select → keep → per-chrom VCF → normalize → PHASE → bref3
#
# The phasing step is the one the rest of the repo was missing. config.yaml
# documents reference_vcf as "must be phased" and trusts the caller; Beagle
# enforces it, aborting a chromosome with "unphased or missing genotype for
# reference sample <id>" the moment it is handed a plain PLINK export.

_rp          = config.get("refpanel", {})
_rp_bfile    = str(_rp.get("bfile", "")).strip()
_rp_name     = str(_rp.get("name", "panel")).strip() or "panel"
_rp_pedigree = str(_rp.get("pedigree", "")).strip()
_rp_max_fam  = int(_rp.get("max_per_family", 0))
_rp_seed     = int(_rp.get("random_seed", 42))
_rp_variants = _rp.get("variants", {"full": []}) or {"full": []}
_rp_out      = config.get("output_dir", "refpanel_output")

if not _rp_bfile:
    raise ValueError(
        "refpanel.bfile is required: the QC'd PLINK prefix for ONE array "
        "(or the flat alias refpanel_bfile=)."
    )


def _rp_excludes(variant):
    """Exclusion files for one variant; a bare string is one file, not a list of chars."""
    value = _rp_variants.get(variant, [])
    if isinstance(value, str):
        return [value] if value.strip() else []
    return [str(v) for v in (value or []) if str(v).strip()]



rule refpanel_select:
    """
    Decide who is in the panel, and record why for everyone.

    Family capping needs a pedigree; an array export has none in its .fam, so
    the script refuses rather than silently building an uncapped panel.
    """
    input:
        fam = _rp_bfile + ".fam",
        pedigree = ([_rp_pedigree] if _rp_pedigree else []),
        excludes = lambda wc: _rp_excludes(wc.variant),
    output:
        ids      = _rp_out + "/{variant}/panel_ids.txt",
        manifest = _rp_out + "/{variant}/manifest.tsv",
    params:
        pedigree = _rp_pedigree,
        max_fam  = _rp_max_fam,
        seed     = _rp_seed,
        excludes = lambda wc: " ".join(_rp_excludes(wc.variant)),
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/refpanel/{variant}_select.log",
    shell:
        """
        python {workflow.basedir}/scripts/select_refpanel.py \
            --fam {input.fam} \
            --pedigree "{params.pedigree}" \
            --exclude {params.excludes} \
            --max-per-family {params.max_fam} \
            --seed {params.seed} \
            --ids-out {output.ids} \
            --manifest-out {output.manifest} \
            &> {log}
        """


rule refpanel_bfile:
    input:
        bed = _rp_bfile + ".bed",
        ids = _rp_out + "/{variant}/panel_ids.txt",
    output:
        bed = temp(_rp_out + "/{variant}/panel.bed"),
        bim = temp(_rp_out + "/{variant}/panel.bim"),
        fam = temp(_rp_out + "/{variant}/panel.fam"),
    params:
        plink = config["plink_path"],
        bfile = _rp_bfile,
        out   = _rp_out + "/{variant}/panel",
        extra = config.get("plink_extra_flags", ""),
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/refpanel/{variant}_bfile.log",
    shell:
        """
        ({params.plink} \
            --bfile {params.bfile} \
            --keep {input.ids} \
            --make-bed \
            --out {params.out} \
            --nonfounders --allow-no-sex {params.extra}) &> {log}
        """


rule refpanel_chrom_vcf:
    input:
        # All three, not just the .bed: the trio is temp(), so anything not
        # named here is swept the moment refpanel_bfile finishes and plink2
        # then fails on a missing .fam.
        bed = _rp_out + "/{variant}/panel.bed",
        bim = _rp_out + "/{variant}/panel.bim",
        fam = _rp_out + "/{variant}/panel.fam",
    output:
        vcf = temp(_rp_out + "/{variant}/raw/chr{chrom}.vcf.gz"),
    params:
        plink = config["plink_path"],
        bfile = _rp_out + "/{variant}/panel",
        out   = _rp_out + "/{variant}/raw/chr{chrom}",
        extra = config.get("plink_extra_flags", ""),
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/refpanel/{variant}_export_chr{chrom}.log",
    shell:
        """
        ({params.plink} \
            --bfile {params.bfile} \
            --chr {wildcards.chrom} \
            --snps-only \
            --export vcf bgz id-paste=iid \
            --out {params.out} \
            --nonfounders --allow-no-sex {params.extra}) &> {log}
        """


rule refpanel_normalize:
    """
    Drop the duplicate physical positions before Beagle sees them.

    Redesigned probes keep both IDs at one coordinate, and real arrays carry
    them -- 31 in step 1's high-density fileset alone.
    """
    input:
        vcf = _rp_out + "/{variant}/raw/chr{chrom}.vcf.gz",
    output:
        vcf = temp(_rp_out + "/{variant}/normalized/chr{chrom}.vcf.gz"),
        tbi = temp(_rp_out + "/{variant}/normalized/chr{chrom}.vcf.gz.tbi"),
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/refpanel/{variant}_normalize_chr{chrom}.log",
    shell:
        """
        (bcftools norm -d snps {input.vcf} | bgzip > {output.vcf}) 2> {log}
        tabix -f -p vcf {output.vcf}
        """


rule refpanel_phase:
    """
    Phase the panel and fill its missing calls -- Beagle needs both.

    Run with gt= and no ref=, so the panel is phased against itself. This is
    the expensive step and the whole reason a panel is worth building once:
    every later imputation reads the bref3 instead of repeating it.
    """
    input:
        vcf = _rp_out + "/{variant}/normalized/chr{chrom}.vcf.gz",
        tbi = _rp_out + "/{variant}/normalized/chr{chrom}.vcf.gz.tbi",
    output:
        vcf = _rp_out + "/{variant}/phased/chr{chrom}.vcf.gz",
        tbi = _rp_out + "/{variant}/phased/chr{chrom}.vcf.gz.tbi",
    params:
        beagle  = config["beagle_jar"],
        window  = config["beagle_params"]["window"],
        overlap = config["beagle_params"]["overlap"],
        ne      = config["beagle_params"]["ne"],
        outbase = _rp_out + "/{variant}/phased/chr{chrom}",
        heap_mb = java_heap_mb(70000),
    threads:
        config["beagle_params"]["nthreads"]
    conda:
        "../envs/workflow_env.yaml"
    resources:
        mem_mb = 70000,
        slurm_partition = "r7i-ondemand-4xlarge",
    log:
        "logs/refpanel/{variant}_phase_chr{chrom}.log",
    shell:
        """
        (java -Xmx{params.heap_mb}m -jar {params.beagle} \
            gt={input.vcf} \
            window={params.window} overlap={params.overlap} \
            out={params.outbase} nthreads={threads} ne={params.ne} \
            chrom={wildcards.chrom}) &> {log}
        tabix -f {output.vcf}
        """


rule refpanel_bref3:
    """Beagle's binary reference format; loads 3-40x faster than the VCF."""
    input:
        vcf = _rp_out + "/{variant}/phased/chr{chrom}.vcf.gz",
    output:
        bref3 = _rp_out + "/{variant}/bref3/chr{chrom}.bref3",
    params:
        bref3_jar = config["bref3_jar"],
        heap_mb   = java_heap_mb(16000),
    conda:
        "../envs/workflow_env.yaml"
    resources:
        mem_mb = 16000,
        slurm_partition = "r7i-ondemand-2xlarge",
    log:
        "logs/refpanel/{variant}_bref3_chr{chrom}.log",
    shell:
        """
        (java -Xmx{params.heap_mb}m -jar {params.bref3_jar} \
            {input.vcf} > {output.bref3}) 2> {log}
        """


rule refpanel_report:
    """
    One file that says what this panel is, so it can be explained a year later.

    A bref3 blob carries no provenance: nothing in it records which array it
    came from, who was held out, or which selection produced it.
    """
    input:
        manifest = _rp_out + "/{variant}/manifest.tsv",
        phased   = lambda wc: expand(
            _rp_out + "/{variant}/phased/chr{chrom}.vcf.gz",
            variant=wc.variant, chrom=rp_get_chromosomes(),
        ),
        bref3    = lambda wc: expand(
            _rp_out + "/{variant}/bref3/chr{chrom}.bref3",
            variant=wc.variant, chrom=rp_get_chromosomes(),
        ),
    output:
        report = _rp_out + "/{variant}/panel_report.tsv",
    params:
        name     = _rp_name,
        source   = _rp_bfile,
        pedigree = _rp_pedigree or "(none)",
        max_fam  = _rp_max_fam,
        seed     = _rp_seed,
        excludes = lambda wc: ",".join(_rp_excludes(wc.variant)) or "(none)",
    conda:
        "../envs/workflow_env.yaml"
    log:
        "logs/refpanel/{variant}_report.log",
    run:
        import pandas as pd

        manifest = pd.read_csv(input.manifest, sep="\t", dtype=str)
        kept = manifest[manifest["kept"].astype(str).str.lower() == "true"]
        markers = 0
        for vcf in input.phased:
            markers += sum(1 for _ in shell(f"bcftools view -H {vcf}", iterable=True))

        rows = [
            ("panel_name", params.name),
            ("variant", wildcards.variant),
            ("source_bfile", params.source),
            ("pedigree", params.pedigree),
            ("max_per_family", params.max_fam),
            ("random_seed", params.seed),
            ("exclusion_lists", params.excludes),
            ("animals_in_source", len(manifest)),
            ("animals_in_panel", len(kept)),
            ("chromosomes", len(input.phased)),
            ("markers_phased", markers),
        ]
        pd.DataFrame(rows, columns=["key", "value"]).to_csv(
            output.report, sep="\t", index=False
        )
