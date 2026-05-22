import gzip
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from make_accuracy_cv_setup import assign_folds, choose_ld_panel
from fimpute_io import (
    fimpute_calls_to_vcf_gt,
    write_fimpute_inputs_from_raw,
    write_fimpute_vcf,
)
from aggregate_cv_metrics import aggregate_summaries
from test_dryrun import REPO_ROOT, REAL_BFILE, requires_snakemake


def test_assign_folds_is_balanced_deterministic_and_complete(synth_plink):
    fam = pd.read_csv(
        synth_plink["fam"],
        sep=r"\s+",
        header=None,
        names=["fid", "iid", "pat", "mat", "sex", "phen"],
    )

    first = assign_folds(fam[["fid", "iid"]], n_folds=10, seed=42)
    second = assign_folds(fam[["fid", "iid"]], n_folds=10, seed=42)

    assert first.equals(second)
    assert len(first) == len(fam)
    assert first[["fid", "iid"]].drop_duplicates().shape[0] == len(fam)
    fold_sizes = first["fold"].value_counts().sort_index()
    assert set(fold_sizes.index) == set(range(1, 11))
    assert fold_sizes.max() - fold_sizes.min() <= 1


def test_choose_ld_panel_uses_same_seeded_panel_for_all_folds(synth_plink):
    bim = pd.read_csv(synth_plink["bim"], sep=r"\s+", header=None)

    first = choose_ld_panel(bim, n_snps=10, seed=99)
    second = choose_ld_panel(bim, n_snps=10, seed=99)

    assert first == second
    assert len(first) == 10
    assert set(first).issubset(set(bim[1]))


def test_choose_ld_panel_accepts_explicit_snp_list(tmp_path, synth_plink):
    snp_list = tmp_path / "panel.txt"
    snp_list.write_text("snp_chr1_001\nsnp_chr2_002\n")
    bim = pd.read_csv(synth_plink["bim"], sep=r"\s+", header=None)

    assert choose_ld_panel(bim, n_snps=10, seed=99, snp_list=snp_list) == [
        "snp_chr1_001",
        "snp_chr2_002",
    ]


def test_write_fimpute_inputs_from_plink_raw_uses_ggf2_layout(tmp_path):
    raw = tmp_path / "chr1.raw"
    raw.write_text(
        "FID IID PAT MAT SEX PHEN snp1_A snp2_A snp3_A\n"
        "F1 A 0 0 1 -9 0 1 NA\n"
        "F1 B A 0 2 -9 2 1 0\n"
    )
    bim = tmp_path / "chr1.bim"
    bim.write_text("1\tsnp1\t0\t100\tA\tG\n1\tsnp2\t0\t200\tA\tG\n1\tsnp3\t0\t300\tA\tG\n")
    fam = tmp_path / "chr1.fam"
    fam.write_text("F1 A 0 0 1 -9\nF1 B A 0 2 -9\n")

    outputs = write_fimpute_inputs_from_raw(
        raw_path=raw,
        bim_path=bim,
        fam_path=fam,
        out_dir=tmp_path / "fimpute",
        chrom="1",
        nthreads=2,
    )

    assert outputs["genos"].read_text().splitlines() == [
        "ID\tChip\tGenotypes",
        "1\t1\t015",
        "2\t1\t210",
    ]
    assert outputs["snps"].read_text().splitlines() == [
        "SNP_ID\tChr\tPos\tChip1",
        "snp1\t1\t100\t1",
        "snp2\t1\t200\t2",
        "snp3\t1\t300\t3",
    ]
    assert outputs["ped"].read_text().splitlines() == [
        "ID\tsire ID\tdam ID\tsex",
        "1\t0\t0\tM",
        "2\t1\t0\tF",
    ]
    assert 'nthread=2;' in outputs["ctrl"].read_text()


def test_fimpute_calls_to_vcf_gt_maps_phased_and_missing_codes():
    assert fimpute_calls_to_vcf_gt("0") == "0/0"
    assert fimpute_calls_to_vcf_gt("2") == "1/1"
    assert fimpute_calls_to_vcf_gt("3") == "0/1"
    assert fimpute_calls_to_vcf_gt("4") == "1/0"
    assert fimpute_calls_to_vcf_gt("5") == "./."


def test_write_fimpute_vcf_keeps_real_sample_ids_and_variant_metadata(tmp_path):
    bim = tmp_path / "chr1.bim"
    bim.write_text("1\tsnp1\t0\t100\tA\tG\n1\tsnp2\t0\t200\tC\tT\n")
    id_map = tmp_path / "id_map.tsv"
    id_map.write_text("short_id\tfid\tiid\n1\tF1\tA\n2\tF1\tB\n")
    imp = tmp_path / "genotypes_imp.txt"
    imp.write_text("ID\tChip\tCalls...\n1\t1\t03\n2\t1\t25\n")
    out_vcf = tmp_path / "chr1.vcf"

    write_fimpute_vcf(imp, bim, id_map, out_vcf)

    lines = [line for line in out_vcf.read_text().splitlines() if not line.startswith("##")]
    assert lines[0].split("\t")[-2:] == ["A", "B"]
    assert lines[1].split("\t") == [
        "1", "100", "snp1", "G", "A", ".", "PASS", ".", "GT", "0/0", "1/1"
    ]
    assert lines[2].split("\t")[-2:] == ["0/1", "./."]


def test_aggregate_summaries_adds_fold_and_imputer_columns(tmp_path):
    for imputer, fold, concordance in [
        ("beagle", 1, "0.91"),
        ("beagle", 2, "0.93"),
        ("fimpute", 1, "0.95"),
    ]:
        path = tmp_path / imputer / f"fold{fold}"
        path.mkdir(parents=True)
        (path / "summary.tsv").write_text(
            "metric\tvalue\n"
            f"mean_concordance\t{concordance}\n"
            "mean_r2\t0.8\n"
        )

    long_df, imputer_df = aggregate_summaries(tmp_path, ["beagle", "fimpute"], [1, 2])

    assert set(long_df.columns) >= {"imputer", "fold", "metric", "value"}
    assert long_df.query("imputer == 'beagle' and metric == 'mean_concordance'").shape[0] == 2
    assert imputer_df.query("imputer == 'beagle' and metric == 'mean_concordance'")[
        "mean"
    ].iloc[0] == pytest.approx(0.92)


@requires_snakemake
def test_accuracy_cv_dryrun_resolves_three_imputer_benchmark(tmp_path):
    out_dir = tmp_path / "accuracy_cv"
    result = subprocess.run(
        [
            "snakemake",
            "--snakefile", "Snakefile_accuracy",
            "--dryrun",
            "--quiet",
            "--config",
            "accuracy_mode=kfold_mask_and_impute",
            f"accuracy_output_dir={out_dir}",
            f"bfile={REAL_BFILE}",
            "plink_path=plink2",
            "beagle_jar=fake.jar",
            "cv_n_folds=2",
            "cv_target_n_snps=10",
            "cv_imputers=beagle alphaimpute2 fimpute",
            "fimpute_params_executable=/mnt/efshome/applications/FImpute3/2026/FImpute3",
        ],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"CV dryrun failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
