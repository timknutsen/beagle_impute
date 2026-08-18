"""
Tests for the cross-array K-fold accuracy mode.

Each test here guards a failure that actually happened, or that the previous
design could not have detected:

  * a quarter of one benchmark step's animals had their two typings taken from
    different physical fish, and nothing checked;
  * an allele swap between two exports mirrors every genotype, and allelic r2
    cannot see it because squaring the correlation hides the sign;
  * a marker that collapses in one fold out of five passes on its mean;
  * an omitted --target-snp-list was read as "the current directory is my SNP
    list", because Path("") is PosixPath(".") and therefore truthy.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from aggregate_snp_reliability import combine_runs, summarise_markers  # noqa: E402
from identity_gate import classify  # noqa: E402
from make_accuracy_cv_setup import choose_ld_panel  # noqa: E402
from mask_validation_genotypes import flip_needed  # noqa: E402
from pair_animals import pair, read_fam  # noqa: E402
from shared_marker_list import read_bim, shared_markers  # noqa: E402

REPO_ROOT = Path(__file__).parent.parent
REAL_BFILE = str(REPO_ROOT / "tests" / "data" / "test_salmon")

requires_snakemake = pytest.mark.skipif(
    not shutil.which("snakemake"), reason="snakemake not on PATH"
)


# ── Pairing ──────────────────────────────────────────────────────────────────

def test_animals_pair_on_iid_even_when_fid_differs(tmp_path):
    """
    Per-array exports carry FID=0 and are produced in separate runs, so an
    (FID, IID) match returns nothing and the run fails with no hint why.
    """
    (tmp_path / "ld.fam").write_text("0 fishA 0 0 1 -9\n0 fishB 0 0 2 -9\n")
    (tmp_path / "hd.fam").write_text("FAM1 fishA 0 0 1 -9\nFAM1 fishC 0 0 2 -9\n")

    paired = pair(read_fam(tmp_path / "ld.fam"), read_fam(tmp_path / "hd.fam"))

    assert paired["iid"].tolist() == ["fishA"]
    # The FID comes from the HD side, so downstream --keep files agree with the
    # fileset the truth is cut from.
    assert paired["fid"].tolist() == ["FAM1"]


# ── Shared marker panel ──────────────────────────────────────────────────────

def test_shared_panel_is_the_intersection_of_the_two_bims(synth_two_chip):
    """
    The panel must be LD n HD exactly. An LD-only probe would be handed to the
    imputer as observed and then never scored, because the truth fileset does
    not carry it.
    """
    ld = read_bim(Path(synth_two_chip["ld_bfile"] + ".bim"))
    hd = read_bim(Path(synth_two_chip["hd_bfile"] + ".bim"))

    markers, stats = shared_markers([ld, hd])

    assert sorted(markers) == sorted(synth_two_chip["shared_ids"])
    assert not set(markers) & set(synth_two_chip["ld_only_ids"])
    assert not set(markers) & set(synth_two_chip["hd_only_ids"])
    assert stats["position_mismatch"] == 0


def test_shared_panel_keeps_markers_whose_alleles_are_merely_reversed(synth_two_chip):
    """
    A1/A2 order is an export artefact, not a difference in the marker, and the
    orientation is forced downstream. Dropping these would silently shrink the
    panel by however many the two exports happened to disagree about.
    """
    ld = read_bim(Path(synth_two_chip["ld_bfile"] + ".bim"))
    hd = read_bim(Path(synth_two_chip["hd_bfile"] + ".bim"))

    markers, stats = shared_markers([ld, hd])

    assert stats["allele_mismatch"] == 0
    assert set(synth_two_chip["swapped_ids"]).issubset(set(markers))


def test_shared_panel_drops_a_marker_that_moved(tmp_path):
    """Same ID at a different position is not the same marker."""
    (tmp_path / "a.bim").write_text("1\tm1\t0\t100\tA\tG\n1\tm2\t0\t200\tA\tG\n")
    (tmp_path / "b.bim").write_text("1\tm1\t0\t100\tA\tG\n1\tm2\t0\t999\tA\tG\n")

    markers, stats = shared_markers(
        [read_bim(tmp_path / "a.bim"), read_bim(tmp_path / "b.bim")]
    )

    assert markers == ["m1"]
    assert stats["position_mismatch"] == 1


# ── Identity gate ────────────────────────────────────────────────────────────

def _verdicts(concordances, threshold=0.90):
    pairs = pd.DataFrame(
        {"fid": ["F"] * len(concordances), "iid": list(concordances)}
    )
    metrics = pd.DataFrame(
        {
            "sample": list(concordances),
            "concordance": list(concordances.values()),
            "n_evaluated": [100] * len(concordances),
        }
    )
    return classify(pairs, metrics, threshold)


def test_identity_gate_separates_the_two_modes():
    """
    The distribution is bimodal -- the same animal lands at 0.95-1.0, two
    different salmon near 0.5 -- so the threshold sits in an empty gap.
    """
    verdicts = _verdicts({"same1": 0.99, "same2": 0.96, "other1": 0.52, "other2": 0.47})

    assert verdicts.set_index("iid")["passed"].to_dict() == {
        "same1": True,
        "same2": True,
        "other1": False,
        "other2": False,
    }


def test_identity_gate_fails_an_animal_it_could_not_score():
    """
    An animal with no overlapping calls has not been shown to be the same fish.
    Defaulting it to pass is how an unverified pair reaches the results.
    """
    pairs = pd.DataFrame({"fid": ["F"], "iid": ["ghost"]})
    metrics = pd.DataFrame({"sample": [], "concordance": [], "n_evaluated": []})

    verdicts = classify(pairs, metrics, 0.90)

    assert verdicts["passed"].tolist() == [False]


# ── Allele orientation ───────────────────────────────────────────────────────

def test_flip_is_detected_only_for_a_genuinely_reversed_pair():
    assert flip_needed("m", ("A", "G"), ("A", "G")) is False
    assert flip_needed("m", ("G", "A"), ("A", "G")) is True


def test_incompatible_alleles_raise_rather_than_being_guessed():
    """
    Two different allele pairs at one ID are two different markers. Splicing
    them would mirror or scramble every call there, and allelic r2 cannot see
    it -- only concordance collapses, much later.
    """
    with pytest.raises(ValueError, match="not the same marker"):
        flip_needed("m", ("C", "T"), ("A", "G"))


# ── Marker reliability ───────────────────────────────────────────────────────

def _fold_metrics(rows):
    return pd.DataFrame(
        rows, columns=["run", "imputer", "variant", "allelic_r2", "concordance",
                       "maf", "n_evaluated", "fold"]
    )


def test_reliability_filters_on_the_worst_fold_not_the_mean():
    """
    A marker that collapses in one fold out of five is not one to build a
    reference panel on, and its mean hides exactly that.
    """
    metrics = _fold_metrics([
        ("runA", "beagle", "1:100:A:G", 0.99, 0.99, 0.3, 100, 1),
        ("runA", "beagle", "1:100:A:G", 0.40, 0.60, 0.3, 100, 2),
        ("runA", "beagle", "1:200:A:G", 0.95, 0.98, 0.3, 100, 1),
        ("runA", "beagle", "1:200:A:G", 0.93, 0.97, 0.3, 100, 2),
    ])

    combined = combine_runs(summarise_markers(metrics), ["variant"])
    by_variant = combined.set_index("variant")

    assert by_variant.loc["1:100:A:G", "mean_r2"] == pytest.approx(0.695)
    assert by_variant.loc["1:100:A:G", "min_r2"] == pytest.approx(0.40)
    assert by_variant.loc["1:200:A:G", "min_r2"] == pytest.approx(0.93)


def test_a_marker_missing_from_one_run_is_not_counted_as_tested_in_it():
    """
    Intersecting runs is the point of passing --root twice: a V4 marker earns
    its place by holding up in both the real array transition and the masked
    test, which stress different things.
    """
    metrics = _fold_metrics([
        ("runA", "beagle", "1:100:A:G", 0.99, 0.99, 0.3, 100, 1),
        ("runB", "beagle", "1:100:A:G", 0.98, 0.99, 0.3, 100, 1),
        ("runA", "beagle", "1:200:A:G", 0.99, 0.99, 0.3, 100, 1),
    ])

    combined = combine_runs(summarise_markers(metrics), ["variant"]).set_index("variant")

    assert combined.loc["1:100:A:G", "n_runs"] == 2
    assert combined.loc["1:200:A:G", "n_runs"] == 1


# ── CV setup ─────────────────────────────────────────────────────────────────

def test_an_omitted_snp_list_samples_instead_of_reading_the_cwd(synth_plink):
    """
    choose_ld_panel used to be handed Path("") -- which is PosixPath("."), and
    truthy -- so every random-panel run died trying to read the working
    directory as a SNP list. Dry-run cannot catch it.
    """
    bim = pd.read_csv(synth_plink["bim"], sep=r"\s+", header=None)

    panel = choose_ld_panel(bim, n_snps=10, seed=1, snp_list=None)

    assert len(panel) == 10
    assert set(panel).issubset(set(bim[1].astype(str)))


# ── DAG ──────────────────────────────────────────────────────────────────────

@requires_snakemake
def test_cross_array_dryrun_resolves_with_the_identity_gate(tmp_path, synth_two_chip):
    """The mode had never produced a valid number; the DAG must at least build."""
    out_dir = tmp_path / "accuracy_cross"
    result = subprocess.run(
        [
            "snakemake",
            "--snakefile", "Snakefile_accuracy",
            "--dryrun", "--quiet", "rules",
            "--config",
            "accuracy_mode=cross_array",
            f"accuracy_output_dir={out_dir}",
            f"cross_array_ld_bfile={synth_two_chip['ld_bfile']}",
            f"cross_array_hd_bfile={synth_two_chip['hd_bfile']}",
            "cross_array_min_shared_markers=10",
            "plink_path=plink2",
            "beagle_jar=fake.jar",
            "cv_n_folds=2",
            "cv_imputers=beagle alphaimpute2 fimpute",
        ],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, (
        f"cross_array dryrun failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
    for rule in (
        "acc_cv_pair_animals",
        "acc_cv_shared_markers",
        "acc_cv_identity_gate",
        "acc_cv_snp_reliability",
    ):
        assert rule in result.stdout, f"{rule} missing from the DAG:\n{result.stdout}"


@requires_snakemake
def test_flat_config_aliases_override_the_yaml_block(tmp_path):
    """
    --config cv_n_folds=2 used to be accepted and ignored: _nested_config read
    the YAML block first, and config_accuracy.yaml defines every cv.* key. The
    run then silently used 10 folds while reporting the command that asked for 2.
    """
    out_dir = tmp_path / "accuracy_cv_folds"
    result = subprocess.run(
        [
            "snakemake",
            "--snakefile", "Snakefile_accuracy",
            "--dryrun", "--quiet", "rules",
            "--config",
            "accuracy_mode=kfold_mask_and_impute",
            f"accuracy_output_dir={out_dir}",
            f"bfile={REAL_BFILE}",
            "plink_path=plink2",
            "beagle_jar=fake.jar",
            "cv_n_folds=2",
            "cv_target_n_snps=10",
            "cv_imputers=beagle",
        ],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    fold_jobs = [
        line for line in result.stdout.splitlines()
        if line.startswith("acc_cv_fold_ids")
    ]
    assert fold_jobs, f"acc_cv_fold_ids missing:\n{result.stdout}"
    assert fold_jobs[0].split()[-1] == "2", fold_jobs[0]


# ── Cross-array combined fileset ─────────────────────────────────────────────

_BED_DECODE = {0b00: 2, 0b01: None, 0b10: 1, 0b11: 0}


def _read_bed_dosages(prefix):
    """Return {iid: {snp_id: a1_dosage or None}} from a PLINK fileset."""
    prefix = Path(prefix)
    iids = [line.split()[1] for line in Path(f"{prefix}.fam").read_text().splitlines()]
    snps = [line.split()[1] for line in Path(f"{prefix}.bim").read_text().splitlines()]
    raw = Path(f"{prefix}.bed").read_bytes()
    assert raw[:3] == b"\x6c\x1b\x01"

    bytes_per_variant = (len(iids) + 3) // 4
    out = {iid: {} for iid in iids}
    for variant_index, snp in enumerate(snps):
        start = 3 + variant_index * bytes_per_variant
        block = raw[start:start + bytes_per_variant]
        for sample_index, iid in enumerate(iids):
            byte_index, offset = divmod(sample_index, 4)
            out[iid][snp] = _BED_DECODE[(block[byte_index] >> (2 * offset)) & 0b11]
    return out


def _run_mask(args):
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "mask_validation_genotypes.py")] + args,
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    return dict(
        line.split("=", 1) for line in result.stdout.strip().splitlines() if "=" in line
    )


def test_replace_from_takes_the_low_density_array_own_calls(synth_two_chip, tmp_path):
    """
    The held-out animals' panel genotypes must come from the array they were
    really typed on. Reusing the high-density calls at the same positions
    deletes the probe differences and between-array error the run exists to
    measure, turning a cross-array benchmark into a masking benchmark.
    """
    fold = [animal[1] for animal in synth_two_chip["animals"][:6]]
    ids = tmp_path / "fold.ids"
    ids.write_text("".join(f"FAM1\t{iid}\n" for iid in fold))
    panel = tmp_path / "panel.snps"
    panel.write_text("\n".join(synth_two_chip["shared_ids"]) + "\n")
    out = tmp_path / "combined"

    _run_mask([
        "--bfile", synth_two_chip["hd_bfile"],
        "--validation-ids", str(ids),
        "--ld-snps", str(panel),
        "--replace-from", synth_two_chip["ld_bfile"],
        "--out", str(out),
    ])

    combined = _read_bed_dosages(out)
    hd = _read_bed_dosages(synth_two_chip["hd_bfile"])
    truth = synth_two_chip["truth"]

    # Held-out animals: panel markers carry the LD array's call, everything
    # else is missing.
    disagreements = 0
    for iid in fold:
        index = [a[1] for a in synth_two_chip["animals"]].index(iid)
        for snp in synth_two_chip["shared_ids"]:
            assert combined[iid][snp] is not None
            if combined[iid][snp] != truth[snp][index]:
                disagreements += 1
        for snp in synth_two_chip["hd_only_ids"]:
            assert combined[iid][snp] is None, f"{iid} kept an HD-only call at {snp}"

    # The fixture plants between-array noise plus two mispaired animals, so the
    # spliced calls must differ from the HD truth somewhere -- otherwise the
    # splice silently did nothing.
    assert disagreements > 0

    # Everyone else is untouched at full density.
    for animal in synth_two_chip["animals"][6:]:
        iid = animal[1]
        assert combined[iid] == hd[iid]


def test_swapped_markers_are_oriented_before_they_are_spliced(synth_two_chip, tmp_path):
    """
    Three fixture markers have A1/A2 reversed on the low-density chip. Splicing
    their raw codes would mirror every call there, and allelic r2 would not
    show it.
    """
    fold = [animal[1] for animal in synth_two_chip["animals"][:6]]
    ids = tmp_path / "fold.ids"
    ids.write_text("".join(f"FAM1\t{iid}\n" for iid in fold))
    panel = tmp_path / "panel.snps"
    panel.write_text("\n".join(synth_two_chip["shared_ids"]) + "\n")
    out = tmp_path / "combined"

    _run_mask([
        "--bfile", synth_two_chip["hd_bfile"],
        "--validation-ids", str(ids),
        "--ld-snps", str(panel),
        "--replace-from", synth_two_chip["ld_bfile"],
        "--out", str(out),
    ])

    combined = _read_bed_dosages(out)
    all_iids = [animal[1] for animal in synth_two_chip["animals"]]
    truth = synth_two_chip["truth"]

    # The two mispaired animals are random by construction, so judge on the rest.
    honest = [iid for iid in fold if iid not in synth_two_chip["mispaired_iids"]]
    for snp in synth_two_chip["swapped_ids"]:
        matches = sum(
            combined[iid][snp] == truth[snp][all_iids.index(iid)] for iid in honest
        )
        assert matches >= len(honest) - 1, (
            f"{snp} looks mirrored: only {matches}/{len(honest)} calls survived the splice"
        )


def test_a_panel_marker_the_ld_array_lacks_is_masked_not_left_as_truth(
    synth_two_chip, tmp_path
):
    """
    Handing back the high-density call for a marker the low-density array never
    measured gives the imputer an answer it was never given in the real
    transition.
    """
    fold = [animal[1] for animal in synth_two_chip["animals"][:6]]
    ids = tmp_path / "fold.ids"
    ids.write_text("".join(f"FAM1\t{iid}\n" for iid in fold))
    panel = tmp_path / "panel.snps"
    # Deliberately include an HD-only marker in the "panel".
    intruder = synth_two_chip["hd_only_ids"][0]
    panel.write_text("\n".join(synth_two_chip["shared_ids"] + [intruder]) + "\n")
    out = tmp_path / "combined"

    counts = _run_mask([
        "--bfile", synth_two_chip["hd_bfile"],
        "--validation-ids", str(ids),
        "--ld-snps", str(panel),
        "--replace-from", synth_two_chip["ld_bfile"],
        "--out", str(out),
    ])

    assert int(counts["panel_variants_without_replacement"]) == 1
    combined = _read_bed_dosages(out)
    assert all(combined[iid][intruder] is None for iid in fold)


# ── plink2 --export A counts a2 ──────────────────────────────────────────────

def test_fimpute_encoder_reads_the_counted_allele_off_the_raw_column(
    synth_raw_counting_a2, synth_two_chip
):
    """
    Our exports count a2, not a1. Assuming otherwise mirrors every genotype,
    which allelic r2 cannot see -- it showed up once as concordance 0.41
    against a perfectly reasonable-looking r2 of 0.932.
    """
    from fimpute_io import _panel_genotype_strings

    records = synth_raw_counting_a2["records"]
    raw = pd.read_csv(synth_raw_counting_a2["raw"], sep=r"\s+")
    bim = pd.DataFrame(
        [(c, s, p, a1, a2) for c, s, p, a1, a2 in records],
        columns=["chrom", "snp", "pos", "a1", "a2"],
    )

    encoded = _panel_genotype_strings(raw, bim)

    ordered = bim.sort_values("pos")["snp"].tolist()
    truth = synth_two_chip["truth"]
    first_animal = encoded[0]
    # FImpute codes are the a1 dosage; the .raw counted a2, so the encoder must
    # have flipped every column back.
    assert [int(ch) for ch in first_animal] == [truth[snp][0] for snp in ordered]


def test_a_duplicate_physical_position_is_dropped_across_both_chips(synth_two_chip):
    """
    FImpute aborts the whole run on these ("SNPs with the same physical
    position found"), and real arrays do carry them: a redesigned probe keeps
    both IDs. The old fixture had none, so nothing exercised the guard.
    """
    from fimpute_io import drop_duplicate_positions

    def as_frame(records):
        return pd.DataFrame(
            [(c, s, p, a1, a2) for c, s, p, a1, a2 in records],
            columns=["chrom", "snp", "pos", "a1", "a2"],
        )

    ld = as_frame(synth_two_chip["ld_records"])
    hd = as_frame(synth_two_chip["hd_records"])
    colliding = synth_two_chip["duplicate_position_ids"]
    assert (
        hd.loc[hd["snp"] == colliding[0], "pos"].iloc[0]
        == hd.loc[hd["snp"] == colliding[1], "pos"].iloc[0]
    )

    deduped_ld, deduped_hd = drop_duplicate_positions([ld, hd])

    survivors = set(deduped_hd["snp"]) & set(colliding)
    assert len(survivors) == 1, f"both markers survived: {survivors}"
    assert not deduped_hd.duplicated(subset=["chrom", "pos"]).any()
    assert not deduped_ld.duplicated(subset=["chrom", "pos"]).any()


# ── Reference panel selection ────────────────────────────────────────────────

from select_refpanel import family_key, select  # noqa: E402


def _panel_fam(n_families=3, per_family=5):
    rows = []
    for f in range(n_families):
        for i in range(per_family):
            rows.append(("FAM1", f"fish_{f}_{i}", f"sire{f}", f"dam{f}"))
    return (
        pd.DataFrame(rows, columns=["fid", "iid", "sire", "dam"]),
        pd.DataFrame(rows, columns=["fid", "iid", "sire", "dam"])[["iid", "sire", "dam"]],
    )


def test_family_cap_thins_each_family_and_keeps_them_all():
    """
    The tenth full sib adds almost no haplotype a panel can use, but costs the
    same phasing time and pulls allele frequencies toward the biggest families.
    Capping must thin every family, not drop whole ones.
    """
    fam, ped = _panel_fam(n_families=3, per_family=5)

    panel = select(fam, ped, excluded=set(), max_per_family=2, seed=1)
    kept = panel.loc[panel["kept"]]

    assert len(kept) == 6
    assert kept.groupby(["sire", "dam"]).size().tolist() == [2, 2, 2]


def test_family_cap_is_deterministic_at_a_seed():
    fam, ped = _panel_fam()
    first = select(fam, ped, set(), 2, seed=7)
    second = select(fam, ped, set(), 2, seed=7)
    assert first.loc[first.kept, "iid"].tolist() == second.loc[second.kept, "iid"].tolist()


def test_animals_without_parents_are_not_capped_as_one_family():
    """
    Founders carry no sire or dam. Treating "unknown" as a single family would
    cap the entire founder set to N animals and throw away the broadest
    haplotype diversity in the panel.
    """
    fam = pd.DataFrame(
        [("FAM1", f"founder{i}", "0", "0") for i in range(10)],
        columns=["fid", "iid", "sire", "dam"],
    )
    ped = fam[["iid", "sire", "dam"]]

    panel = select(fam, ped, excluded=set(), max_per_family=2, seed=1)

    assert panel["kept"].all(), "founders were capped as though they were full sibs"
    assert family_key({"iid": "a", "sire": "0", "dam": "0"}) != family_key(
        {"iid": "b", "sire": "0", "dam": "0"}
    )


def test_exclusions_are_applied_before_the_cap_and_recorded():
    """
    Held-out animals must leave before phasing. Subsetting a phased panel later
    is cheaper but leaks: their genotypes would already have informed the phase
    of the relatives that remain.
    """
    fam, ped = _panel_fam(n_families=2, per_family=4)
    excluded = {"fish_0_0", "fish_0_1"}

    panel = select(fam, ped, excluded=excluded, max_per_family=4, seed=1)

    assert set(panel.loc[~panel["kept"], "iid"]) == excluded
    reasons = dict(zip(panel["iid"], panel["reason"]))
    assert reasons["fish_0_0"] == "held out by exclusion list"
    assert reasons["fish_1_0"] == "selected"


@requires_snakemake
def test_refpanel_dryrun_resolves(tmp_path, synth_two_chip):
    """The panel workflow is the piece that lets reference_vcf ever be phased."""
    result = subprocess.run(
        [
            "snakemake",
            "--snakefile", "Snakefile_refpanel",
            "--dryrun", "--quiet", "rules",
            "--config",
            f"refpanel_bfile={synth_two_chip['hd_bfile']}",
            "refpanel_name=testarray",
            "refpanel_max_per_family=0",
            f"output_dir={tmp_path / 'refpanel'}",
            "plink_path=plink2",
            "beagle_jar=fake.jar",
            "bref3_jar=fake_bref3.jar",
        ],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"refpanel dryrun failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
    )
    for rule in ("refpanel_select", "refpanel_phase", "refpanel_bref3", "refpanel_report"):
        assert rule in result.stdout, f"{rule} missing:\n{result.stdout}"
