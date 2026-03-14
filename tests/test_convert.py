"""
Unit tests for AlphaImpute2 format conversion helpers.

Tests are split into two groups:
  - Pure-Python (no external tools required): always run.
  - bgzip/tabix integration (requires htslib): marked with @pytest.mark.slow
    and skipped automatically when bgzip is not on PATH.
"""

import gzip
import shutil
import subprocess

import pytest

from alphaimpute2_to_vcf import GT_MAP, read_ai2_genotypes, read_bim

# Skip the integration tests when bgzip isn't installed
_bgzip_available = shutil.which("bgzip") is not None
_tabix_available = shutil.which("tabix") is not None
_htslib_available = _bgzip_available and _tabix_available

requires_htslib = pytest.mark.skipif(
    not _htslib_available,
    reason="bgzip/tabix (htslib) not on PATH",
)


# ── GT_MAP ────────────────────────────────────────────────────────────────────

class TestGtMap:
    def test_hom_ref(self):
        assert GT_MAP["0"] == "0/0"

    def test_het(self):
        assert GT_MAP["1"] == "0/1"

    def test_hom_alt(self):
        assert GT_MAP["2"] == "1/1"

    def test_missing(self):
        assert GT_MAP["9"] == "./."

    def test_covers_all_valid_codes(self):
        assert set(GT_MAP.keys()) == {"0", "1", "2", "9"}

    def test_no_unknown_gt_string(self):
        """All GT values must be valid VCF genotype strings."""
        valid = {"0/0", "0/1", "1/1", "./."}
        assert set(GT_MAP.values()) == valid


# ── read_bim ─────────────────────────────────────────────────────────────────

class TestReadBim:
    def test_record_count(self, synth_plink):
        snps = read_bim(synth_plink["bim"])
        expected = len(synth_plink["snp_records"])
        assert len(snps) == expected

    def test_tuple_length(self, synth_plink):
        snps = read_bim(synth_plink["bim"])
        assert all(len(s) == 5 for s in snps)

    def test_chrom_values(self, synth_plink):
        snps = read_bim(synth_plink["bim"])
        chroms = {s[0] for s in snps}
        assert chroms == {"1", "2"}

    def test_positions_are_positive_ints(self, synth_plink):
        snps = read_bim(synth_plink["bim"])
        assert all(isinstance(s[2], int) and s[2] > 0 for s in snps)

    def test_allele_columns(self, synth_plink):
        """a1 (ALT) and a2 (REF) must be non-empty strings."""
        snps = read_bim(synth_plink["bim"])
        assert all(s[3] and s[4] for s in snps)

    def test_alleles_match_fixture(self, synth_plink):
        snps = read_bim(synth_plink["bim"])
        # Fixture uses A=ALT, G=REF throughout
        assert all(s[3] == "A" and s[4] == "G" for s in snps)


# ── read_ai2_genotypes ────────────────────────────────────────────────────────

class TestReadAi2Genotypes:
    def test_individual_count(self, synth_ai2):
        indivs, rows = read_ai2_genotypes(synth_ai2["genotypes"])
        assert len(indivs) == len(synth_ai2["indivs"])

    def test_row_count(self, synth_ai2, synth_plink):
        indivs, rows = read_ai2_genotypes(synth_ai2["genotypes"])
        assert len(rows) == len(synth_plink["indivs"])

    def test_snp_count_per_row(self, synth_ai2, synth_plink):
        indivs, rows = read_ai2_genotypes(synth_ai2["genotypes"])
        expected_snps = len(synth_plink["snp_records"])
        assert all(len(r) == expected_snps for r in rows)

    def test_iid_order(self, synth_ai2):
        indivs, _ = read_ai2_genotypes(synth_ai2["genotypes"])
        expected_iids = [ind[1] for ind in synth_ai2["indivs"]]
        assert indivs == expected_iids

    def test_valid_genotype_codes(self, synth_ai2):
        _, rows = read_ai2_genotypes(synth_ai2["genotypes"])
        valid = {"0", "1", "2", "9"}
        for row in rows:
            assert set(row) <= valid, f"Unexpected codes: {set(row) - valid}"

    def test_founders_have_no_missing(self, synth_ai2):
        """The 4 founder individuals must be fully genotyped (no 9s)."""
        indivs, rows = read_ai2_genotypes(synth_ai2["genotypes"])
        founder_iids = {"sire1", "sire2", "dam1", "dam2"}
        for iid, row in zip(indivs, rows):
            if iid in founder_iids:
                assert "9" not in row, f"Founder {iid} has missing genotypes"

    def test_offspring_have_some_missing(self, synth_ai2):
        """At least one offspring row must contain a missing genotype."""
        indivs, rows = read_ai2_genotypes(synth_ai2["genotypes"])
        offspring_rows = [row for iid, row in zip(indivs, rows) if iid.startswith("off")]
        assert any("9" in row for row in offspring_rows)

    def test_roundtrip_with_fixture_matrix(self, synth_ai2, synth_plink):
        """Values read back must match the original synthetic matrix."""
        indivs, rows = read_ai2_genotypes(synth_ai2["genotypes"])
        original = synth_ai2["geno_by_indiv"]
        for s_idx, (row, orig_row) in enumerate(zip(rows, original)):
            for snp_idx, (code, orig_g) in enumerate(zip(row, orig_row)):
                assert code == str(orig_g), (
                    f"Mismatch at sample {s_idx}, SNP {snp_idx}: "
                    f"got {code!r}, expected {orig_g!r}"
                )


# ── Pedigree file format ──────────────────────────────────────────────────────

class TestPedigreeFile:
    def test_line_count(self, synth_ai2, synth_plink):
        with open(synth_ai2["pedigree"]) as fh:
            lines = [l.strip() for l in fh if l.strip()]
        assert len(lines) == len(synth_plink["indivs"])

    def test_three_columns(self, synth_ai2):
        with open(synth_ai2["pedigree"]) as fh:
            for line in fh:
                parts = line.strip().split()
                assert len(parts) == 3, f"Expected 3 columns, got: {line!r}"

    def test_founders_have_zero_parents(self, synth_ai2):
        with open(synth_ai2["pedigree"]) as fh:
            rows = [l.strip().split() for l in fh if l.strip()]
        founders = {r[0]: (r[1], r[2]) for r in rows if r[0] in {"sire1", "sire2", "dam1", "dam2"}}
        assert len(founders) == 4
        for iid, (pat, mat) in founders.items():
            assert pat == "0" and mat == "0", f"Founder {iid} has non-zero parent"

    def test_offspring_have_nonzero_parents(self, synth_ai2):
        with open(synth_ai2["pedigree"]) as fh:
            rows = [l.strip().split() for l in fh if l.strip()]
        offspring = [r for r in rows if r[0].startswith("off")]
        assert len(offspring) == 8
        for iid, pat, mat in offspring:
            assert pat != "0" and mat != "0", f"Offspring {iid} has unknown parent"


# ── Full convert() round-trip (requires bgzip + tabix) ───────────────────────

@requires_htslib
class TestConvertRoundtrip:
    """End-to-end test: AI2 genotypes + BIM → bgzipped VCF → verify content."""

    def test_vcf_is_valid_bgzip(self, synth_ai2, synth_plink, tmp_path):
        from alphaimpute2_to_vcf import convert

        out_vcf = str(tmp_path / "out.vcf.gz")
        convert(synth_ai2["genotypes"], synth_plink["bim"], out_vcf)

        # bgzip-compressed files start with the gzip magic bytes
        with open(out_vcf, "rb") as fh:
            magic = fh.read(2)
        assert magic == b"\x1f\x8b", "Output is not gzip-compressed"

    def test_tbi_index_created(self, synth_ai2, synth_plink, tmp_path):
        from alphaimpute2_to_vcf import convert

        out_vcf = str(tmp_path / "out.vcf.gz")
        convert(synth_ai2["genotypes"], synth_plink["bim"], out_vcf)
        assert (tmp_path / "out.vcf.gz.tbi").exists()

    def test_vcf_has_correct_sample_count(self, synth_ai2, synth_plink, tmp_path):
        from alphaimpute2_to_vcf import convert

        out_vcf = str(tmp_path / "out.vcf.gz")
        convert(synth_ai2["genotypes"], synth_plink["bim"], out_vcf)

        header_line = _read_vcf_header(out_vcf)
        samples = header_line.split("\t")[9:]
        assert len(samples) == len(synth_plink["indivs"])

    def test_vcf_has_correct_variant_count(self, synth_ai2, synth_plink, tmp_path):
        from alphaimpute2_to_vcf import convert

        out_vcf = str(tmp_path / "out.vcf.gz")
        convert(synth_ai2["genotypes"], synth_plink["bim"], out_vcf)

        data_lines = _read_vcf_data_lines(out_vcf)
        assert len(data_lines) == len(synth_plink["snp_records"])

    def test_vcf_ref_alt_alleles(self, synth_ai2, synth_plink, tmp_path):
        """REF must be A2 (G) and ALT must be A1 (A) from the BIM fixture."""
        from alphaimpute2_to_vcf import convert

        out_vcf = str(tmp_path / "out.vcf.gz")
        convert(synth_ai2["genotypes"], synth_plink["bim"], out_vcf)

        for line in _read_vcf_data_lines(out_vcf):
            fields = line.split("\t")
            ref, alt = fields[3], fields[4]
            assert ref == "G", f"Expected REF=G (A2), got {ref!r}"
            assert alt == "A", f"Expected ALT=A (A1), got {alt!r}"

    def test_vcf_gt_values_are_valid(self, synth_ai2, synth_plink, tmp_path):
        """Every GT call must be one of 0/0, 0/1, 1/1, ./."""
        from alphaimpute2_to_vcf import convert

        out_vcf = str(tmp_path / "out.vcf.gz")
        convert(synth_ai2["genotypes"], synth_plink["bim"], out_vcf)

        valid_gts = {"0/0", "0/1", "1/1", "./."}
        for line in _read_vcf_data_lines(out_vcf):
            fields = line.split("\t")
            for gt in fields[9:]:
                assert gt in valid_gts, f"Invalid GT: {gt!r}"

    def test_missing_encodes_to_dot_slash_dot(self, synth_ai2, synth_plink, tmp_path):
        """Offspring positions that were missing (9) must appear as ./. in VCF."""
        from alphaimpute2_to_vcf import convert

        out_vcf = str(tmp_path / "out.vcf.gz")
        convert(synth_ai2["genotypes"], synth_plink["bim"], out_vcf)

        data_lines = _read_vcf_data_lines(out_vcf)
        geno_matrix = synth_ai2["geno_by_indiv"]  # [sample][snp]
        n_founders = 4

        for snp_idx, line in enumerate(data_lines):
            gts = line.split("\t")[9:]
            for s_idx in range(n_founders, len(gts)):
                if geno_matrix[s_idx][snp_idx] == 9:
                    assert gts[s_idx] == "./.", (
                        f"SNP {snp_idx}, sample {s_idx}: expected './.' for missing"
                    )

    def test_column_count_mismatch_raises(self, synth_plink, tmp_path):
        """convert() must raise ValueError when BIM and genotype file disagree."""
        from alphaimpute2_to_vcf import convert

        # Build a 1-SNP genotype file but the BIM has 80 SNPs → mismatch
        bad_geno = tmp_path / "bad_genotypes.txt"
        with open(bad_geno, "w") as fh:
            fh.write("ind1 0\n")  # only 1 SNP column
            fh.write("ind2 1\n")

        with pytest.raises(ValueError, match="column count mismatch"):
            convert(str(bad_geno), synth_plink["bim"], str(tmp_path / "out.vcf.gz"))


# ── Helpers ───────────────────────────────────────────────────────────────────

def _read_vcf_header(vcf_gz_path: str) -> str:
    """Return the #CHROM header line from a bgzipped VCF."""
    result = subprocess.run(
        ["bgzip", "-cd", vcf_gz_path],
        capture_output=True, text=True, check=True,
    )
    for line in result.stdout.splitlines():
        if line.startswith("#CHROM"):
            return line
    raise AssertionError("No #CHROM header found in VCF")


def _read_vcf_data_lines(vcf_gz_path: str) -> list:
    """Return non-header lines from a bgzipped VCF."""
    result = subprocess.run(
        ["bgzip", "-cd", vcf_gz_path],
        capture_output=True, text=True, check=True,
    )
    return [l for l in result.stdout.splitlines() if l and not l.startswith("#")]
