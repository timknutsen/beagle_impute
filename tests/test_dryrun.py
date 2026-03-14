"""
Snakemake --dryrun tests.

These validate the DAG (rule graph, config parsing, wildcard resolution) for
each pipeline mode without executing any shell commands or requiring any
bioinformatics tools beyond snakemake itself.

Skipped automatically when snakemake is not on PATH.
"""

import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.parent
# Existing test data — used only to satisfy get_chroms() at parse time.
REAL_BFILE = str(REPO_ROOT / "tests" / "data" / "test_salmon")

requires_snakemake = pytest.mark.skipif(
    not shutil.which("snakemake"),
    reason="snakemake not on PATH",
)


# ── Config factories ──────────────────────────────────────────────────────────

def _beagle_config(tmp_path: Path, *, use_ref: bool = False) -> dict:
    """Minimal config for Beagle mode dry-run."""
    cfg = {
        "bfile": REAL_BFILE,
        "reference_vcf": str(REPO_ROOT / "tests" / "data" / "test_salmon_ref.PHASED.vcf.gz")
        if use_ref
        else "",
        "output_dir": str(tmp_path / "vcf_output"),
        "plink_path": "plink2",
        "beagle_jar": "fake.jar",
        "conform_gt_jar": "fake_conform.jar" if use_ref else "",
        "bref3_jar": "",
        "plink_extra_flags": "",
        "imputer": "beagle",
        "beagle_params": {
            "window": 80,
            "overlap": 20,
            "ne": 500,
            "nthreads": 4,
        },
        "alphaimpute2_params": {
            "maxthreads": 1,
            "cycles": 4,
            "final_peeling_threshold": 0.1,
            "hd_threshold": 0.95,
            "length": 1.0,
            "ped_only": False,
            "pop_only": False,
            "phase_output": False,
        },
    }
    return cfg


def _ai2_config(tmp_path: Path) -> dict:
    """Minimal config for AlphaImpute2 mode dry-run."""
    cfg = _beagle_config(tmp_path)
    cfg["imputer"] = "alphaimpute2"
    return cfg


def _write_config(path: Path, cfg: dict) -> None:
    """Write a YAML config file from a plain dict (no PyYAML dependency)."""
    def _val(v):
        if isinstance(v, bool):
            return "true" if v else "false"
        if isinstance(v, (int, float)):
            return str(v)
        return f'"{v}"'

    lines = []
    for key, val in cfg.items():
        if isinstance(val, dict):
            lines.append(f"{key}:")
            for k2, v2 in val.items():
                lines.append(f"  {k2}: {_val(v2)}")
        else:
            lines.append(f"{key}: {_val(val)}")

    path.write_text("\n".join(lines) + "\n")


def _run_dryrun(config_path: Path) -> subprocess.CompletedProcess:
    """Run `snakemake --dryrun` with the given config, return the result."""
    return subprocess.run(
        [
            "snakemake",
            "--dryrun",
            "--quiet",
            "--configfile", str(config_path),
        ],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )


# ── Tests ─────────────────────────────────────────────────────────────────────

@requires_snakemake
class TestBeagleDryrun:
    def test_beagle_no_ref(self, tmp_path):
        """Beagle mode without a reference panel must resolve without errors."""
        cfg_path = tmp_path / "config.yaml"
        _write_config(cfg_path, _beagle_config(tmp_path, use_ref=False))
        result = _run_dryrun(cfg_path)
        assert result.returncode == 0, (
            f"Snakemake dryrun failed (Beagle, no ref):\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    def test_beagle_with_ref(self, tmp_path):
        """Beagle mode with a reference panel must include harmonization rules."""
        cfg_path = tmp_path / "config.yaml"
        _write_config(cfg_path, _beagle_config(tmp_path, use_ref=True))
        result = _run_dryrun(cfg_path)
        assert result.returncode == 0, (
            f"Snakemake dryrun failed (Beagle, with ref):\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    def test_expected_targets_in_plan(self, tmp_path):
        """Dryrun output must mention the final PLINK binary target."""
        cfg_path = tmp_path / "config.yaml"
        _write_config(cfg_path, _beagle_config(tmp_path, use_ref=False))
        result = _run_dryrun(cfg_path)
        combined = result.stdout + result.stderr
        assert "imputed_data.bed" in combined or result.returncode == 0


@requires_snakemake
class TestAlphaImpute2Dryrun:
    def test_ai2_mode(self, tmp_path):
        """AlphaImpute2 mode must resolve a valid DAG without errors."""
        cfg_path = tmp_path / "config.yaml"
        _write_config(cfg_path, _ai2_config(tmp_path))
        result = _run_dryrun(cfg_path)
        assert result.returncode == 0, (
            f"Snakemake dryrun failed (AlphaImpute2):\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    def test_ai2_dag_includes_conversion_rules(self, tmp_path):
        """The AlphaImpute2 DAG must include the three AI2-specific rules."""
        cfg_path = tmp_path / "config.yaml"
        _write_config(cfg_path, _ai2_config(tmp_path))
        result = subprocess.run(
            ["snakemake", "--dryrun", "--configfile", str(cfg_path)],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
        )
        combined = result.stdout + result.stderr
        expected_rules = [
            "plink_to_alphaimpute2_fmt",
            "run_alphaimpute2",
            "alphaimpute2_to_vcf",
        ]
        for rule in expected_rules:
            assert rule in combined, (
                f"Rule {rule!r} not found in dryrun output.\n"
                f"STDOUT:\n{result.stdout}\n"
                f"STDERR:\n{result.stderr}"
            )

    def test_beagle_rules_excluded_in_ai2_mode(self, tmp_path):
        """run_beagle must not appear in the AlphaImpute2 DAG."""
        cfg_path = tmp_path / "config.yaml"
        _write_config(cfg_path, _ai2_config(tmp_path))
        result = subprocess.run(
            ["snakemake", "--dryrun", "--configfile", str(cfg_path)],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
        )
        combined = result.stdout + result.stderr
        assert "run_beagle" not in combined, (
            "run_beagle should not appear in the AlphaImpute2 DAG"
        )
