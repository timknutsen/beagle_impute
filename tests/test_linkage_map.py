"""scripts/linkage_map_to_plink.py: the Beagle map= file built from a linkage map table."""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from linkage_map_to_plink import build_map, main  # noqa: E402


def _table():
    # Brekke layout: positions in a column called Mb, three cM columns.
    return pd.DataFrame({
        "chr": [10, 2, 2, 2, 2, 1],
        "snpid": ["s6", "s3", "s4", "s5", "s5b", "s1"],
        "Mb": [500, 300, 100, 200, 200, 50],
        "male_cM": [0.0, 1.0, 0.0, 0.5, 0.5, 0.0],
        "female_cM": [0.0, 2.0, 0.0, 1.5, 1.5, 0.0],
        "avg_cM": [0.0, 1.5, 0.0, 1.6, 1.6, 0.0],
    })


def test_sorted_numerically_by_chrom_then_position():
    d, _ = build_map(_table(), "chr", "snpid", "Mb", "avg_cM")
    assert list(d.chrom) == ["1", "2", "2", "2", "10"]
    assert list(d[d.chrom == "2"].bp) == [100, 200, 300]


def test_duplicate_position_dropped_and_cm_dip_raised():
    d, stats = build_map(_table(), "chr", "snpid", "Mb", "avg_cM")
    assert stats["duplicate_positions_dropped"] == 1
    # chr2 avg_cM by position: 0.0, 1.6, 1.5 -> the dip at 300 bp is raised to 1.6
    assert stats["cm_dips_raised"] == 1
    assert list(d[d.chrom == "2"].cm) == [0.0, 1.6, 1.6]


def test_writes_plink_map_columns(tmp_path):
    src = tmp_path / "map.txt"; out = tmp_path / "out.map"
    _table().to_csv(src, sep=" ", index=False)
    assert main([str(src), str(out), "--cm-col", "female_cM"]) == 0
    rows = [line.split("\t") for line in out.read_text().splitlines()]
    assert rows[0] == ["1", "s1", "0.000000", "50"]
    assert all(len(r) == 4 for r in rows)
