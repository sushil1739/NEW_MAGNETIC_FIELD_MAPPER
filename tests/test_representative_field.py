
from pathlib import Path

import numpy as np

from representative_field import load_bmag_matrix, profile_statistics, summarize_region


def write_mapper(path: Path):
    rows = []
    for x in (0.0, 100.0):
        for y in (0.0, 100.0):
            for z, b in zip((-1000.0, -500.0, 0.0, 500.0), (1.0, 2.0, 0.0, 3.0)):
                rows.append(f"{x} {y} {z} 0 {b} 0 {abs(b)}")
    path.write_text(
        "# shape 2 2 4\n"
        "0 0 -1000 100 100 500\n"
        + "\n".join(rows)
        + "\n",
        encoding="utf-8",
    )


def test_profile_reads_z_fastest_map_and_excludes_zero(tmp_path: Path):
    mapper = tmp_path / "Mapper.txt"
    write_mapper(mapper)
    bmag, z_mm, _ = load_bmag_matrix(mapper)
    assert bmag.shape == (4, 4)
    np.testing.assert_allclose(z_mm, [-1000, -500, 0, 500])
    stats = profile_statistics(bmag, 0.01)
    np.testing.assert_allclose(stats["mean"][[0, 1, 3]], [1.0, 2.0, 3.0])
    assert np.isnan(stats["mean"][2])
    np.testing.assert_array_equal(stats["n_active"], [4, 4, 0, 4])


def test_region_summary_is_active_cell_mean(tmp_path: Path):
    mapper = tmp_path / "Mapper.txt"
    write_mapper(mapper)
    bmag, z_mm, _ = load_bmag_matrix(mapper)
    row = summarize_region("test", z_mm, bmag, 0.01, -1000, 0)
    assert row["mean_active_B_T"] == 1.5
    assert row["median_active_B_T"] == 1.5
    assert row["active_fraction"] == 1.0
