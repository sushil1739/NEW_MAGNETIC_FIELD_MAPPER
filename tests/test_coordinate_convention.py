import numpy as np

from tms_mapper.core import TMS_COORDINATE_CONVENTION, source_to_tms_axes
from plot_mapper import _slice_for_plane


def test_source_to_tms_axis_mapping_is_explicit_identity():
    xyz = np.array([[10.0, 20.0, -4000.0], [30.0, 40.0, 3375.0]])
    b = np.array([[0.1, 1.2, -0.2], [-0.3, 1.5, 0.04]])
    got_xyz, got_b = source_to_tms_axes(xyz, b)
    np.testing.assert_allclose(got_xyz, xyz)
    np.testing.assert_allclose(got_b, b)
    assert "Z=beam-direction" in TMS_COORDINATE_CONVENTION
    assert "Y=up-opposite-gravity" in TMS_COORDINATE_CONVENTION


def test_longitudinal_plot_axes_put_beam_horizontal_and_y_vertical():
    rows = []
    for x in (0.0, 100.0):
        for y in (0.0, 100.0):
            for z in (-100.0, 0.0, 100.0):
                rows.append([x, y, z, 0.0, 1.0, 0.0, 1.0])
    data = np.asarray(rows)

    _, _, xz_dims, _ = _slice_for_plane(data, "xz", 0.0)
    _, _, yz_dims, _ = _slice_for_plane(data, "yz", 0.0)

    assert xz_dims == (2, 0)  # horizontal Z (beam), vertical X
    assert yz_dims == (2, 1)  # horizontal Z (beam), vertical Y (up)
