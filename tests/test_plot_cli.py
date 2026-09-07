import numpy as np

from plot_mapper import (
    _slice_request_to_mm,
    _xyz_for_3d,
    make_parser,
)


def test_scatter_mode_and_metre_units_are_supported():
    args = make_parser().parse_args([
        "Mapper.txt",
        "--kind", "scatter",
        "--plane", "xy",
        "--slice", "0",
        "--component", "mag",
        "--geometry", "full",
        "--units", "m",
    ])
    assert args.kind == "scatter"
    assert args.units == "m"
    assert args.plane == "xy"


def test_slice_values_in_metres_are_converted_to_stored_mm():
    assert _slice_request_to_mm(1.25, "m") == 1250.0
    assert _slice_request_to_mm(-4000.0, "mm") == -4000.0
    assert _slice_request_to_mm(None, "m") is None


def test_3d_display_is_x_zbeam_yvertical():
    data = np.array([
        [1000.0, 2000.0, -3000.0, 0.0, 0.0, 0.0, 0.0],
        [2000.0, 2500.0, -1000.0, 0.0, 0.0, 0.0, 0.0],
    ])
    x_plot, z_plot, y_plot = _xyz_for_3d(data, "m")
    np.testing.assert_allclose(x_plot, [1.0, 2.0])
    np.testing.assert_allclose(z_plot, [-3.0, -1.0])
    np.testing.assert_allclose(y_plot, [2.0, 2.5])
