import numpy as np

from tms_mapper.core import reflect_axial_from_quarter


def test_axial_reflection_all_quadrants():
    bq = np.array([[1.0, 2.0, 3.0]] * 4)
    x = np.array([1.0, -1.0, 1.0, -1.0])
    y = np.array([1.0, 1.0, -1.0, -1.0])

    got = reflect_axial_from_quarter(bq, x, y)

    expected = np.array([
        [ 1.0,  2.0,  3.0],  # +x +y
        [ 1.0, -2.0, -3.0],  # -x +y
        [-1.0,  2.0, -3.0],  # +x -y
        [-1.0, -2.0,  3.0],  # -x -y
    ])
    np.testing.assert_allclose(got, expected)


def test_symmetry_plane_uses_positive_side_parity():
    bq = np.array([[1.0, 2.0, 3.0]])
    got = reflect_axial_from_quarter(bq, np.array([0.0]), np.array([0.0]))
    np.testing.assert_allclose(got, bq)
