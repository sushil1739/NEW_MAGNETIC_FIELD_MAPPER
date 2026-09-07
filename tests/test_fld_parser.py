from pathlib import Path

import numpy as np

from tms_mapper.core import read_fld


def test_read_fld_ignores_headers_and_converts_metres(tmp_path: Path):
    f = tmp_path / "Plate1_15.fld"
    f.write_text(
        'Vector data "VolumeValue(Volume(Plate1_15), <Bx,By,Bz>)"\n'
        "NumElems 2\n"
        "1.0 2.0 -4.0 0.1 0.2 0.3\n"
        "1.1 2.1 -3.985 0.4 0.5 0.6\n"
    )
    xyz, b = read_fld(f, source_length_unit="m")
    np.testing.assert_allclose(
        xyz,
        [[1000.0, 2000.0, -4000.0], [1100.0, 2100.0, -3985.0]],
    )
    np.testing.assert_allclose(b, [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
