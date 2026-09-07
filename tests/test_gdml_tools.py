
from pathlib import Path

import numpy as np
import pytest

from geometry.gdml_tools import (
    compute_field_alignment,
    globalize_mapper,
    inject_arb_bfield,
    inspect_tms_geometry,
    read_mapper_header,
    validate_injected_gdml,
)


FIELD_VOLUMES = (
    "thinvolTMS",
    "thinvol2TMS",
    "thickvolTMS",
    "thickvol2TMS",
    "doublevolTMS",
    "doublevol2TMS",
)


def synthetic_gdml(path: Path, *, pitch80=135.0):
    thin = [-3650.0 + i * 65.0 for i in range(34)]
    thick = [-1422.5 + i * 90.0 for i in range(22)]
    double = [587.5 + i * pitch80 for i in range(24)]

    placements = []
    for i, z in enumerate(thin):
        placements.append(
            f'<physvol><volumeref ref="thinlayervol"/>'
            f'<position name="thin{i}" x="0" y="850" z="{z}" unit="mm"/></physvol>'
        )
    for i, z in enumerate(thick):
        placements.append(
            f'<physvol><volumeref ref="thicklayervol"/>'
            f'<position name="thick{i}" x="0" y="850" z="{z}" unit="mm"/></physvol>'
        )
    for i, z in enumerate(double):
        placements.append(
            f'<physvol><volumeref ref="doublelayervol"/>'
            f'<position name="double{i}" x="0" y="850" z="{z}" unit="mm"/></physvol>'
        )

    steel_volumes = "\n".join(
        f'<volume name="{name}"><materialref ref="Steel"/>'
        f'<solidref ref="dummy"/><auxiliary auxtype="BField" '
        f'auxvalue="(0 T, 1 T, 0 T)"/></volume>'
        for name in FIELD_VOLUMES
    )

    xml = """<?xml version="1.0"?>
<gdml>
 <define/>
 <solids><box name="dummy" x="1" y="1" z="1" lunit="mm"/></solids>
 <structure>
  {steel_volumes}
  <volume name="thinlayervol"/>
  <volume name="thicklayervol"/>
  <volume name="doublelayervol"/>
  <volume name="volTMS">
   {placements}
  </volume>
  <volume name="hall">
   <physvol>
    <volumeref ref="volTMS"/>
    <position name="tms_global" x="18484" y="-3207.5" z="5260" unit="mm"/>
   </physvol>
  </volume>
  <volume name="world">
   <physvol><volumeref ref="hall"/></physvol>
  </volume>
 </structure>
 <setup name="Default" version="1.0"><world ref="world"/></setup>
</gdml>
""".format(steel_volumes=steel_volumes, placements="".join(placements))
    path.write_text(xml, encoding="utf-8")


def mapper_file(path: Path):
    path.write_text(
        "# z_placement=zcoord base_z_mm=-4000\n"
        "# shape 2 1 2\n"
        "-3800 -2500 -4000 100 100 10\n"
        "-3800 -2500 -4000 0 1 0 1\n"
        "-3800 -2500 -3990 0 1 0 1\n"
        "-3700 -2500 -4000 0 1 0 1\n"
        "-3700 -2500 -3990 0 1 0 1\n",
        encoding="utf-8",
    )


def fake_field_dir(path: Path):
    path.mkdir()
    rows = ["Name,z[mm]"]
    for part in (1, 2):
        for thickness, count, first, pitch in (
            (15, 34, 0, 65),
            (40, 22, 2210, 90),
            (80, 24, 4190, 135),
        ):
            for i in range(count):
                suffix = "" if i == 0 else f"_{i}"
                stem = f"Plate{part}_{thickness}{suffix}"
                rows.append(f"{stem},{first + i*pitch}")
                (path / f"{stem}.fld").write_text(
                    "0 0 -4 0 1 0\n0.1 0.1 -3.9 0 1 0\n",
                    encoding="utf-8",
                )
    (path / "z_coord.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")


def test_inspect_tms_global_translation_and_layers(tmp_path: Path):
    gdml = tmp_path / "raw.gdml"
    synthetic_gdml(gdml)

    info = inspect_tms_geometry(gdml)
    assert info.axis_aligned
    np.testing.assert_allclose(info.global_translation_mm, [18484, -3207.5, 5260])
    assert len(info.layer_centers_mm[15.0]) == 34
    assert len(info.layer_centers_mm[40.0]) == 22
    assert len(info.layer_centers_mm[80.0]) == 24


def test_globalize_mapper_changes_offsets_and_xyz_only(tmp_path: Path):
    src = tmp_path / "Mapper.txt"
    dst = tmp_path / "Mapper_global.txt"
    mapper_file(src)

    globalize_mapper(src, dst, [1000, 2000, 3000])
    h = read_mapper_header(dst)
    assert h.offset_mm == (-2800.0, -500.0, -1000.0)
    assert h.spacing_mm == (100.0, 100.0, 10.0)

    numeric = [
        np.fromstring(line, sep=" ")
        for line in dst.read_text().splitlines()
        if line and not line.startswith("#")
    ]
    np.testing.assert_allclose(numeric[1][:6], [-2800, -500, -1000, 0, 1, 0])


def test_inject_replaces_constant_field_on_six_volumes(tmp_path: Path):
    raw = tmp_path / "raw.gdml"
    out = tmp_path / "mapper.gdml"
    synthetic_gdml(raw)

    assert inject_arb_bfield(raw, out, "/pnfs/test/Mapper_global.txt") == 6
    result = validate_injected_gdml(
        out,
        mapper_runtime_path="/pnfs/test/Mapper_global.txt",
    )
    assert result["valid"]
    assert result["count"] == 6


def test_alignment_accepts_translation_compatible_geometry(tmp_path: Path):
    gdml = tmp_path / "raw.gdml"
    mapper = tmp_path / "Mapper.txt"
    fields = tmp_path / "Field_maps"
    synthetic_gdml(gdml, pitch80=135.0)
    mapper_file(mapper)
    fake_field_dir(fields)

    _, report = compute_field_alignment(
        gdml, mapper, fields, tolerance_mm=20.0
    )
    assert report.geometry_compatible
    assert 300 < report.field_to_tms_local_mm[2] < 400
    assert report.field_to_tms_local_mm[1] == pytest.approx(850.0)


def test_alignment_rejects_pitch_mismatch(tmp_path: Path):
    gdml = tmp_path / "raw.gdml"
    mapper = tmp_path / "Mapper.txt"
    fields = tmp_path / "Field_maps"
    synthetic_gdml(gdml, pitch80=130.0)
    mapper_file(mapper)
    fake_field_dir(fields)

    _, report = compute_field_alignment(
        gdml, mapper, fields, tolerance_mm=5.0
    )
    assert not report.geometry_compatible
    assert report.field_pitch_mm[80.0] == pytest.approx(135.0)
    assert report.gdml_pitch_mm[80.0] == pytest.approx(130.0)
    assert report.max_z_residual_mm > 5.0
