
#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from geometry.gdml_tools import (
    DUNENDGGD_COMMIT,
    build_upstream_gdml,
    compute_field_alignment,
    ensure_dunendggd,
    globalize_mapper,
    inject_arb_bfield,
    validate_injected_gdml,
)

TARGETS = (
    "tms_nosand",
    "tms",
    "tms_drift1",
    "prism_nosand",
    "prism",
    "prism_drift1",
)


def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Generate DUNE ND GDML with pinned dunendggd, align the TMS field "
            "map to GDML global coordinates, and inject ArbBField into the six "
            "TMS steel logical volumes."
        )
    )
    p.add_argument("--mapper", default="Mapper.txt")
    p.add_argument("--field-dir", default="Field_maps")
    p.add_argument("--target", choices=TARGETS, default="tms_nosand")
    p.add_argument(
        "--dunendggd-dir",
        default="~/.cache/tms-mapper/dunendggd",
        help="Pinned upstream checkout location",
    )
    p.add_argument(
        "--bootstrap",
        action="store_true",
        help="Clone/fetch pinned dunendggd and pip-install it in the active environment",
    )
    p.add_argument("--raw-gdml", default="output/gdml/raw_tms.gdml")
    p.add_argument("--output-gdml", default="output/gdml/tms_mapper.gdml")
    p.add_argument("--output-mapper", default="output/Mapper_global.txt")
    p.add_argument(
        "--mapper-runtime-path",
        default=None,
        help=(
            "Path stored in the GDML ArbBField auxiliary. Default is the "
            "absolute path to --output-mapper."
        ),
    )
    p.add_argument("--tms-volume", default="volTMS")
    p.add_argument("--geometry-tolerance-mm", type=float, default=5.0)
    p.add_argument(
        "--allow-geometry-mismatch",
        action="store_true",
        help="Proceed only as an explicit diagnostic despite plate-layout mismatch",
    )
    p.add_argument("--source-length-unit", choices=("m", "mm"), default="m")
    p.add_argument("--tms-shift-mm", type=float, default=10000.0)
    p.add_argument("--lar-shift-mm", type=float, default=10000.0)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument(
        "--inspect-only",
        action="store_true",
        help="Build/inspect raw GDML and alignment without globalizing/injecting",
    )
    return p


def main() -> int:
    args = make_parser().parse_args()
    upstream = Path(args.dunendggd_dir).expanduser()

    if args.bootstrap:
        print(f"Bootstrapping dunendggd at pinned commit {DUNENDGGD_COMMIT}")
        ensure_dunendggd(upstream, install=True)
    elif not upstream.exists():
        raise SystemExit(
            f"dunendggd checkout not found at {upstream}\n"
            "Run again with --bootstrap."
        )

    raw_gdml = Path(args.raw_gdml)
    print(f"Building upstream GDML target: {args.target}")
    build_upstream_gdml(
        upstream,
        args.target,
        raw_gdml,
        tms_shift_mm=args.tms_shift_mm,
        lar_shift_mm=args.lar_shift_mm,
    )

    geom, report = compute_field_alignment(
        raw_gdml,
        args.mapper,
        args.field_dir,
        tms_volume=args.tms_volume,
        source_length_unit=args.source_length_unit,
        tolerance_mm=args.geometry_tolerance_mm,
    )

    print("\nGDML / field alignment")
    print(f"  TMS global origin [mm] : {geom.global_translation_mm}")
    print(f"  TMS axis aligned       : {geom.axis_aligned}")
    print(f"  field -> TMS local [mm]: {report.field_to_tms_local_mm}")
    print(f"  field -> global [mm]   : {report.field_to_global_mm}")
    print(f"  max Z residual [mm]    : {report.max_z_residual_mm:.3f}")
    print(f"  RMS Z residual [mm]    : {report.rms_z_residual_mm:.3f}")

    for thickness in (15.0, 40.0, 80.0):
        fc, gc = report.family_counts.get(thickness, (0, 0))
        fp = report.field_pitch_mm.get(thickness)
        gp = report.gdml_pitch_mm.get(thickness)
        print(
            f"  {thickness:g} mm family: count field/gdml={fc}/{gc}, "
            f"pitch field/gdml={fp}/{gp} mm"
        )

    for message in report.messages:
        print(f"  - {message}")

    if not report.geometry_compatible and not args.allow_geometry_mismatch:
        raise SystemExit(
            "\nSTOPPED: field map and generated GDML are not translation-compatible.\n"
            "This protects against applying a field map to a different plate geometry.\n"
            "Resolve the mismatch, or use --allow-geometry-mismatch only for a "
            "deliberate diagnostic build."
        )

    if args.inspect_only:
        print("\nInspection complete; no mapper/GDML injection was written.")
        return 0

    global_mapper = Path(args.output_mapper)
    global_header = globalize_mapper(
        args.mapper,
        global_mapper,
        report.field_to_global_mm,
        overwrite=args.overwrite,
    )

    runtime_path = args.mapper_runtime_path or str(global_mapper.resolve())
    output_gdml = Path(args.output_gdml)

    count = inject_arb_bfield(
        raw_gdml,
        output_gdml,
        runtime_path,
        overwrite=args.overwrite,
    )
    validation = validate_injected_gdml(
        output_gdml,
        mapper_runtime_path=runtime_path,
    )

    print("\nGenerated simulation inputs")
    print(f"  global mapper          : {global_mapper}")
    print(f"  global mapper offset   : {global_header.offset_mm} mm")
    print(f"  mapper spacing         : {global_header.spacing_mm} mm")
    print(f"  mapper runtime path    : {runtime_path}")
    print(f"  injected GDML          : {output_gdml}")
    print(f"  ArbBField volumes      : {count}")
    print(f"  GDML validation        : {validation['valid']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
