#!/usr/bin/env python3
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from tms_mapper.core import build_uniform_mapper


def optional_float(text: str):
    if text.lower() in {"auto", "none"}:
        return text.lower()
    return float(text)


def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Build a uniformly spaced TMS magnetic-field map from Marco's "
            "+x,+y quarter-geometry .fld files. The output is directly "
            "compatible with ClarkMcGrew/edep-sim ArbBField."
        )
    )
    p.add_argument("--field-dir", default="Field_maps",
                   help="Directory containing *.fld and z_coord.csv")
    p.add_argument("-o", "--output", default="Mapper.txt",
                   help="Output edep-sim mapper text file")

    p.add_argument("--grid-x-size", type=float, required=True, metavar="MM",
                   help="Uniform TMS X spacing in mm (transverse; +X left of +Z)")
    p.add_argument("--grid-y-size", type=float, required=True, metavar="MM",
                   help="Uniform TMS Y spacing in mm (+Y opposite gravity)")
    p.add_argument("--grid-z-size", type=float, required=True, metavar="MM",
                   help="Uniform TMS Z spacing in mm (+Z along beam direction)")

    p.add_argument("--geometry", choices=("full", "quarter"), default="full",
                   help="Build full four-quadrant geometry (default) or source quarter")
    p.add_argument("--source-length-unit", choices=("m", "mm"), default="m",
                   help="Coordinate unit used inside .fld files (default: m)")
    p.add_argument("--z-placement", choices=("auto", "source", "zcoord"), default="auto",
                   help="Use raw .fld z, z_coord.csv placement, or auto-detect")
    p.add_argument("--z-auto-tolerance", type=float, default=2.0, metavar="MM",
                   help="Tolerance used by --z-placement auto")

    p.add_argument("--interpolation", choices=("idw", "nearest"), default="idw")
    p.add_argument("-k", "--k-neighbors", type=int, default=8,
                   help="Neighbors for IDW interpolation (default: 8)")
    p.add_argument("--idw-power", type=float, default=2.0)
    p.add_argument(
        "--support-radius-mm",
        type=optional_float,
        default="auto",
        metavar="AUTO|NONE|MM",
        help=(
            "Zero target points farther than this from source mesh. "
            "'auto' estimates a radius from the quarter XY mesh; "
            "'none' allows extrapolation. Default: auto."
        ),
    )
    p.add_argument(
        "--reflection",
        choices=("axial", "polar", "none"),
        default="axial",
        help=(
            "Field parity used when reflecting +x,+y into other quadrants. "
            "Magnetic field is an axial vector, so 'axial' is the default."
        ),
    )

    for axis in ("x", "y", "z"):
        p.add_argument(f"--{axis}-min", type=float, default=None, metavar="MM")
        p.add_argument(f"--{axis}-max", type=float, default=None, metavar="MM")

    p.add_argument("--max-points", type=int, default=50_000_000,
                   help="Safety limit before generation (default: 50 million)")
    p.add_argument("--force", action="store_true",
                   help="Allow output larger than --max-points")
    p.add_argument("--workers", type=int, default=-1,
                   help="cKDTree query workers; -1 uses all cores")
    p.add_argument("--chunk-size", type=int, default=200_000)
    p.add_argument("--dry-run", action="store_true",
                   help="Inspect geometry/grid size without writing Mapper.txt")
    p.add_argument("-v", "--verbose", action="count", default=0)
    return p


def main() -> int:
    args = make_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose >= 2 else
              logging.INFO if args.verbose >= 1 else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )

    bounds = (
        args.x_min, args.x_max,
        args.y_min, args.y_max,
        args.z_min, args.z_max,
    )
    summary = build_uniform_mapper(
        field_dir=args.field_dir,
        output=args.output,
        grid_x_size_mm=args.grid_x_size,
        grid_y_size_mm=args.grid_y_size,
        grid_z_size_mm=args.grid_z_size,
        geometry=args.geometry,
        source_length_unit=args.source_length_unit,
        z_placement=args.z_placement,
        z_auto_tolerance_mm=args.z_auto_tolerance,
        interpolation=args.interpolation,
        k_neighbors=args.k_neighbors,
        idw_power=args.idw_power,
        support_radius_mm=args.support_radius_mm,
        reflection=args.reflection,
        bounds=bounds,
        max_points=args.max_points,
        force=args.force,
        workers=args.workers,
        chunk_size=args.chunk_size,
        dry_run=args.dry_run,
    )

    print("\nTMS mapper plan" if args.dry_run else "\nTMS mapper written")
    print("  coordinates   : X transverse/left of +Z; Y up (-gravity); Z beam")
    print(f"  geometry      : {summary.geometry}")
    print(f"  source maps   : {summary.source_maps}")
    print(f"  plate groups  : {summary.plate_groups}")
    print(f"  z placement   : {summary.z_placement}")
    print(f"  base z [mm]   : {summary.base_z_mm:g}")
    print(f"  spacing [mm]  : {summary.grid_spacing_mm}")
    print(f"  shape          : {summary.grid_shape}")
    print(f"  total points   : {summary.total_points:,}")
    print(f"  bounds [mm]    : {summary.grid_bounds_mm}")
    if not args.dry_run:
        print(f"  output         : {summary.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
