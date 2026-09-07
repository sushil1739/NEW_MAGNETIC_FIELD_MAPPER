#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

from tms_mapper.edep import load_mapper_data, read_mapper_header


_COMPONENT_INDEX = {
    "bx": 3,
    "by": 4,
    "bz": 5,
    "mag": 6,
}


def _component(data: np.ndarray, name: str) -> np.ndarray:
    return data[:, _COMPONENT_INDEX[name]]


def _quarter_filter(data: np.ndarray) -> np.ndarray:
    eps = 1e-9
    return data[(data[:, 0] >= -eps) & (data[:, 1] >= -eps)]


def _expand_quarter_for_plot(data: np.ndarray) -> np.ndarray:
    """Reflect quarter data into four quadrants with axial-vector parity."""
    pieces = []
    for sx, sy in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
        d = data.copy()
        d[:, 0] *= sx
        d[:, 1] *= sy
        # Axial-vector parity:
        d[:, 3] *= sy       # Bx
        d[:, 4] *= sx       # By
        d[:, 5] *= sx * sy  # Bz
        d[:, 6] = np.linalg.norm(d[:, 3:6], axis=1)
        pieces.append(d)
    out = np.concatenate(pieces, axis=0)
    # Remove exact duplicate symmetry-plane points.
    key = np.round(out[:, :3], 9)
    _, idx = np.unique(key, axis=0, return_index=True)
    return out[np.sort(idx)]


def _apply_geometry(data: np.ndarray, geometry: str) -> np.ndarray:
    if geometry == "stored":
        return data
    if geometry == "quarter":
        return _quarter_filter(data)
    if geometry == "full":
        if np.min(data[:, 0]) >= -1e-9 and np.min(data[:, 1]) >= -1e-9:
            return _expand_quarter_for_plot(data)
        return data
    raise ValueError(geometry)


def _nearest_slice_value(data: np.ndarray, axis_col: int, requested: float | None) -> float:
    vals = np.unique(data[:, axis_col])
    if requested is None:
        requested = float(vals[len(vals) // 2])
    return float(vals[np.argmin(np.abs(vals - requested))])


def _slice_for_plane(
    data: np.ndarray,
    plane: str,
    requested: float | None,
) -> Tuple[np.ndarray, float, Tuple[int, int], Tuple[str, str]]:
    # TMS physical display convention:
    #   +Z is the beam direction.
    #   +Y is vertical, opposite gravity.
    #   +X is transverse (left of +Z).
    #
    # For the longitudinal views we deliberately draw Z horizontally so the
    # beam direction reads left-to-right, while Y is vertical in the YZ view.
    if plane == "xy":
        fixed_col, dims = 2, (0, 1)
        labels = ("X [mm] — transverse", "Y [mm] — up / -gravity")
    elif plane == "xz":
        fixed_col, dims = 1, (2, 0)
        labels = ("Z [mm] — beam", "X [mm] — transverse")
    elif plane == "yz":
        fixed_col, dims = 0, (2, 1)
        labels = ("Z [mm] — beam", "Y [mm] — up / -gravity")
    else:
        raise ValueError(plane)

    fixed = _nearest_slice_value(data, fixed_col, requested)
    mask = np.isclose(data[:, fixed_col], fixed, atol=1e-8)
    return data[mask], fixed, dims, labels


def _grid_matrix(slice_data: np.ndarray, a_col: int, b_col: int, values: np.ndarray):
    a = np.unique(slice_data[:, a_col])
    b = np.unique(slice_data[:, b_col])
    matrix = np.full((len(b), len(a)), np.nan, dtype=float)
    amap = {round(v, 9): i for i, v in enumerate(a)}
    bmap = {round(v, 9): i for i, v in enumerate(b)}
    for row, value in zip(slice_data, values):
        matrix[bmap[round(row[b_col], 9)], amap[round(row[a_col], 9)]] = value
    return a, b, matrix


def plot_2d(args, data: np.ndarray):
    s, fixed, dims, labels = _slice_for_plane(data, args.plane, args.slice)
    values = _component(s, args.component)
    a, b, matrix = _grid_matrix(s, dims[0], dims[1], values)

    fig, ax = plt.subplots(figsize=(9, 7))
    mesh = ax.pcolormesh(a, b, matrix, shading="auto")
    fig.colorbar(mesh, ax=ax, label=f"{args.component.upper()} [T]")
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_aspect("equal" if args.plane == "xy" else "auto")
    ax.set_title(
        f"TMS {args.component.upper()} — {args.plane.upper()} slice at "
        f"{ {'xy':'z','xz':'y','yz':'x'}[args.plane] }={fixed:g} mm"
    )
    fig.tight_layout()
    return fig


def plot_quiver(args, data: np.ndarray):
    s, fixed, dims, labels = _slice_for_plane(data, args.plane, args.slice)
    stride = max(1, args.stride)

    # Thin on the regular slice grid rather than by arbitrary row order.
    a_vals = np.unique(s[:, dims[0]])[::stride]
    b_vals = np.unique(s[:, dims[1]])[::stride]
    keep_a = np.isin(s[:, dims[0]], a_vals)
    keep_b = np.isin(s[:, dims[1]], b_vals)
    s = s[keep_a & keep_b]

    # Field columns corresponding to position columns X,Y,Z are Bx,By,Bz.
    # Use the same displayed axis order as the position slice so quiver arrows
    # are never silently transposed.
    field_col_for_position_col = {0: 3, 1: 4, 2: 5}
    u = s[:, field_col_for_position_col[dims[0]]]
    v = s[:, field_col_for_position_col[dims[1]]]
    magnitude = np.hypot(u, v)

    fig, ax = plt.subplots(figsize=(9, 7))
    q = ax.quiver(
        s[:, dims[0]], s[:, dims[1]], u, v, magnitude,
        angles="xy", scale_units="xy", scale=args.quiver_scale,
        pivot="mid",
    )
    fig.colorbar(q, ax=ax, label="In-plane |B| [T]")
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_aspect("equal" if args.plane == "xy" else "auto")
    ax.set_title(
        f"TMS B-field quiver — {args.plane.upper()} slice at "
        f"{ {'xy':'z','xz':'y','yz':'x'}[args.plane] }={fixed:g} mm"
    )
    fig.tight_layout()
    return fig


def plot_3d(args, data: np.ndarray):
    stride = max(1, args.stride)
    d = data[::stride]
    if args.min_field > 0:
        d = d[d[:, 6] >= args.min_field]
    if len(d) > args.max_plot_points:
        idx = np.linspace(0, len(d) - 1, args.max_plot_points, dtype=int)
        d = d[idx]

    values = _component(d, args.component)
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    # Display the detector in its physical TMS orientation.  Matplotlib calls
    # the third display axis "z", but that does not mean it must represent the
    # physical TMS Z coordinate.  We draw physical Y vertically, physical Z
    # along the beam direction, and physical X transverse to it.
    sc = ax.scatter(
        d[:, 2], d[:, 0], d[:, 1],
        c=values,
        s=args.point_size,
        alpha=args.alpha,
    )
    fig.colorbar(sc, ax=ax, shrink=0.7, label=f"{args.component.upper()} [T]")
    ax.set_xlabel("Z [mm] — beam")
    ax.set_ylabel("X [mm] — transverse / left of +Z")
    ax.set_zlabel("Y [mm] — up / -gravity")

    z_span = max(np.ptp(d[:, 2]), 1.0)
    x_span = max(np.ptp(d[:, 0]), 1.0)
    y_span = max(np.ptp(d[:, 1]), 1.0)
    ax.set_box_aspect((z_span, x_span, y_span))
    ax.view_init(elev=22, azim=-58)
    ax.set_title(
        f"TMS 3D magnetic field — {args.component.upper()}\n"
        "Z = beam, Y = up (-gravity), X = transverse"
    )
    fig.tight_layout()
    return fig


def make_parser():
    p = argparse.ArgumentParser(
        description=(
            "CLI plotting utility for TMS Mapper.txt files. Supports 2D slices, "
            "3D scatter, and vector quiver plots of the stored/full/quarter geometry. "
            "Plots use TMS axes: +Z beam, +Y opposite gravity, +X transverse."
        )
    )
    p.add_argument("mapper", help="Mapper.txt generated by build_mapper.py")
    p.add_argument("--kind", choices=("2d", "3d", "quiver"), required=True)
    p.add_argument("--geometry", choices=("stored", "quarter", "full"), default="stored")
    p.add_argument("--component", choices=("mag", "bx", "by", "bz"), default="mag")
    p.add_argument("--plane", choices=("xy", "xz", "yz"), default="xy",
                   help="Slice plane for --kind 2d/quiver")
    p.add_argument("--slice", type=float, default=None, metavar="MM",
                   help="Coordinate of orthogonal slice; nearest grid plane is used")
    p.add_argument("--stride", type=int, default=3,
                   help="Grid thinning for quiver/3D (default: 3)")
    p.add_argument("--min-field", type=float, default=0.0,
                   help="3D: omit points with |B| below this value in tesla")
    p.add_argument("--max-plot-points", type=int, default=200_000,
                   help="3D rendering cap after stride/filter")
    p.add_argument("--point-size", type=float, default=2.0)
    p.add_argument("--alpha", type=float, default=0.6)
    p.add_argument("--quiver-scale", type=float, default=None,
                   help="Matplotlib quiver scale; default lets matplotlib choose")
    p.add_argument("-o", "--output", default=None,
                   help="PNG/PDF/SVG output. If omitted, display interactively.")
    p.add_argument("--dpi", type=int, default=180)
    return p


def main() -> int:
    args = make_parser().parse_args()
    header, data = load_mapper_data(args.mapper)
    if len(data) == 0:
        raise SystemExit("Mapper contains no data rows")

    data = _apply_geometry(data, args.geometry)

    if args.kind == "2d":
        fig = plot_2d(args, data)
    elif args.kind == "quiver":
        fig = plot_quiver(args, data)
    else:
        fig = plot_3d(args, data)

    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
        print(f"Wrote {out}")
    else:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
