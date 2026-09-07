#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

from tms_mapper.edep import iter_mapper_rows, read_mapper_header


DEFAULT_ACTIVE_THRESHOLD_T = 0.01
DEFAULT_SMOOTH_MM = 150.0
DEFAULT_REGION1_END_MM = -500.0
DEFAULT_REGION2_END_MM = 2000.0
DEFAULT_REGION3_END_MM = 2900.0


def moving_average_ignore_nan(values: np.ndarray, window: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if window <= 1:
        return values.copy()
    valid = np.isfinite(values)
    filled = np.where(valid, values, 0.0)
    kernel = np.ones(int(window), dtype=float)
    numerator = np.convolve(filled, kernel, mode="same")
    denominator = np.convolve(valid.astype(float), kernel, mode="same")
    out = numerator / np.maximum(denominator, 1.0)
    out[denominator == 0] = np.nan
    return out


def load_bmag_matrix(mapper_path: Path) -> Tuple[np.ndarray, np.ndarray, object]:
    """Load Bmag only and reshape as (Nxy, Nz), using z-fastest ordering."""
    header = read_mapper_header(mapper_path)
    if header.shape is None:
        raise ValueError(
            "Representative-field analysis requires '# shape NX NY NZ' metadata. "
            "Regenerate Mapper.txt with the current build_mapper.py."
        )

    nx, ny, nz = header.shape
    expected = nx * ny * nz
    bmag = np.empty(expected, dtype=np.float64)

    count = 0
    for row in iter_mapper_rows(mapper_path):
        if count >= expected:
            raise ValueError(
                f"Mapper contains more rows than # shape predicts ({expected})"
            )
        bmag[count] = float(row[6])
        count += 1

    if count != expected:
        raise ValueError(
            f"Mapper row count {count} does not match # shape product {expected}"
        )

    z_mm = header.offset_z_mm + np.arange(nz, dtype=float) * header.dz_mm
    return bmag.reshape(nx * ny, nz), z_mm, header


def profile_statistics(
    bmag: np.ndarray,
    threshold_t: float,
) -> Dict[str, np.ndarray]:
    active = bmag > threshold_t
    masked = np.where(active, bmag, np.nan)

    n_active = np.sum(active, axis=0)
    n_total = np.full(bmag.shape[1], bmag.shape[0], dtype=int)

    # Entire z slices can legitimately contain no active field.  NumPy emits
    # RuntimeWarning for those all-NaN reductions; the slices are explicitly
    # marked NaN below, so suppress only those expected warnings here.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.nanmean(masked, axis=0)
            median = np.nanmedian(masked, axis=0)
            p16 = np.nanpercentile(masked, 16, axis=0)
            p84 = np.nanpercentile(masked, 84, axis=0)

    empty = n_active == 0
    for arr in (mean, median, p16, p84):
        arr[empty] = np.nan

    return {
        "mean": mean,
        "median": median,
        "p16": p16,
        "p84": p84,
        "n_active": n_active,
        "n_total": n_total,
        "active_fraction": n_active / n_total,
    }


def summarize_region(
    name: str,
    z_mm: np.ndarray,
    bmag: np.ndarray,
    threshold_t: float,
    zlo: float,
    zhi: float,
    *,
    include_upper: bool = False,
) -> Dict[str, float]:
    if include_upper:
        zmask = (z_mm >= zlo) & (z_mm <= zhi)
    else:
        zmask = (z_mm >= zlo) & (z_mm < zhi)

    values = bmag[:, zmask]
    active_values = values[values > threshold_t]
    n_cells = int(values.size)
    n_active = int(active_values.size)

    if n_active == 0:
        return {
            "region": name,
            "z_start_mm": float(zlo),
            "z_end_mm": float(zhi),
            "mean_active_B_T": np.nan,
            "median_active_B_T": np.nan,
            "p16_active_B_T": np.nan,
            "p84_active_B_T": np.nan,
            "active_fraction": 0.0,
            "n_active": 0,
            "n_cells": n_cells,
        }

    return {
        "region": name,
        "z_start_mm": float(zlo),
        "z_end_mm": float(zhi),
        "mean_active_B_T": float(np.mean(active_values)),
        "median_active_B_T": float(np.median(active_values)),
        "p16_active_B_T": float(np.percentile(active_values, 16)),
        "p84_active_B_T": float(np.percentile(active_values, 84)),
        "active_fraction": float(n_active / n_cells),
        "n_active": n_active,
        "n_cells": n_cells,
    }


def write_profile_csv(
    path: Path,
    z_mm: np.ndarray,
    raw: Dict[str, np.ndarray],
    smooth: Dict[str, np.ndarray],
) -> None:
    columns = np.column_stack([
        z_mm,
        raw["mean"],
        smooth["mean"],
        raw["median"],
        smooth["median"],
        raw["p16"],
        smooth["p16"],
        raw["p84"],
        smooth["p84"],
        raw["n_active"],
        raw["n_total"],
        raw["active_fraction"],
        smooth["active_fraction"],
    ])
    np.savetxt(
        path,
        columns,
        delimiter=",",
        comments="",
        header=(
            "z_mm,mean_B_raw_T,mean_B_smooth_T,"
            "median_B_raw_T,median_B_smooth_T,"
            "p16_B_raw_T,p16_B_smooth_T,"
            "p84_B_raw_T,p84_B_smooth_T,"
            "n_active,n_total,active_fraction,active_fraction_smooth"
        ),
    )


def write_region_csv(path: Path, rows: List[Dict[str, float]]) -> None:
    fields = [
        "region",
        "z_start_mm",
        "z_end_mm",
        "mean_active_B_T",
        "median_active_B_T",
        "p16_active_B_T",
        "p84_active_B_T",
        "active_fraction",
        "n_active",
        "n_cells",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot_profile(
    path: Path,
    z_mm: np.ndarray,
    raw: Dict[str, np.ndarray],
    smooth: Dict[str, np.ndarray],
    regions: List[Dict[str, float]],
    threshold_t: float,
) -> None:
    fig, ax = plt.subplots(figsize=(13, 6))
    ax.fill_between(
        z_mm, smooth["p16"], smooth["p84"], alpha=0.18,
        label="16–84% |B| spread across active x-y points",
    )
    ax.plot(z_mm, raw["mean"], linewidth=0.8, alpha=0.25, label="Mean |B| raw")
    ax.plot(z_mm, smooth["mean"], linewidth=2.4, label="Mean |B| smoothed")
    ax.plot(
        z_mm, smooth["median"], linewidth=2.0, linestyle="--",
        label="Median |B| smoothed",
    )

    xmin, xmax = float(z_mm.min()), float(z_mm.max())
    for row in regions[:3]:
        ax.axvspan(row["z_start_mm"], row["z_end_mm"], alpha=0.08)
        x0 = max(0.0, (row["z_start_mm"] - xmin) / (xmax - xmin))
        x1 = min(1.0, (row["z_end_mm"] - xmin) / (xmax - xmin))
        ax.axhline(
            row["mean_active_B_T"], xmin=x0, xmax=x1,
            linewidth=1.4, linestyle=":",
        )
        mid = 0.5 * (row["z_start_mm"] + row["z_end_mm"])
        ax.text(
            mid, row["mean_active_B_T"] + 0.015,
            f'{row["region"]}: {row["mean_active_B_T"]:.3f} T',
            ha="center", va="bottom", fontsize=9,
        )

    for boundary in (
        regions[0]["z_end_mm"],
        regions[1]["z_end_mm"],
        regions[2]["z_end_mm"],
    ):
        ax.axvline(boundary, linestyle=":", linewidth=1.2)

    ax.set_xlabel("TMS local z [mm]")
    ax.set_ylabel("Representative |B| [T]")
    ax.set_title(
        "Representative TMS magnetic field vs local z "
        f"(active points: |B| > {threshold_t:g} T)"
    )
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)


def plot_region_summary(path: Path, rows: List[Dict[str, float]]) -> None:
    labels = [row["region"] for row in rows]
    values = np.asarray([row["mean_active_B_T"] for row in rows], dtype=float)
    p16 = np.asarray([row["p16_active_B_T"] for row in rows], dtype=float)
    p84 = np.asarray([row["p84_active_B_T"] for row in rows], dtype=float)
    lower = np.maximum(0.0, values - p16)
    upper = np.maximum(0.0, p84 - values)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(len(labels))
    bars = ax.bar(x, values)
    ax.errorbar(
        x, values, yerr=np.vstack([lower, upper]),
        fmt="none", capsize=5, linewidth=1.2,
    )
    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.02,
            f"{value:.3f} T",
            ha="center", va="bottom", fontsize=10,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Mean active |B| [T]")
    ax.set_title("Representative magnetic field by working TMS region")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)


def plot_active_fraction(
    path: Path,
    z_mm: np.ndarray,
    raw: Dict[str, np.ndarray],
    smooth: Dict[str, np.ndarray],
    threshold_t: float,
) -> None:
    fig, ax = plt.subplots(figsize=(13, 4.5))
    ax.plot(
        z_mm, raw["active_fraction"], linewidth=0.8,
        alpha=0.25, label="Raw active fraction",
    )
    ax.plot(
        z_mm, smooth["active_fraction"], linewidth=2.0,
        label="Smoothed active fraction",
    )
    ax.set_xlabel("TMS local z [mm]")
    ax.set_ylabel("Fraction of x-y grid points")
    ax.set_ylim(0, 1.02)
    ax.set_title(f"Diagnostic active-field coverage (|B| > {threshold_t:g} T)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=250)
    plt.close(fig)


def make_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Compute expert-facing representative TMS magnetic-field profiles "
            "and region summaries from the new uniform Mapper.txt."
        )
    )
    p.add_argument("mapper", nargs="?", default="Mapper.txt")
    p.add_argument("-o", "--output-dir", default="plots/representative_field")
    p.add_argument("--active-threshold", type=float, default=DEFAULT_ACTIVE_THRESHOLD_T)
    p.add_argument("--smooth-mm", type=float, default=DEFAULT_SMOOTH_MM)
    p.add_argument("--region1-end", type=float, default=DEFAULT_REGION1_END_MM)
    p.add_argument("--region2-end", type=float, default=DEFAULT_REGION2_END_MM)
    p.add_argument("--region3-end", type=float, default=DEFAULT_REGION3_END_MM)
    return p


def main() -> int:
    args = make_parser().parse_args()
    mapper = Path(args.mapper)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    bmag, z_mm, header = load_bmag_matrix(mapper)
    raw = profile_statistics(bmag, args.active_threshold)

    window = max(1, int(round(args.smooth_mm / header.dz_mm)))
    if window % 2 == 0:
        window += 1

    smooth = {
        key: moving_average_ignore_nan(raw[key], window)
        if key in {"mean", "median", "p16", "p84", "active_fraction"}
        else raw[key]
        for key in raw
    }

    r1_start = float(z_mm.min())
    r1_end = float(args.region1_end)
    r2_end = float(args.region2_end)
    r3_end = float(args.region3_end)

    if not (r1_start < r1_end < r2_end < r3_end <= z_mm.max()):
        raise ValueError(
            "Working region boundaries must satisfy "
            "z_min < region1_end < region2_end < region3_end <= z_max"
        )

    regions = [
        summarize_region("Region 1", z_mm, bmag, args.active_threshold, r1_start, r1_end),
        summarize_region("Region 2", z_mm, bmag, args.active_threshold, r1_end, r2_end),
        summarize_region(
            "Region 3", z_mm, bmag, args.active_threshold,
            r2_end, r3_end, include_upper=True,
        ),
        summarize_region(
            "Whole active TMS", z_mm, bmag, args.active_threshold,
            float(z_mm.min()), float(z_mm.max()), include_upper=True,
        ),
    ]

    profile_csv = output_dir / "representative_B_profile_vs_z.csv"
    regions_csv = output_dir / "representative_B_regions.csv"
    profile_png = output_dir / "representative_B_profile_vs_z.png"
    regions_png = output_dir / "representative_B_regions.png"
    active_png = output_dir / "active_field_fraction_vs_z.png"

    write_profile_csv(profile_csv, z_mm, raw, smooth)
    write_region_csv(regions_csv, regions)
    plot_profile(profile_png, z_mm, raw, smooth, regions, args.active_threshold)
    plot_region_summary(regions_png, regions)
    plot_active_fraction(active_png, z_mm, raw, smooth, args.active_threshold)

    print("\nRepresentative magnetic-field summary")
    print(
        f"  grid spacing [mm] : "
        f"({header.dx_mm:g}, {header.dy_mm:g}, {header.dz_mm:g})"
    )
    print(f"  z range [mm]      : {z_mm.min():g} to {z_mm.max():g}")
    print(f"  active threshold  : {args.active_threshold:g} T")
    print(f"  smoothing         : {window} slices ≈ {window * header.dz_mm:g} mm")

    for row in regions:
        print(
            f"  {row['region']:<16} "
            f"mean={row['mean_active_B_T']:.4f} T  "
            f"median={row['median_active_B_T']:.4f} T  "
            f"p16-p84=({row['p16_active_B_T']:.4f}, "
            f"{row['p84_active_B_T']:.4f}) T  "
            f"active={100.0 * row['active_fraction']:.1f}%"
        )

    print("\nSaved:")
    for path in (profile_png, regions_png, active_png, profile_csv, regions_csv):
        print(f"  {path}")

    print(
        "\nNote: these are field-map representative |B| values. "
        "They are not a track-weighted effective field or integral(B dl)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
