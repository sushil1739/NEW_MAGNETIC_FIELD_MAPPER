#!/usr/bin/env python3
"""
TMS_RAW_PHYSICAL_PLATE_CROSSCHECK.py

Canonical CROSS-CHECK analysis for the DUNE TMS representative-field study.

PURPOSE
-------
This script is NOT one of the three main representative-field definitions.

Main analysis chain:
    1. Full steel plate
    2. Central steel region
    3. Flat central steel region

This script performs two independent checks directly from the ORIGINAL
Plate1 + Plate2 .fld maps and z_coord.csv:

A) AUTOMATIC RAW-PHYSICAL-PLATE CHECK
   - preserves the earlier physical-plate method
   - uses the original raw-map transverse selection |y| <= 1780 mm
   - finds the longest locally consistent sequence in R1/R2/R3
   - historically gives a combined value near 1.343 T

B) SAME-WINDOW RAW-PHYSICAL-PLATE CHECK
   - uses the same central transverse bounds as the main analysis when
     applied to the raw-map coordinates:
         10 <= |x| <= 3730 mm
         |y| <= 1780 mm
   - uses the EXACT fixed flat-z windows from the canonical flat-central
     mapper analysis:
         R1: -3300 < z_center < -2300 mm
         R2: -1800 < z_center <   200 mm
         R3:   300 < z_center <  3000 mm
   - combines selected physical plates using the same thickness-weighted
     formula as the main flat-central analysis

The same-window result is therefore the more direct apples-to-apples raw-map
cross-check of the mapper-based flat-central result.

IMPORTANT
---------
z_coord.csv entries are treated as plate START positions, matching the
canonical mapper analyses.  The flat-window selection is made on the physical
plate CENTRE:
    z_center = z_start + thickness/2

No magnetic-field threshold and no smoothing are used for the raw per-plate
means themselves.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Canonical analysis definitions
# ---------------------------------------------------------------------------

BASE_Z_MM = -4000.0

# Main central-steel transverse selection.
CENTRAL_X_MIN_MM = 10.0
CENTRAL_X_MAX_MM = 3730.0
CENTRAL_Y_MAX_MM = 1780.0

# Automatic-cross-check settings from the earlier physical-plate method.
ROLLING_WINDOW = 5
MAX_LOCAL_DEVIATION_T = 0.020
MIN_STABLE_PLATES = 5

# Same fixed flat windows as TMS_FLAT_CENTRAL_STEEL_REPRESENTATIVE_B.py
FLAT_Z_WINDOWS = {
    "R1": (-3300.0, -2300.0),
    "R2": (-1800.0, 200.0),
    "R3": (300.0, 3000.0),
}

REGION_LABELS = {
    "R1": "R1 Thin (15 mm)",
    "R2": "R2 Thick (40 mm)",
    "R3": "R3 Double-thick (80 mm)",
}

REGION_MARKERS = {
    "R1": "o",
    "R2": "s",
    "R3": "^",
}


# ---------------------------------------------------------------------------
# File discovery / parsing
# ---------------------------------------------------------------------------

def locate_field_dir() -> Path:
    if Path("Field_maps").is_dir():
        return Path("Field_maps")

    if Path("z_coord.csv").exists():
        return Path(".")

    raise FileNotFoundError(
        "Could not find Field_maps/ or z_coord.csv. "
        "Run from the project root or the field-map directory."
    )


def read_zcsv(path: Path):
    rows = []

    with path.open() as f:
        for row in csv.reader(f):
            if len(row) < 2:
                continue

            name = row[0].strip()

            try:
                z_start_source_mm = float(row[1])
            except ValueError:
                continue

            rows.append(
                (name, z_start_source_mm)
            )

    return rows


def canonical_key(name: str) -> str:
    s = Path(name).stem
    s = re.sub(
        r"^Plate[12]_",
        "",
        s,
        flags=re.I,
    )
    return s.lower()


def find_fld_files(field_dir: Path):
    files = list(
        field_dir.glob("*.fld")
    )

    if not files:
        files = list(
            field_dir.rglob("*.fld")
        )

    return files


def load_raw_map(path: Path):
    """
    Read numeric raw-field rows:
        x y z Bx By Bz

    Raw source coordinates are assumed to be in metres.
    """
    rows = []

    with path.open(
        errors="ignore"
    ) as f:
        for line in f:
            s = line.strip()

            if (
                not s
                or s.startswith(("#", "%", "//"))
            ):
                continue

            parts = (
                s.replace(",", " ")
                .split()
            )

            if len(parts) < 6:
                continue

            try:
                vals = [
                    float(v)
                    for v in parts[:6]
                ]
            except ValueError:
                continue

            rows.append(vals)

    if not rows:
        raise RuntimeError(
            f"No numeric field rows found in {path}"
        )

    return np.asarray(
        rows,
        dtype=float,
    )


# ---------------------------------------------------------------------------
# Raw-map transverse selections
# ---------------------------------------------------------------------------

def raw_piece_means(path: Path):
    """
    Return two means for one Plate1/Plate2 raw piece:

      legacy_mean:
          |y| <= 1780 mm
          (preserves the previous automatic raw-plate cross-check)

      central_mean:
          10 <= |x| <= 3730 mm
          |y| <= 1780 mm
          (matches the canonical central transverse bounds)

    Also return surviving point counts for correct Plate1+Plate2 combination.
    """
    a = load_raw_map(path)

    x_mm = a[:, 0] * 1000.0
    y_mm = a[:, 1] * 1000.0

    bx = a[:, 3]
    by = a[:, 4]
    bz = a[:, 5]

    bmag = np.sqrt(
        bx * bx
        + by * by
        + bz * bz
    )

    legacy_mask = (
        np.abs(y_mm)
        <= CENTRAL_Y_MAX_MM
    )

    central_mask = (
        (np.abs(x_mm) >= CENTRAL_X_MIN_MM)
        & (np.abs(x_mm) <= CENTRAL_X_MAX_MM)
        & (np.abs(y_mm) <= CENTRAL_Y_MAX_MM)
    )

    if not np.any(legacy_mask):
        raise RuntimeError(
            f"No points survive legacy |y| mask in {path}"
        )

    if not np.any(central_mask):
        raise RuntimeError(
            "\nNo points survive the canonical central x/y mask in:\n"
            f"  {path}\n"
            f"Raw x range = {x_mm.min():.3f} ... {x_mm.max():.3f} mm\n"
            f"Raw y range = {y_mm.min():.3f} ... {y_mm.max():.3f} mm\n"
            "Inspect the raw coordinate convention before proceeding."
        )

    return {
        "legacy_mean": float(
            np.mean(
                bmag[legacy_mask]
            )
        ),
        "legacy_n": int(
            np.sum(
                legacy_mask
            )
        ),
        "central_mean": float(
            np.mean(
                bmag[central_mask]
            )
        ),
        "central_n": int(
            np.sum(
                central_mask
            )
        ),
        "x_min_mm": float(
            np.min(x_mm)
        ),
        "x_max_mm": float(
            np.max(x_mm)
        ),
        "y_min_mm": float(
            np.min(y_mm)
        ),
        "y_max_mm": float(
            np.max(y_mm)
        ),
    }


# ---------------------------------------------------------------------------
# Geometry / selection helpers
# ---------------------------------------------------------------------------

def region_and_thickness(
    z_start_source_mm: float
):
    if z_start_source_mm <= 2145.0:
        return "R1", 15.0

    if z_start_source_mm <= 4100.0:
        return "R2", 40.0

    return "R3", 80.0


def rolling_median(
    values,
    window=ROLLING_WINDOW,
):
    half = window // 2

    return np.array(
        [
            np.median(
                values[
                    max(0, i - half):
                    min(
                        len(values),
                        i + half + 1,
                    )
                ]
            )
            for i in range(
                len(values)
            )
        ],
        dtype=float,
    )


def longest_true_run(mask):
    best = None
    start = None

    for i, good in enumerate(mask):
        if good and start is None:
            start = i

        if (
            start is not None
            and (
                not good
                or i == len(mask) - 1
            )
        ):
            end = (
                i
                if (
                    good
                    and i == len(mask) - 1
                )
                else i - 1
            )

            if (
                best is None
                or end - start
                > best[1] - best[0]
            ):
                best = (
                    start,
                    end,
                )

            start = None

    return best


def thickness_weighted_mean(
    rows,
    field,
):
    if not rows:
        raise RuntimeError(
            "Cannot calculate weighted mean of an empty selection."
        )

    numerator = sum(
        r[field]
        * r["thickness_mm"]
        for r in rows
    )

    denominator = sum(
        r["thickness_mm"]
        for r in rows
    )

    return (
        numerator
        / denominator
    )


def in_fixed_flat_window(row):
    lo, hi = FLAT_Z_WINDOWS[
        row["region"]
    ]

    return (
        lo
        < row["z_center_mapper_mm"]
        < hi
    )


# ---------------------------------------------------------------------------
# Build 80 physical raw-plate means
# ---------------------------------------------------------------------------

def build_physical_plate_rows(
    field_dir: Path,
):
    zcsv = (
        field_dir
        / "z_coord.csv"
    )

    if not zcsv.exists():
        raise FileNotFoundError(
            zcsv
        )

    zrows = read_zcsv(
        zcsv
    )

    flds = find_fld_files(
        field_dir
    )

    by_key = {}

    for p in flds:
        by_key.setdefault(
            canonical_key(
                p.name
            ),
            [],
        ).append(p)

    groups = {}

    for name, z_start_source_mm in zrows:
        key = canonical_key(
            name
        )

        groups.setdefault(
            (
                key,
                z_start_source_mm,
            ),
            [],
        ).append(name)

    physical = []
    missing = []

    first_range_report = None

    for (
        key,
        z_start_source_mm,
    ), _csv_names in sorted(
        groups.items(),
        key=lambda item: item[0][1],
    ):
        candidates = by_key.get(
            key,
            [],
        )

        p1 = next(
            (
                p
                for p in candidates
                if re.search(
                    r"plate1_",
                    p.stem,
                    re.I,
                )
            ),
            None,
        )

        p2 = next(
            (
                p
                for p in candidates
                if re.search(
                    r"plate2_",
                    p.stem,
                    re.I,
                )
            ),
            None,
        )

        if (
            p1 is None
            or p2 is None
        ):
            missing.append(
                (
                    key,
                    z_start_source_mm,
                    p1,
                    p2,
                )
            )
            continue

        m1 = raw_piece_means(
            p1
        )

        m2 = raw_piece_means(
            p2
        )

        if first_range_report is None:
            first_range_report = (
                p1,
                m1,
                p2,
                m2,
            )

        region, thickness_mm = (
            region_and_thickness(
                z_start_source_mm
            )
        )

        # z_coord.csv entry is a physical plate START position.
        z_start_mapper_mm = (
            BASE_Z_MM
            + z_start_source_mm
        )

        z_center_mapper_mm = (
            z_start_mapper_mm
            + 0.5 * thickness_mm
        )

        legacy_combined = (
            m1["legacy_mean"]
            * m1["legacy_n"]
            + m2["legacy_mean"]
            * m2["legacy_n"]
        ) / (
            m1["legacy_n"]
            + m2["legacy_n"]
        )

        central_combined = (
            m1["central_mean"]
            * m1["central_n"]
            + m2["central_mean"]
            * m2["central_n"]
        ) / (
            m1["central_n"]
            + m2["central_n"]
        )

        physical.append(
            {
                "key": key,
                "region": region,
                "thickness_mm": thickness_mm,
                "z_start_source_mm": z_start_source_mm,
                "z_start_mapper_mm": z_start_mapper_mm,
                "z_center_mapper_mm": z_center_mapper_mm,
                "B_legacy_raw_T": legacy_combined,
                "B_central_raw_T": central_combined,
                "legacy_npoints": (
                    m1["legacy_n"]
                    + m2["legacy_n"]
                ),
                "central_npoints": (
                    m1["central_n"]
                    + m2["central_n"]
                ),
            }
        )

    if missing:
        print(
            "\nWARNING: unmatched raw-map physical positions:"
        )

        for (
            key,
            z_start_source_mm,
            p1,
            p2,
        ) in missing[:10]:
            print(
                f"  {key:15s} "
                f"z_start={z_start_source_mm:7.1f}: "
                f"Plate1={p1}, Plate2={p2}"
            )

        if len(missing) > 10:
            print(
                f"  ... plus {len(missing) - 10} more"
            )

    if len(physical) != 80:
        raise RuntimeError(
            f"Matched {len(physical)} physical positions; expected 80."
        )

    return (
        physical,
        first_range_report,
    )


# ---------------------------------------------------------------------------
# Automatic legacy cross-check
# ---------------------------------------------------------------------------

def automatic_crosscheck(
    physical,
):
    selected = {}
    diagnostics = {}

    for region in (
        "R1",
        "R2",
        "R3",
    ):
        rows = [
            r
            for r in physical
            if r["region"] == region
        ]

        values = np.array(
            [
                r["B_legacy_raw_T"]
                for r in rows
            ],
            dtype=float,
        )

        med = rolling_median(
            values
        )

        dev = np.abs(
            values - med
        )

        stable = (
            dev
            <= MAX_LOCAL_DEVIATION_T
        )

        run = longest_true_run(
            stable
        )

        if (
            run is None
            or (
                run[1]
                - run[0]
                + 1
                < MIN_STABLE_PLATES
            )
        ):
            raise RuntimeError(
                f"No >= {MIN_STABLE_PLATES}-plate "
                f"automatic sequence found in {region}"
            )

        start, stop = run

        selected[region] = (
            rows[
                start:
                stop + 1
            ]
        )

        diagnostics[region] = (
            rows,
            med,
            dev,
            stable,
        )

    all_selected = (
        selected["R1"]
        + selected["R2"]
        + selected["R3"]
    )

    region_means = {
        region: thickness_weighted_mean(
            selected[region],
            "B_legacy_raw_T",
        )
        for region in (
            "R1",
            "R2",
            "R3",
        )
    }

    combined = thickness_weighted_mean(
        all_selected,
        "B_legacy_raw_T",
    )

    return (
        selected,
        diagnostics,
        region_means,
        combined,
    )


# ---------------------------------------------------------------------------
# Same-window raw physical-plate cross-check
# ---------------------------------------------------------------------------

def same_window_crosscheck(
    physical,
):
    selected = [
        r
        for r in physical
        if in_fixed_flat_window(
            r
        )
    ]

    selected_by_region = {
        region: [
            r
            for r in selected
            if r["region"] == region
        ]
        for region in (
            "R1",
            "R2",
            "R3",
        )
    }

    region_means = {
        region: thickness_weighted_mean(
            selected_by_region[region],
            "B_central_raw_T",
        )
        for region in (
            "R1",
            "R2",
            "R3",
        )
    }

    combined = thickness_weighted_mean(
        selected,
        "B_central_raw_T",
    )

    return (
        selected,
        selected_by_region,
        region_means,
        combined,
    )


# ---------------------------------------------------------------------------
# Optional mapper-reference loading
# ---------------------------------------------------------------------------

def load_mapper_flat_reference():
    """
    Read the canonical mapper flat-central output if present.

    Returns:
        {"R1": ..., "R2": ..., "R3": ..., "Full": ...}
    or None if unavailable.
    """
    path = Path(
        "plots/"
        "TMS_flat_central_steel_representative_B/"
        "flat_central_section_results.csv"
    )

    if not path.exists():
        return None

    result = {}

    with path.open() as f:
        for row in csv.DictReader(f):
            section = row.get(
                "section",
                "",
            )

            try:
                value = float(
                    row["Brep_T"]
                )
            except (
                KeyError,
                TypeError,
                ValueError,
            ):
                continue

            if section.startswith(
                "R1"
            ):
                result["R1"] = value

            elif section.startswith(
                "R2"
            ):
                result["R2"] = value

            elif section.startswith(
                "R3"
            ):
                result["R3"] = value

            elif section == (
                "Flat central steel"
            ):
                result["Full"] = value

    if set(result) == {
        "R1",
        "R2",
        "R3",
        "Full",
    }:
        return result

    return None


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_all_plate_csv(
    path,
    physical,
    auto_selected,
    same_selected,
):
    auto_keys = {
        r["key"]
        for region in auto_selected.values()
        for r in region
    }

    same_keys = {
        r["key"]
        for r in same_selected
    }

    fields = [
        "key",
        "region",
        "thickness_mm",
        "z_start_source_mm",
        "z_start_mapper_mm",
        "z_center_mapper_mm",
        "B_legacy_raw_T",
        "B_central_raw_T",
        "legacy_npoints",
        "central_npoints",
        "automatic_selected",
        "same_flat_window_selected",
    ]

    with path.open(
        "w",
        newline="",
    ) as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fields,
        )

        writer.writeheader()

        for r in physical:
            out = {
                field: r.get(
                    field,
                    "",
                )
                for field in fields
            }

            out[
                "automatic_selected"
            ] = int(
                r["key"]
                in auto_keys
            )

            out[
                "same_flat_window_selected"
            ] = int(
                r["key"]
                in same_keys
            )

            writer.writerow(
                out
            )


def plot_automatic_crosscheck(
    path,
    physical,
    auto_selected,
    auto_region_means,
    auto_combined,
):
    fig, ax = plt.subplots(
        figsize=(13.5, 7.0)
    )

    ax.plot(
        [
            r["z_center_mapper_mm"]
            for r in physical
        ],
        [
            r["B_legacy_raw_T"]
            for r in physical
        ],
        marker="o",
        markersize=4,
        linewidth=1.0,
        alpha=0.28,
        label="All 80 raw physical plate means",
    )

    for region in (
        "R1",
        "R2",
        "R3",
    ):
        rows = auto_selected[
            region
        ]

        ax.plot(
            [
                r["z_center_mapper_mm"]
                for r in rows
            ],
            [
                r["B_legacy_raw_T"]
                for r in rows
            ],
            marker=REGION_MARKERS[
                region
            ],
            linewidth=1.8,
            label=(
                f"{REGION_LABELS[region]} automatic sequence: "
                f"{auto_region_means[region]:.3f} T"
            ),
        )

        ax.hlines(
            auto_region_means[
                region
            ],
            min(
                r["z_center_mapper_mm"]
                for r in rows
            ),
            max(
                r["z_center_mapper_mm"]
                for r in rows
            ),
            linestyles="--",
            linewidth=1.3,
        )

    ax.axhline(
        auto_combined,
        linestyle="-.",
        linewidth=2.0,
        label=(
            "Automatic raw-plate cross-check "
            f"= {auto_combined:.3f} T"
        ),
    )

    ax.set_xlabel(
        "Physical plate centre z [mm]"
    )

    ax.set_ylabel(
        "Mean |B| from original Plate1 + Plate2 maps [T]"
    )

    ax.set_title(
        "TMS Raw Physical-Plate Cross-Check — Automatic Selection\n"
        "Independent legacy-style local-consistency check"
    )

    ax.grid(
        alpha=0.25
    )

    ax.legend(
        fontsize=8.5
    )

    fig.tight_layout()

    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def plot_same_window_crosscheck(
    path,
    physical,
    selected_by_region,
    region_means,
    combined,
):
    fig, ax = plt.subplots(
        figsize=(13.5, 7.0)
    )

    ax.plot(
        [
            r["z_center_mapper_mm"]
            for r in physical
        ],
        [
            r["B_central_raw_T"]
            for r in physical
        ],
        marker="o",
        markersize=4,
        linewidth=1.0,
        alpha=0.25,
        label=(
            "All 80 raw physical plate means "
            "with canonical central x/y mask"
        ),
    )

    for region in (
        "R1",
        "R2",
        "R3",
    ):
        rows = selected_by_region[
            region
        ]

        lo, hi = FLAT_Z_WINDOWS[
            region
        ]

        ax.axvspan(
            lo,
            hi,
            alpha=0.06,
        )

        ax.plot(
            [
                r["z_center_mapper_mm"]
                for r in rows
            ],
            [
                r["B_central_raw_T"]
                for r in rows
            ],
            marker=REGION_MARKERS[
                region
            ],
            linewidth=1.8,
            label=(
                f"{REGION_LABELS[region]} same flat window: "
                f"{region_means[region]:.3f} T"
            ),
        )

        ax.hlines(
            region_means[
                region
            ],
            lo,
            hi,
            linestyles="--",
            linewidth=1.3,
        )

    ax.axhline(
        combined,
        linestyle="-.",
        linewidth=2.0,
        label=(
            "Same-window raw-plate cross-check "
            f"= {combined:.3f} T"
        ),
    )

    ax.set_xlabel(
        "Physical plate centre z [mm]"
    )

    ax.set_ylabel(
        "Mean |B| from original Plate1 + Plate2 maps [T]"
    )

    ax.set_title(
        "TMS Raw Physical-Plate Cross-Check — Same Flat Windows\n"
        "Same central transverse bounds and fixed longitudinal windows as the main flat-central analysis"
    )

    ax.text(
        0.01,
        0.02,
        (
            "R1 (-3300,-2300), "
            "R2 (-1800,200), "
            "R3 (300,3000) mm"
        ),
        transform=ax.transAxes,
        fontsize=9,
    )

    ax.grid(
        alpha=0.25
    )

    ax.legend(
        fontsize=8.3
    )

    fig.tight_layout()

    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def plot_raw_vs_mapper(
    path,
    raw_region_means,
    raw_combined,
    mapper_ref,
):
    labels = [
        "R1",
        "R2",
        "R3",
        "Combined",
    ]

    raw_vals = [
        raw_region_means[
            "R1"
        ],
        raw_region_means[
            "R2"
        ],
        raw_region_means[
            "R3"
        ],
        raw_combined,
    ]

    mapper_vals = [
        mapper_ref[
            "R1"
        ],
        mapper_ref[
            "R2"
        ],
        mapper_ref[
            "R3"
        ],
        mapper_ref[
            "Full"
        ],
    ]

    x = np.arange(
        len(labels),
        dtype=float,
    )

    width = 0.35

    fig, ax = plt.subplots(
        figsize=(10.5, 6.2)
    )

    bars1 = ax.bar(
        x - width / 2,
        mapper_vals,
        width,
        label="Mapper flat-central main result",
    )

    bars2 = ax.bar(
        x + width / 2,
        raw_vals,
        width,
        label="Raw Plate1+Plate2 same-window cross-check",
    )

    for bars in (
        bars1,
        bars2,
    ):
        for bar in bars:
            value = bar.get_height()

            ax.text(
                bar.get_x()
                + bar.get_width() / 2,
                value,
                f"{value:.3f}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    ax.set_xticks(
        x
    )

    ax.set_xticklabels(
        labels
    )

    ax.set_ylabel(
        "Representative |B| [T]"
    )

    ax.set_title(
        "Flat-Central Main Result vs Raw Physical-Plate Same-Window Cross-Check"
    )

    ax.grid(
        True,
        axis="y",
        alpha=0.25,
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def write_summary(
    path,
    physical,
    auto_selected,
    auto_region_means,
    auto_combined,
    same_selected,
    same_selected_by_region,
    same_region_means,
    same_combined,
    mapper_ref,
):
    lines = [
        "TMS RAW PHYSICAL-PLATE CROSS-CHECKS",
        "=" * 40,
        "",
        "ROLE",
        "  Independent cross-checks only.",
        "  These do not redefine the three main representative-field results.",
        "",
        "MAIN ANALYSIS REFERENCE",
        "  Full steel -> Central steel -> Flat central steel",
        "",
        "CROSS-CHECK A: AUTOMATIC RAW PHYSICAL-PLATE SELECTION",
        "  Raw Plate1+Plate2 maps",
        "  Legacy transverse selection: |y| <= 1780 mm",
        f"  Rolling median window: {ROLLING_WINDOW} plates",
        f"  Local deviation criterion: <= {MAX_LOCAL_DEVIATION_T:.3f} T",
        "",
    ]

    for region in (
        "R1",
        "R2",
        "R3",
    ):
        rows = auto_selected[
            region
        ]

        lines.append(
            f"  {region}: "
            f"n={len(rows):2d}, "
            f"z_center={rows[0]['z_center_mapper_mm']:.1f} to "
            f"{rows[-1]['z_center_mapper_mm']:.1f} mm, "
            f"B={auto_region_means[region]:.6f} T"
        )

    lines.extend(
        [
            "",
            (
                "  Automatic selected-plate thickness-weighted "
                f"cross-check = {auto_combined:.6f} T"
            ),
            "",
            "CROSS-CHECK B: SAME FIXED FLAT WINDOWS",
            (
                "  Raw central transverse selection: "
                "10 <= |x| <= 3730 mm, |y| <= 1780 mm"
            ),
            "  z-window selection is on PHYSICAL PLATE CENTRES.",
            "  R1: -3300 < z_center < -2300 mm",
            "  R2: -1800 < z_center <   200 mm",
            "  R3:   300 < z_center <  3000 mm",
            "",
        ]
    )

    for region in (
        "R1",
        "R2",
        "R3",
    ):
        rows = same_selected_by_region[
            region
        ]

        lines.append(
            f"  {region}: "
            f"n={len(rows):2d}, "
            f"selected steel="
            f"{sum(r['thickness_mm'] for r in rows):.0f} mm, "
            f"B={same_region_means[region]:.6f} T"
        )

    lines.extend(
        [
            "",
            (
                "  Same-window selected plates="
                f"{len(same_selected)}"
            ),
            (
                "  Same-window selected steel="
                f"{sum(r['thickness_mm'] for r in same_selected):.0f} mm"
            ),
            (
                "  Same-window raw physical-plate "
                f"cross-check = {same_combined:.6f} T"
            ),
        ]
    )

    if mapper_ref is not None:
        lines.extend(
            [
                "",
                "MAPPER FLAT-CENTRAL REFERENCE",
                f"  R1={mapper_ref['R1']:.6f} T",
                f"  R2={mapper_ref['R2']:.6f} T",
                f"  R3={mapper_ref['R3']:.6f} T",
                f"  Combined={mapper_ref['Full']:.6f} T",
                "",
                "RAW SAME-WINDOW MINUS MAPPER",
                (
                    f"  R1: "
                    f"{same_region_means['R1'] - mapper_ref['R1']:+.6f} T"
                ),
                (
                    f"  R2: "
                    f"{same_region_means['R2'] - mapper_ref['R2']:+.6f} T"
                ),
                (
                    f"  R3: "
                    f"{same_region_means['R3'] - mapper_ref['R3']:+.6f} T"
                ),
                (
                    f"  Combined: "
                    f"{same_combined - mapper_ref['Full']:+.6f} T"
                ),
            ]
        )

    lines.extend(
        [
            "",
            "IMPORTANT INTERPRETATION",
            "  * Automatic selection is an independent legacy-style raw-map check.",
            "  * Same-window selection is the direct raw-map cross-check of the",
            "    canonical flat-central mapper analysis.",
            "  * Maxwell/AEDT comparison remains a separate field-profile/scale",
            "    validation and is not a detector-average B_rep definition.",
        ]
    )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Canonical raw physical-plate cross-checks for the "
            "TMS representative-field analysis."
        )
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default="plots/TMS_raw_physical_plate_crosschecks",
    )

    args = parser.parse_args()

    outdir = Path(
        args.output_dir
    )

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    field_dir = locate_field_dir()

    print(
        "TMS RAW PHYSICAL-PLATE CROSS-CHECKS"
    )

    print(
        "==================================="
    )

    print()

    print(
        f"Field directory: {field_dir.resolve()}"
    )

    physical, first_range_report = (
        build_physical_plate_rows(
            field_dir
        )
    )

    print(
        f"Matched physical positions: {len(physical)}"
    )

    if first_range_report is not None:
        p1, m1, p2, m2 = first_range_report

        print()

        print(
            "Raw coordinate check from first Plate1/Plate2 pair:"
        )

        print(
            f"  {p1.name}: "
            f"x={m1['x_min_mm']:.1f}..{m1['x_max_mm']:.1f} mm, "
            f"y={m1['y_min_mm']:.1f}..{m1['y_max_mm']:.1f} mm"
        )

        print(
            f"  {p2.name}: "
            f"x={m2['x_min_mm']:.1f}..{m2['x_max_mm']:.1f} mm, "
            f"y={m2['y_min_mm']:.1f}..{m2['y_max_mm']:.1f} mm"
        )

    (
        auto_selected,
        _auto_diagnostics,
        auto_region_means,
        auto_combined,
    ) = automatic_crosscheck(
        physical
    )

    (
        same_selected,
        same_selected_by_region,
        same_region_means,
        same_combined,
    ) = same_window_crosscheck(
        physical
    )

    mapper_ref = (
        load_mapper_flat_reference()
    )

    print()

    print(
        "CROSS-CHECK A — AUTOMATIC RAW-PHYSICAL-PLATE SELECTION"
    )

    print(
        "-------------------------------------------------------"
    )

    for region in (
        "R1",
        "R2",
        "R3",
    ):
        rows = auto_selected[
            region
        ]

        print(
            f"{region}: "
            f"n={len(rows):2d}, "
            f"B={auto_region_means[region]:.6f} T"
        )

    print(
        "Automatic raw-plate thickness-weighted "
        f"cross-check = {auto_combined:.6f} T"
    )

    print()

    print(
        "CROSS-CHECK B — SAME FLAT WINDOWS"
    )

    print(
        "---------------------------------"
    )

    for region in (
        "R1",
        "R2",
        "R3",
    ):
        rows = same_selected_by_region[
            region
        ]

        lo, hi = FLAT_Z_WINDOWS[
            region
        ]

        print(
            f"{region}: "
            f"z_center=({lo:.0f},{hi:.0f}) mm, "
            f"n={len(rows):2d}, "
            f"steel={sum(r['thickness_mm'] for r in rows):.0f} mm, "
            f"B={same_region_means[region]:.6f} T"
        )

    print(
        "Same-window raw-plate thickness-weighted "
        f"cross-check = {same_combined:.6f} T"
    )

    if mapper_ref is not None:
        print()

        print(
            "MAPPER FLAT-CENTRAL REFERENCE"
        )

        print(
            "-----------------------------"
        )

        for region in (
            "R1",
            "R2",
            "R3",
        ):
            print(
                f"{region}: "
                f"mapper={mapper_ref[region]:.6f} T, "
                f"raw={same_region_means[region]:.6f} T, "
                f"delta="
                f"{same_region_means[region] - mapper_ref[region]:+.6f} T"
            )

        print(
            "Combined: "
            f"mapper={mapper_ref['Full']:.6f} T, "
            f"raw={same_combined:.6f} T, "
            f"delta="
            f"{same_combined - mapper_ref['Full']:+.6f} T"
        )

    all_csv = (
        outdir
        / "raw_physical_plate_crosscheck_results.csv"
    )

    write_all_plate_csv(
        all_csv,
        physical,
        auto_selected,
        same_selected,
    )

    plot_automatic_crosscheck(
        outdir
        / "01_raw_physical_automatic_crosscheck.png",
        physical,
        auto_selected,
        auto_region_means,
        auto_combined,
    )

    plot_same_window_crosscheck(
        outdir
        / "02_raw_physical_same_flat_windows.png",
        physical,
        same_selected_by_region,
        same_region_means,
        same_combined,
    )

    if mapper_ref is not None:
        plot_raw_vs_mapper(
            outdir
            / "03_same_window_raw_vs_mapper.png",
            same_region_means,
            same_combined,
            mapper_ref,
        )

    summary = (
        outdir
        / "raw_physical_plate_crosscheck_summary.txt"
    )

    write_summary(
        summary,
        physical,
        auto_selected,
        auto_region_means,
        auto_combined,
        same_selected,
        same_selected_by_region,
        same_region_means,
        same_combined,
        mapper_ref,
    )

    print()

    print(
        "Saved:"
    )

    print(
        f"  {outdir / '01_raw_physical_automatic_crosscheck.png'}"
    )

    print(
        f"  {outdir / '02_raw_physical_same_flat_windows.png'}"
    )

    if mapper_ref is not None:
        print(
            f"  {outdir / '03_same_window_raw_vs_mapper.png'}"
        )

    print(
        f"  {all_csv}"
    )

    print(
        f"  {summary}"
    )


if __name__ == "__main__":
    main()
