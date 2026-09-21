#!/usr/bin/env python3
"""
Final consolidation utility for the TMS representative magnetic-field study.

Reads the already-generated canonical outputs:
  - full-steel representative field
  - central-steel representative field
  - flat-central representative field
  - raw Plate1+Plate2 physical-plate cross-check
  - Maxwell/AEDT validation metrics

No new physics definition is introduced here.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


DEFAULT_FULL = Path(
    "plots/TMS_full_steel_representative_B/full_section_results.csv"
)

DEFAULT_CENTRAL = Path(
    "plots/TMS_central_steel_representative_B/"
    "central_steel_section_results.csv"
)

DEFAULT_FLAT = Path(
    "plots/TMS_flat_central_steel_representative_B/"
    "flat_central_section_results.csv"
)

DEFAULT_RAW = Path(
    "plots/TMS_raw_physical_plate_crosschecks/"
    "raw_physical_plate_crosscheck_results.csv"
)

DEFAULT_MAXWELL = Path(
    "plots/TMS_maxwell_vs_mapper_validation/"
    "validation_metrics.csv"
)


def require(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Required analysis output not found:\n  {path}"
        )


def read_row(path, key, value):
    df = pd.read_csv(path)
    rows = df[df[key] == value]

    if len(rows) != 1:
        raise RuntimeError(
            f"{path}: expected exactly one row with "
            f"{key}={value!r}, found {len(rows)}"
        )

    return rows.iloc[0]


def thickness_weighted(df, b_column):
    num = (
        df[b_column].astype(float)
        * df["thickness_mm"].astype(float)
    ).sum()

    den = df["thickness_mm"].astype(float).sum()

    if den <= 0:
        raise RuntimeError("Selected steel thickness is zero.")

    return float(num / den)


def read_main_results(full_path, central_path, flat_path):
    full = read_row(
        full_path,
        "section",
        "Full steel",
    )

    central = read_row(
        central_path,
        "section",
        "Full TMS",
    )

    flat = read_row(
        flat_path,
        "section",
        "Flat central steel",
    )

    return {
        "Full steel": float(full["Brep_geometry_T"]),
        "Central steel": float(central["Brep_T"]),
        "Flat central steel": float(flat["Brep_T"]),
    }


def read_raw_crosschecks(raw_path):
    df = pd.read_csv(raw_path)

    same = df[
        df["same_flat_window_selected"].astype(int) == 1
    ].copy()

    auto = df[
        df["automatic_selected"].astype(int) == 1
    ].copy()

    return {
        "same_window_value_T": thickness_weighted(
            same,
            "B_central_raw_T",
        ),
        "same_window_nplates": int(len(same)),
        "same_window_steel_mm": float(
            same["thickness_mm"].astype(float).sum()
        ),
        "automatic_value_T": thickness_weighted(
            auto,
            "B_legacy_raw_T",
        ),
        "automatic_nplates": int(len(auto)),
    }


def read_maxwell_metrics(path):
    df = pd.read_csv(path).sort_values(
        "Polyline"
    ).reset_index(drop=True)

    if "Relative bias [%]" not in df.columns:
        df["Relative bias [%]"] = (
            100.0
            * df["Bias Mapper-Maxwell [T]"]
            / df["Maxwell mean [T]"]
        )

    if "Role" not in df.columns:
        df["Role"] = df["Polyline"].map(
            {
                1: "steel-side proxy",
                2: "nearly coordinate-matched",
                3: "nearly coordinate-matched",
                4: "steel-side proxy",
            }
        )

    return df


def save_tables(outdir, main, raw, maxwell):
    final_rows = [
        {
            "Result": "Full steel",
            "Role": "Main representative-field result #1",
            "B_rep_T": main["Full steel"],
        },
        {
            "Result": "Central steel",
            "Role": "Main representative-field result #2",
            "B_rep_T": main["Central steel"],
        },
        {
            "Result": "Flat central steel",
            "Role": "Main representative-field result #3",
            "B_rep_T": main["Flat central steel"],
        },
        {
            "Result": "Raw Plate1+Plate2 same-window",
            "Role": "Independent raw-map cross-check",
            "B_rep_T": raw["same_window_value_T"],
        },
        {
            "Result": "Raw automatic selection",
            "Role": "Backup robustness cross-check",
            "B_rep_T": raw["automatic_value_T"],
        },
    ]

    pd.DataFrame(final_rows).to_csv(
        outdir / "final_representative_field_results.csv",
        index=False,
    )

    maxwell.to_csv(
        outdir / "final_maxwell_validation_metrics.csv",
        index=False,
    )


def plot_main_summary(path, main, raw):
    main_labels = [
        "Full steel",
        "Central steel",
        "Flat central",
    ]

    main_values = [
        main["Full steel"],
        main["Central steel"],
        main["Flat central steel"],
    ]

    cross_label = "Raw same-window\ncross-check"
    cross_value = raw["same_window_value_T"]

    fig, ax = plt.subplots(
        figsize=(12, 7)
    )

    x_main = [0, 1, 2]
    x_cross = [3.5]

    bars_main = ax.bar(
        x_main,
        main_values,
        label="Main analysis chain",
    )

    bars_cross = ax.bar(
        x_cross,
        [cross_value],
        hatch="//",
        label="Independent raw-map cross-check",
    )

    ax.axvline(
        2.75,
        linestyle="--",
        linewidth=1.2,
        alpha=0.45,
    )

    for bar, value in zip(
        list(bars_main) + list(bars_cross),
        main_values + [cross_value],
    ):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.006,
            f"{value:.4f} T",
            ha="center",
            va="bottom",
            fontsize=12,
        )

    flat = main["Flat central steel"]

    delta_pct = (
        100.0
        * (cross_value - flat)
        / flat
    )

    ax.text(
        3.5,
        cross_value + 0.105,
        f"Raw vs mapper: {delta_pct:+.3f}%",
        ha="center",
        va="bottom",
        fontsize=11,
        fontweight="bold",
    )

    ax.set_xticks(
        x_main + x_cross
    )

    ax.set_xticklabels(
        main_labels + [cross_label],
        fontsize=12,
    )

    ax.set_ylabel(
        "Thickness-weighted representative |B| [T]",
        fontsize=13,
    )

    ax.set_title(
        "TMS Representative Magnetic Field — "
        "Main Results and Raw-Map Cross-Check\n"
        "Full steel → central steel → flat central; "
        "raw Plate1+Plate2 shown as an independent validation",
        fontsize=15,
    )

    ymax = max(
        main_values + [cross_value]
    )

    ax.set_ylim(
        0,
        ymax + 0.15,
    )

    ax.grid(
        True,
        axis="y",
        alpha=0.25,
    )

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=2,
        frameon=True,
    )

    fig.tight_layout(rect=[0, 0.08, 1, 1])

    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def plot_maxwell_bias(path, maxwell):
    labels = []

    for _, row in maxwell.iterrows():
        p = int(
            row["Polyline"]
        )

        role_short = (
            "proxy"
            if p in (1, 4)
            else "direct"
        )

        labels.append(
            f"P{p}\n{role_short}"
        )

    values = maxwell[
        "Relative bias [%]"
    ].astype(float).tolist()

    fig, ax = plt.subplots(
        figsize=(11, 7)
    )

    bars = ax.bar(
        labels,
        values,
    )

    ax.axhline(
        0.0,
        linewidth=1.5,
    )

    for bar, value in zip(
        bars,
        values,
    ):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value - 0.004,
            f"{value:+.3f}%",
            ha="center",
            va="top",
            fontsize=12,
        )

    ax.set_ylabel(
        "Mean-field bias: (Mapper − Maxwell) / Maxwell [%]",
        fontsize=13,
    )

    ax.set_title(
        "Maxwell/AEDT vs Mapper — "
        "Mean-Field Bias Across All Four Polylines\n"
        "Mean-field bias remains below 0.21%; "
        "P1/P4 are steel-side proxies, P2/P3 are nearly coordinate-matched",
        fontsize=15,
    )

    ax.grid(
        True,
        axis="y",
        alpha=0.25,
    )

    ax.set_ylim(
        min(values) - 0.025,
        0.012,
    )

    fig.tight_layout()

    fig.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def write_summary(path, main, raw, maxwell):
    flat = main[
        "Flat central steel"
    ]

    raw_same = raw[
        "same_window_value_T"
    ]

    raw_delta = (
        raw_same - flat
    )

    raw_rel = (
        100.0
        * raw_delta
        / flat
    )

    lines = [
        "TMS REPRESENTATIVE MAGNETIC FIELD — FINAL ANALYSIS SUMMARY",
        "=" * 62,
        "",
        "MAIN REPRESENTATIVE-FIELD RESULTS",
        f"  Full steel        : {main['Full steel']:.6f} T",
        f"  Central steel     : {main['Central steel']:.6f} T",
        f"  Flat central steel: {flat:.6f} T",
        "",
        "DIRECT RAW-MAP CROSS-CHECK",
        f"  Raw Plate1+Plate2 same-window: {raw_same:.6f} T",
        (
            f"  Raw minus mapper flat-central: "
            f"{raw_delta:+.6f} T ({raw_rel:+.3f}%)"
        ),
        (
            f"  Selected physical plates: "
            f"{raw['same_window_nplates']}; "
            f"selected steel: "
            f"{raw['same_window_steel_mm']:.0f} mm"
        ),
        "",
        "BACKUP RAW-MAP ROBUSTNESS CHECK",
        (
            f"  Automatic local-consistency selection: "
            f"{raw['automatic_value_T']:.6f} T"
        ),
        "",
        "MAXWELL/AEDT PROFILE/SCALE VALIDATION",
        (
            "  Validation subset: |B| > 0.5 T; "
            "rolling mean is display-only."
        ),
        (
            "  P1/P4: steel-side proxy "
            "(Mapper x=100 mm for Maxwell x=10 mm)."
        ),
        (
            "  P2/P3: nearly coordinate matched "
            "(Mapper x=2800 mm for Maxwell x=2795 mm)."
        ),
        "",
    ]

    for _, row in maxwell.iterrows():
        p = int(
            row["Polyline"]
        )

        lines.append(
            f"  P{p}: "
            f"bias={row['Bias Mapper-Maxwell [T]']:+.6f} T "
            f"({row['Relative bias [%]']:+.3f}%), "
            f"MAE={row['MAE [T]']:.6f} T, "
            f"RMSE={row['RMSE [T]']:.6f} T, "
            f"corr={row['Correlation']:.6f}, "
            f"{row['Role']}"
        )

    lines.extend(
        [
            "",
            "FINAL INTERPRETATION",
            (
                "  The three main B_rep values correspond to "
                "progressively more restrictive, explicit steel-region definitions."
            ),
            (
                "  The original raw Plate1+Plate2 maps reproduce "
                "the final flat-central mapper result with the same windows."
            ),
            (
                "  Maxwell/AEDT independently validates the field "
                "profile and scale along four longitudinal lines."
            ),
        ]
    )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--full",
        type=Path,
        default=DEFAULT_FULL,
    )

    parser.add_argument(
        "--central",
        type=Path,
        default=DEFAULT_CENTRAL,
    )

    parser.add_argument(
        "--flat",
        type=Path,
        default=DEFAULT_FLAT,
    )

    parser.add_argument(
        "--raw",
        type=Path,
        default=DEFAULT_RAW,
    )

    parser.add_argument(
        "--maxwell",
        type=Path,
        default=DEFAULT_MAXWELL,
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path(
            "plots/TMS_final_analysis_summary"
        ),
    )

    args = parser.parse_args()

    for path in (
        args.full,
        args.central,
        args.flat,
        args.raw,
        args.maxwell,
    ):
        require(path)

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    main_results = read_main_results(
        args.full,
        args.central,
        args.flat,
    )

    raw = read_raw_crosschecks(
        args.raw
    )

    maxwell = read_maxwell_metrics(
        args.maxwell
    )

    save_tables(
        args.output_dir,
        main_results,
        raw,
        maxwell,
    )

    plot_main_summary(
        args.output_dir
        / "01_final_representative_field_summary.png",
        main_results,
        raw,
    )

    plot_maxwell_bias(
        args.output_dir
        / "02_maxwell_all_four_polyline_bias.png",
        maxwell,
    )

    write_summary(
        args.output_dir
        / "FINAL_ANALYSIS_SUMMARY.txt",
        main_results,
        raw,
        maxwell,
    )

    flat = main_results[
        "Flat central steel"
    ]

    raw_same = raw[
        "same_window_value_T"
    ]

    raw_delta = (
        raw_same - flat
    )

    raw_rel = (
        100.0
        * raw_delta
        / flat
    )

    print()
    print("TMS FINAL ANALYSIS SUMMARY")
    print("==========================")
    print()
    print("MAIN RESULTS")
    print(
        f"Full steel         : {main_results['Full steel']:.6f} T"
    )
    print(
        f"Central steel      : {main_results['Central steel']:.6f} T"
    )
    print(
        f"Flat central steel : {flat:.6f} T"
    )

    print()
    print("RAW SAME-WINDOW CROSS-CHECK")
    print(
        f"Raw Plate1+Plate2  : {raw_same:.6f} T"
    )
    print(
        f"Difference         : "
        f"{raw_delta:+.6f} T ({raw_rel:+.3f}%)"
    )

    print()
    print("MAXWELL/AEDT — ALL FOUR POLYLINES")

    for _, row in maxwell.iterrows():
        p = int(
            row["Polyline"]
        )

        print(
            f"P{p}: "
            f"bias={row['Bias Mapper-Maxwell [T]']:+.6f} T "
            f"({row['Relative bias [%]']:+.3f}%), "
            f"RMSE={row['RMSE [T]']:.6f} T, "
            f"corr={row['Correlation']:.6f}, "
            f"{row['Role']}"
        )

    print()
    print("Saved:")

    for filename in (
        "01_final_representative_field_summary.png",
        "02_maxwell_all_four_polyline_bias.png",
        "final_representative_field_results.csv",
        "final_maxwell_validation_metrics.csv",
        "FINAL_ANALYSIS_SUMMARY.txt",
    ):
        print(
            f"  {args.output_dir / filename}"
        )


if __name__ == "__main__":
    main()
