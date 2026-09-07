# TMS Magnetic Field Mapper

Build a **full-geometry, uniformly spaced magnetic-field map** for the TMS from the supplied quarter-geometry finite-element field solution, validate the result for `ClarkMcGrew/edep-sim`, and make publication/debugging plots from the generated map.

The repository is designed around one practical workflow:

```text
supplied +X,+Y quarter field maps
            │
            ▼
   parse Plate1 / Plate2 maps
            │
            ▼
 place plates along the TMS beam axis
            │
            ▼
 reflect the quarter solution into all 4 XY quadrants
            │
            ▼
 interpolate onto an arbitrary uniform X/Y/Z grid
            │
            ▼
        Mapper.txt
            │
      ┌─────┼─────────┐
      ▼     ▼         ▼
 validate  plot     edep-sim
```

The raw files in `Field_maps/` are treated as immutable source data. The mapper creates new output files and does not rewrite the supplied field solution.

---

## Table of contents

- [What this repository provides](#what-this-repository-provides)
- [Repository layout](#repository-layout)
- [Original field dataset](#original-field-dataset)
- [TMS coordinate convention](#tms-coordinate-convention)
- [How the mapper works](#how-the-mapper-works)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Building a map with arbitrary grid spacing](#building-a-map-with-arbitrary-grid-spacing)
- [Builder CLI reference](#builder-cli-reference)
- [Output format and edep-sim compatibility](#output-format-and-edep-sim-compatibility)
- [Validation](#validation)
- [Plotting](#plotting)
- [Plot CLI reference](#plot-cli-reference)
- [Reference build](#reference-build)
- [Using the map in edep-sim](#using-the-map-in-edep-sim)
- [Tests](#tests)
- [Important assumptions and physics checks](#important-assumptions-and-physics-checks)
- [Troubleshooting](#troubleshooting)
- [Recommended workflow](#recommended-workflow)

---

# What this repository provides

The repository currently provides:

- parsing of the supplied `.fld` magnetic-field files;
- pairing of the two transverse field-map pieces (`Plate1` and `Plate2`);
- placement of all plate solutions along the TMS beam direction using `z_coord.csv`;
- reconstruction of the full detector from a supplied `+X,+Y` quarter geometry;
- correct axial-vector sign handling for magnetic-field reflection;
- arbitrary user-defined **uniform grid spacing** in X, Y, and Z;
- inverse-distance-weighted or nearest-neighbor interpolation;
- protection against uncontrolled extrapolation outside the supplied field mesh;
- memory-conscious generation of large field maps;
- plain-text output compatible with `edep-sim` `ArbBField`;
- validation of grid spacing, ordering, completeness, and field magnitude;
- 2D heatmaps;
- field-colored regular-grid scatter plots;
- quiver/vector plots;
- 3D field samples;
- quarter/full/stored geometry plotting;
- display in either millimetres or metres;
- reusable plotting controls such as color scale, stride, cropping, titles, and camera angle.

The main command-line programs are:

```text
build_mapper.py       build a uniform TMS magnetic-field map
validate_mapper.py    validate the generated ArbBField text file
plot_mapper.py        make 2D, scatter, quiver, and 3D plots
scatter_xy.py         simple coordinate-only XY helper retained for convenience
```

---

# Repository layout

```text
NEW_MAGNETIC_FIELD_MAPPER/
│
├── Field_maps/
│   ├── Plate1_15.fld
│   ├── Plate1_15_1.fld
│   ├── ...
│   ├── Plate2_80_23.fld
│   └── z_coord.csv
│
├── tms_mapper/
│   ├── __init__.py
│   ├── core.py
│   └── edep.py
│
├── tests/
│   ├── test_coordinate_convention.py
│   ├── test_edep.py
│   ├── test_fld_parser.py
│   ├── test_plot_cli.py
│   └── test_reflection.py
│
├── build_mapper.py
├── plot_mapper.py
├── scatter_xy.py
├── validate_mapper.py
├── requirements.txt
├── README.md
└── .gitignore
```

Generated files such as `Mapper.txt`, `Mapper_*.txt`, `plots/`, and `output/` are intentionally ignored by Git because they can become large.

---

# Original field dataset

## What the supplied dataset represents

The repository contains a finite-element magnetic-field solution for **one quarter of the transverse TMS geometry**.

The supplied source quadrant is:

```text
X >= 0
Y >= 0
```

The full detector is reconstructed using the two symmetry planes:

```text
                  +Y
                   ^
                   |
        reflected  |  supplied quarter
                   |
   ----------------+-----------------> +X
                   |
        reflected  |  reflected
                   |
```

Only the transverse geometry is quartered. The detector still has many different longitudinal plate positions along Z.

## Plate1 and Plate2

At each longitudinal plate position, the quarter cross-section is split into two field-map pieces:

- `Plate1`: the part closer to the transverse symmetry line;
- `Plate2`: the outer part of the quarter geometry.

The mapper loads both pieces at the same Z position and treats them as one physical plate solution before interpolation.

## Dataset inventory

There are **160 `.fld` files** plus one `z_coord.csv`.

| Steel thickness | Plate positions | Plate1 maps | Plate2 maps | Total `.fld` maps |
|---|---:|---:|---:|---:|
| 15 mm | 34 | 34 | 34 | 68 |
| 40 mm | 22 | 22 | 22 | 44 |
| 80 mm | 24 | 24 | 24 | 48 |
| **Total** | **80** | **80** | **80** | **160** |

Thus the field solution contains **80 physical longitudinal plate locations**, each represented by a `Plate1`/`Plate2` pair.

## `z_coord.csv`

`Field_maps/z_coord.csv` contains the longitudinal offsets for every map:

```text
Name,z[mm]
Plate1_15,0
Plate1_15_1,65
...
Plate1_40,2210
...
Plate1_80,4190
...
Plate1_80_23,7295
```

`Plate1` and `Plate2` have matching Z offsets.

The current groups are:

| Plate type | Number of positions | Z-offset range | Pitch |
|---|---:|---:|---:|
| 15 mm | 34 | 0 to 2145 mm | 65 mm |
| 40 mm | 22 | 2210 to 4100 mm | 90 mm |
| 80 mm | 24 | 4190 to 7295 mm | 135 mm |

These values are longitudinal **offsets**. The mapper infers the base source Z coordinate from the first field solution and uses the CSV offsets when required.

## `.fld` format

The `.fld` files contain text headers followed by six-column numerical field samples:

```text
x  y  z  Bx  By  Bz
```

For the supplied dataset:

- source spatial coordinates are interpreted as metres by default;
- the mapper converts them to millimetres internally;
- magnetic-field components are used in tesla;
- text/header lines are ignored by the parser.

A different source length unit can be selected with:

```bash
--source-length-unit mm
```

The raw `.fld` sampling is **not a uniform Cartesian grid**. It is the source finite-element solution that must be resampled before use by `edep-sim`.

---

# TMS coordinate convention

The repository uses one physical coordinate convention everywhere:

```text
+Z : along the neutrino-beam direction
+Y : upward, opposite gravity
+X : transverse; +X is to the left of +Z
```

The supplied field columns are interpreted directly as:

```text
x  y  z  Bx  By  Bz
│  │  │   │   │   │
▼  ▼  ▼   ▼   ▼   ▼
X  Y  Z  Bx  By  Bz
```

There is **no numerical axis permutation** in `Mapper.txt`.

The plotting program changes only the visual orientation where useful. It does not rewrite the physical coordinates.

For 2D plots:

```text
XY : transverse view, looking along the beam
     horizontal = X
     vertical   = Y
     fixed      = Z

XZ : top/longitudinal view
     horizontal = Z (beam)
     vertical   = X
     fixed      = Y

YZ : side/longitudinal view
     horizontal = Z (beam)
     vertical   = Y (up)
     fixed      = X
```

For 3D display:

```text
display horizontal axis = X
display depth/beam axis  = Z
display vertical axis    = Y
```

---

# How the mapper works

## 1. Discover the field-map files

`build_mapper.py` calls the mapper library in `tms_mapper/core.py`.

The code scans `Field_maps/` for names matching the current TMS dataset convention:

```text
Plate1_15.fld
Plate1_15_1.fld
...
Plate2_40_12.fld
...
Plate2_80_23.fld
```

The current implementation recognizes the 15, 40, and 80 mm TMS plate families.

## 2. Read `z_coord.csv`

Each `.fld` stem is matched to its longitudinal Z offset in `z_coord.csv`.

The code then groups:

```text
Plate1_<type>_<index>
Plate2_<type>_<index>
```

into one physical longitudinal plate group.

Incomplete or duplicate Plate1/Plate2 pairs are rejected.

## 3. Parse the source field

For each source file, numerical rows are read as:

```text
x y z Bx By Bz
```

Coordinates are converted to millimetres according to `--source-length-unit`.

## 4. Place each plate along Z

The builder supports three placement modes:

```text
auto
source
zcoord
```

### `--z-placement auto` — default

The code compares representative raw `.fld` Z displacements against the offsets in `z_coord.csv`.

- If the raw source positions already follow the expected longitudinal placement, source Z is retained.
- Otherwise the plate maps are shifted using `z_coord.csv`.

### `--z-placement source`

Always trust the raw `.fld` Z values.

### `--z-placement zcoord`

Always position each plate using `z_coord.csv`.

For the supplied dataset, `auto` is the recommended mode.

## 5. Build the full transverse geometry

The source exists only in the `+X,+Y` quadrant.

For a requested full-detector point:

```text
(X, Y, Z)
```

the mapper queries the supplied quarter field at:

```text
(|X|, |Y|, Z)
```

and then applies the correct magnetic-field reflection signs.

This avoids creating a huge duplicated copy of the unstructured source data.

## 6. Reflect the magnetic field as an axial vector

Magnetic field is an **axial vector**, not an ordinary polar vector.

Define:

```text
sx = +1 if X >= 0, -1 if X < 0
sy = +1 if Y >= 0, -1 if Y < 0
```

For a source-quarter field

```text
Bq = (Bxq, Byq, Bzq)
```

the default full-geometry reflection is:

```text
Bx = sy      * Bxq
By = sx      * Byq
Bz = sx * sy * Bzq
```

Equivalent quadrant table:

| Quadrant | Reflected field |
|---|---|
| `+X,+Y` | `( Bx,  By,  Bz)` |
| `-X,+Y` | `( Bx, -By, -Bz)` |
| `+X,-Y` | `(-Bx,  By, -Bz)` |
| `-X,-Y` | `(-Bx, -By,  Bz)` |

The builder exposes diagnostic alternatives:

```bash
--reflection axial    # default / intended magnetic-field behavior
--reflection polar    # diagnostic comparison only
--reflection none     # diagnostic comparison only
```

## 7. Construct an arbitrary uniform output grid

The user chooses:

```text
GRID_X_SIZE = ΔX in mm
GRID_Y_SIZE = ΔY in mm
GRID_Z_SIZE = ΔZ in mm
```

Example:

```bash
--grid-x-size 100
--grid-y-size 100
--grid-z-size 10
```

The grid coordinates are generated from the requested spacing itself.

The code does **not** choose a point count and then stretch the requested spacing with `linspace`.

Conceptually:

```text
X_i = X_min + i ΔX
Y_j = Y_min + j ΔY
Z_k = Z_min + k ΔZ
```

If the source boundary is not an exact multiple of the requested spacing, the output bounds are expanded outward so the exact requested grid spacing is preserved.

Example:

```text
source X max = 3730 mm
requested ΔX = 100 mm

output X max = 3800 mm
```

The extra grid layer is intentional.

## 8. Interpolate the source field

The source finite-element points are unstructured, so the mapper interpolates them onto the uniform output lattice.

Default:

```bash
--interpolation idw
--k-neighbors 8
--idw-power 2
```

This uses inverse-distance weighting over nearby source points found using `scipy.spatial.cKDTree`.

Nearest-neighbor interpolation is also available:

```bash
--interpolation nearest
```

## 9. Limit extrapolation outside the source mesh

By default:

```bash
--support-radius-mm auto
```

The mapper estimates a conservative support radius from the source XY sampling.

Uniform-grid points too far from the supplied source field are assigned:

```text
B = (0, 0, 0)
```

This prevents uncontrolled extrapolation far outside the finite-element mesh.

Other options:

```bash
--support-radius-mm none
```

to permit extrapolation, or:

```bash
--support-radius-mm 200
```

to use an explicit 200 mm support radius.

## 10. Keep inter-plate regions at zero

The output array is initialized to zero.

Each longitudinal field solution fills only Z grid planes lying inside that plate's source slab.

Therefore air/inter-plate regions remain zero unless represented by source field data.

This is why an XY plot taken in a gap can legitimately show:

```text
|B| = 0
```

across the entire slice.

## 11. Write the edep-sim map

The generated field is stored temporarily in a memory-mapped array so large maps do not require keeping multiple full copies in RAM.

The final text file is written in strict `edep-sim` ordering:

```text
Z changes fastest
Y changes next
X changes slowest
```

---

# Installation

## Clone

```bash
git clone https://github.com/sushil1739/NEW_MAGNETIC_FIELD_MAPPER.git
cd NEW_MAGNETIC_FIELD_MAPPER
```

## Create a Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
```

## Install dependencies

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The repository uses:

```text
numpy
scipy
matplotlib
pytest
```

## Verify the environment

```bash
python -m pytest -q
```

Prefer:

```bash
python -m pytest
```

instead of plain `pytest` if multiple Python installations exist on the machine.

---

# Quick start

## 1. Inspect the planned grid first

```bash
python build_mapper.py \
  --grid-x-size 100 \
  --grid-y-size 100 \
  --grid-z-size 10 \
  --dry-run \
  -v
```

This reports:

- coordinate convention;
- source-map count;
- physical plate-group count;
- selected Z-placement mode;
- inferred base Z;
- requested spacing;
- `NX`, `NY`, `NZ`;
- total number of grid points;
- output bounds;
- approximate output size.

## 2. Generate the full map

```bash
python build_mapper.py \
  --grid-x-size 100 \
  --grid-y-size 100 \
  --grid-z-size 10 \
  --geometry full \
  --output Mapper.txt \
  -v
```

## 3. Validate it

```bash
python validate_mapper.py Mapper.txt
```

## 4. Make a field-colored XY scatter plot

Choose a Z slice inside a steel plate:

```bash
mkdir -p plots

python plot_mapper.py Mapper.txt \
  --kind scatter \
  --plane xy \
  --slice 3.32 \
  --component mag \
  --geometry full \
  --units m \
  -o plots/full_xy_scatter_z3320.png
```

## 5. Make the matching vector plot

```bash
python plot_mapper.py Mapper.txt \
  --kind quiver \
  --plane xy \
  --slice 3.32 \
  --component mag \
  --geometry full \
  --units m \
  --stride 1 \
  -o plots/full_xy_quiver_z3320.png
```

---

# Building a map with arbitrary grid spacing

The primary design goal of this repository is that the **output spacing is user controlled**.

The three required spacing arguments are independent:

```text
--grid-x-size
--grid-y-size
--grid-z-size
```

All three are specified in **millimetres**.

## Example: 100 x 100 x 10 mm

```bash
python build_mapper.py \
  --grid-x-size 100 \
  --grid-y-size 100 \
  --grid-z-size 10 \
  -o Mapper_100_100_10.txt \
  -v
```

## Example: 50 x 50 x 5 mm

```bash
python build_mapper.py \
  --grid-x-size 50 \
  --grid-y-size 50 \
  --grid-z-size 5 \
  -o Mapper_50_50_5.txt \
  -v
```

## Example: intentionally coarse debugging grid

```bash
python build_mapper.py \
  --grid-x-size 200 \
  --grid-y-size 200 \
  --grid-z-size 20 \
  -o Mapper_test.txt \
  -v
```

## Always dry-run very fine grids

A small spacing change can increase the number of points dramatically because:

```text
N_total = NX * NY * NZ
```

For example, halving all three spacings increases the point count by roughly a factor of eight.

Use:

```bash
python build_mapper.py \
  --grid-x-size 25 \
  --grid-y-size 25 \
  --grid-z-size 5 \
  --dry-run \
  -v
```

before generating a large map.

The default safety limit is:

```text
50,000,000 points
```

Override intentionally with:

```bash
--max-points <N>
```

or:

```bash
--force
```

---

# Builder CLI reference

Show all current options with:

```bash
python build_mapper.py --help
```

## Input/output

| Option | Meaning | Default |
|---|---|---|
| `--field-dir DIR` | Directory containing `.fld` files and `z_coord.csv` | `Field_maps` |
| `-o`, `--output FILE` | Output mapper text file | `Mapper.txt` |

## Uniform grid

| Option | Meaning |
|---|---|
| `--grid-x-size MM` | Uniform X spacing |
| `--grid-y-size MM` | Uniform Y spacing |
| `--grid-z-size MM` | Uniform Z spacing |

These are required.

## Geometry

| Option | Values | Meaning |
|---|---|---|
| `--geometry` | `full`, `quarter` | Full reflected TMS or supplied quarter only |
| `--reflection` | `axial`, `polar`, `none` | Vector parity used for full reflection |

Default:

```text
geometry   = full
reflection = axial
```

## Source-coordinate interpretation

| Option | Values | Default |
|---|---|---|
| `--source-length-unit` | `m`, `mm` | `m` |
| `--z-placement` | `auto`, `source`, `zcoord` | `auto` |
| `--z-auto-tolerance MM` | tolerance for auto placement | `2.0` |

## Interpolation

| Option | Meaning | Default |
|---|---|---|
| `--interpolation` | `idw` or `nearest` | `idw` |
| `-k`, `--k-neighbors` | IDW neighbor count | `8` |
| `--idw-power` | IDW exponent | `2` |
| `--support-radius-mm` | `auto`, `none`, or explicit distance | `auto` |

## Custom bounds

The output extent can be constrained or extended with:

```text
--x-min MM
--x-max MM
--y-min MM
--y-max MM
--z-min MM
--z-max MM
```

Build-time bounds are always specified in **millimetres**.

Bounds are snapped outward to preserve the exact requested grid spacing.

## Performance/safety

| Option | Meaning | Default |
|---|---|---|
| `--max-points N` | maximum allowed output points | `50,000,000` |
| `--force` | bypass the point-count safety check | off |
| `--workers N` | KD-tree workers; `-1` uses all cores | `-1` |
| `--chunk-size N` | interpolation query chunk size | `200000` |
| `--dry-run` | inspect only, do not write map | off |
| `-v` | verbose logging | off |

---

# Output format and edep-sim compatibility

The generated file is a plain-text `ArbBField` map.

Comments begin with:

```text
#
```

The first non-comment line contains:

```text
OFFSET_X OFFSET_Y OFFSET_Z GRID_X_SIZE GRID_Y_SIZE GRID_Z_SIZE
```

Example:

```text
-3800 -2500 -4000 100 100 10
```

Every subsequent row contains:

```text
X Y Z Bx By Bz Bmag
```

Units are:

```text
X,Y,Z       = mm
Bx,By,Bz    = tesla
Bmag        = tesla
```

where:

```text
Bmag = sqrt(Bx^2 + By^2 + Bz^2)
```

Required ordering:

```text
Z fastest
Y next
X slowest
```

Conceptually:

```python
for X in grid_x:
    for Y in grid_y:
        for Z in grid_z:
            write(X, Y, Z, Bx, By, Bz, Bmag)
```

This ordering is important for the `edep-sim` arbitrary magnetic-field reader.

---

# Validation

Run:

```bash
python validate_mapper.py Mapper.txt
```

The validator checks:

- a valid six-value grid header;
- positive X/Y/Z spacing;
- every coordinate lies on the declared regular lattice;
- strict `Z-fastest, Y-next, X-slowest` row ordering;
- complete Cartesian coverage;
- `rows = NX * NY * NZ`;
- optional `# shape` metadata consistency;
- `Bmag` consistency with `Bx`, `By`, and `Bz`.

For the current 100 x 100 x 10 mm reference map, a successful result is:

```text
VALID edep-sim ArbBField map
  rows          : 2,902,053
  shape         : (77, 51, 739)
  offset [mm]   : (-3800.0, -2500.0, -4000.0)
  spacing [mm]  : (100.0, 100.0, 10.0)
```

## Explicit uniform-spacing check

The validator already verifies the lattice, but the unique coordinate spacing can also be inspected directly:

```bash
python - <<'PY'
import numpy as np
from tms_mapper.edep import load_mapper_data

_, data = load_mapper_data("Mapper.txt")

for name, col in [("X", 0), ("Y", 1), ("Z", 2)]:
    values = np.unique(data[:, col])
    spacing = np.unique(np.round(np.diff(values), 9))
    print(name, spacing)
PY
```

For the reference map:

```text
X [100.]
Y [100.]
Z [10.]
```

---

# Plotting

All primary plots should be made with:

```bash
python plot_mapper.py ...
```

Show the complete CLI:

```bash
python plot_mapper.py --help
```

The plotter reads the generated `Mapper.txt` directly.

It supports:

```text
--kind 2d
--kind scatter
--kind quiver
--kind 3d
```

## Important: plotting units vs stored units

`Mapper.txt` is always stored in **mm**.

The plotting program can display:

```text
--units mm
```

or:

```text
--units m
```

`--slice` uses the same unit selected by `--units`.

Therefore these are equivalent:

```bash
--slice 3320 --units mm
```

and:

```bash
--slice 3.32 --units m
```

Do **not** use:

```bash
--slice 3320 --units m
```

unless you actually mean 3320 metres.

The plotter selects the nearest available regular-grid plane to the requested slice.

---

## Plot geometry modes

```text
--geometry stored
--geometry quarter
--geometry full
```

### `stored`

Plot exactly what is contained in the input map.

### `quarter`

Filter to:

```text
X >= 0
Y >= 0
```

Useful for comparing the full map against the original source quadrant.

### `full`

Use the full geometry.

If the input file itself contains only a quarter map, the plotting code can reflect it for visualization using the same axial-vector parity.

---

## Field components

All color-based plots support:

```text
--component mag
--component bx
--component by
--component bz
```

Meaning:

```text
mag = |B|
bx  = Bx
by  = By
bz  = Bz
```

For direct comparison of multiple plots, fix the same color range with:

```bash
--vmin 0
--vmax 2.2
```

You can also choose a Matplotlib colormap:

```bash
--cmap viridis
```

---

# 2D heatmaps

A 2D heatmap shows one field component on a regular slice.

## XY magnitude slice

```bash
python plot_mapper.py Mapper.txt \
  --kind 2d \
  --plane xy \
  --slice 3.32 \
  --component mag \
  --geometry full \
  --units m \
  -o plots/full_xy_heatmap_z3320.png
```

## XZ slice of `By`

```bash
python plot_mapper.py Mapper.txt \
  --kind 2d \
  --plane xz \
  --slice 0 \
  --component by \
  --geometry full \
  --units m \
  -o plots/xz_by.png
```

## YZ side view

```bash
python plot_mapper.py Mapper.txt \
  --kind 2d \
  --plane yz \
  --slice 0 \
  --component mag \
  --geometry full \
  --units m \
  -o plots/yz_mag.png
```

---

# Field-colored uniform-grid scatter plots

This mode makes the regular grid itself visible.

Every plotted point is one uniform-grid coordinate and its color represents the selected field quantity.

## Full XY scatter colored by `|B|`

```bash
python plot_mapper.py Mapper.txt \
  --kind scatter \
  --plane xy \
  --slice 3.32 \
  --component mag \
  --geometry full \
  --units m \
  -o plots/full_xy_scatter_z3320.png
```

## Quarter-only scatter

```bash
python plot_mapper.py Mapper.txt \
  --kind scatter \
  --plane xy \
  --slice 3.32 \
  --component mag \
  --geometry quarter \
  --units m \
  -o plots/quarter_xy_scatter_z3320.png
```

## Scatter colored by one component

```bash
python plot_mapper.py Mapper.txt \
  --kind scatter \
  --plane xy \
  --slice 3.32 \
  --component by \
  --geometry full \
  --units m \
  -o plots/full_xy_by_z3320.png
```

## Change marker size

```bash
--scatter-size 8
```

---

# Quiver / vector plots

Quiver plots show the **direction of the magnetic field projected into the chosen plane**.

The arrows use:

```text
XY -> Bx, By
XZ -> Bz horizontally, Bx vertically
YZ -> Bz horizontally, By vertically
```

The XZ and YZ ordering follows the chosen display convention in which Z/beam is horizontal.

Arrow color is controlled independently by:

```bash
--component
```

With:

```bash
--component mag
```

the arrows point according to the in-plane vector, while their color represents total `|B|`.

## Full XY quiver

```bash
python plot_mapper.py Mapper.txt \
  --kind quiver \
  --plane xy \
  --slice 3.32 \
  --component mag \
  --geometry full \
  --units m \
  --stride 1 \
  -o plots/full_xy_quiver_z3320.png
```

## Quarter XY quiver

```bash
python plot_mapper.py Mapper.txt \
  --kind quiver \
  --plane xy \
  --slice 3.32 \
  --component mag \
  --geometry quarter \
  --units m \
  --stride 1 \
  -o plots/quarter_xy_quiver_z3320.png
```

## Reduce arrow density

For a cleaner figure:

```bash
--stride 2
```

or:

```bash
--stride 3
```

## Normalize arrow lengths

If direction is more important than vector magnitude:

```bash
--normalize-arrows
```

Color can still represent the field strength.

## Change arrow width

```bash
--quiver-width 0.0035
```

## Manual Matplotlib quiver scaling

```bash
--quiver-scale <value>
```

If omitted, Matplotlib chooses the scale.

---

# 3D field plots

3D plots show a sampled subset of the map, colored by the selected field quantity.

The physical display is:

```text
X = transverse horizontal axis
Z = beam/depth axis
Y = vertical axis
```

## Full 3D magnitude plot

```bash
python plot_mapper.py Mapper.txt \
  --kind 3d \
  --geometry full \
  --component mag \
  --units m \
  --stride 5 \
  --min-field 0.05 \
  -o plots/full_3d.png
```

## Plot more/fewer points

```bash
--stride 3
```

and:

```bash
--max-plot-points 200000
```

control rendering cost.

## Hide very small field values

```bash
--min-field 0.05
```

is useful for reducing visually uninteresting near-zero regions.

## Change the camera

```bash
--elev 22
--azim -58
```

## Crop to a region

All geometry is still stored in the input file; plotting can select a smaller region:

```bash
python plot_mapper.py Mapper.txt \
  --kind 3d \
  --geometry full \
  --component mag \
  --units m \
  --x-range -4 4 \
  --y-range -2.5 2.5 \
  --z-range 3.20 3.35 \
  --stride 2 \
  -o plots/last_plate_region.png
```

The crop ranges use `--units`.

---

# Common plotting controls

These controls are shared where applicable:

```text
--units mm|m
--component mag|bx|by|bz
--geometry stored|quarter|full
--cmap <matplotlib_colormap>
--vmin <tesla>
--vmax <tesla>
--alpha <0..1>
--title "custom title"
--no-grid
--dpi <value>
-o output.png
```

For comparable scatter/quiver/heatmap figures, explicitly use the same:

```bash
--vmin 0 --vmax 2.2
```

rather than allowing each plot to autoscale independently.

---

# Coordinate-only helper

`scatter_xy.py` remains as a small helper for plotting only the XY grid coordinates at a fixed Z value.

Example:

```bash
python scatter_xy.py Mapper.txt \
  --z 3320 \
  -o plots/grid_xy_z3320.png
```

This helper uses Z in millimetres and does not color points by field strength.

For normal analysis, prefer the unified:

```bash
plot_mapper.py --kind scatter
```

mode.

---

# Understanding zero-field slices

A valid mapper can contain whole XY slices with:

```text
B = 0
```

This does not automatically indicate a failure.

Reasons include:

- the requested Z is between steel plates;
- the requested grid layer lies just outside the source boundary after the bounds were snapped to the uniform lattice;
- the target point lies outside the allowed source support radius.

For the current geometry, choose a Z value known to lie inside a plate when inspecting the steel field.

Example:

```bash
--slice 3.32 --units m
```

selects a downstream steel slice.

---

# Plot CLI reference

Show all options:

```bash
python plot_mapper.py --help
```

## Core plot selection

| Option | Values |
|---|---|
| `--kind` | `2d`, `scatter`, `quiver`, `3d` |
| `--geometry` | `stored`, `quarter`, `full` |
| `--component` | `mag`, `bx`, `by`, `bz` |
| `--plane` | `xy`, `xz`, `yz` |
| `--slice` | orthogonal coordinate in `--units` |
| `--units` | `mm`, `m` |

## Cropping

```text
--x-range MIN MAX
--y-range MIN MAX
--z-range MIN MAX
```

## Colors

```text
--cmap
--vmin
--vmax
--alpha
```

## Scatter

```text
--scatter-size
```

## Quiver

```text
--stride
--normalize-arrows
--quiver-scale
--quiver-width
```

## 3D

```text
--stride
--point-size
--min-field
--max-plot-points
--elev
--azim
```

## Presentation/output

```text
--title
--no-grid
--dpi
-o / --output
```

---

# Reference build

A useful reference configuration for the current source dataset is:

```text
ΔX = 100 mm
ΔY = 100 mm
ΔZ = 10 mm
```

Build command:

```bash
python build_mapper.py \
  --grid-x-size 100 \
  --grid-y-size 100 \
  --grid-z-size 10 \
  --geometry full \
  --output Mapper.txt \
  -v
```

For the current dataset, the resulting uniform grid is:

```text
X bounds = -3800 to +3800 mm
Y bounds = -2500 to +2500 mm
Z bounds = -4000 to +3380 mm

NX = 77
NY = 51
NZ = 739

total points = 2,902,053
```

Validation:

```text
rows        = 2,902,053
shape       = (77, 51, 739)
offset      = (-3800, -2500, -4000) mm
spacing     = (100, 100, 10) mm
```

The source field itself extends slightly less far in some directions; the final boundaries are snapped outward so the requested spacing stays exact.

This reference output is useful as a sanity check when setting up the repository on a new machine.

---

# Using the map in edep-sim

Once validated, point the TMS logical volume to the generated map with the GDML auxiliary field entry:

```xml
<auxiliary
    auxtype="ArbBField"
    auxvalue="/absolute/path/to/Mapper.txt"/>
```

The generated map is already written in the coordinate units, field units, and row ordering required by the mapper reader.

Before a large simulation campaign:

1. run `validate_mapper.py`;
2. inspect XY scatter and quiver plots;
3. inspect one longitudinal XZ or YZ plot;
4. confirm the reflected field directions against the expected magnet convention;
5. run a small edep-sim smoke test.

---

# Tests

Run:

```bash
python -m pytest -q
```

The test suite covers:

- `.fld` parsing;
- source coordinate conversion;
- coordinate convention;
- source-to-TMS axis mapping;
- axial-vector quadrant reflection;
- edep-sim header/order validation;
- uniform-grid assumptions;
- plotting CLI parsing;
- metre/mm slice conversion;
- 3D display-axis convention.

When debugging import problems, verify that the virtual environment is active:

```bash
which python
python -c "import tms_mapper; print(tms_mapper.__file__)"
```

---

# Important assumptions and physics checks

## 1. Reflection signs

The current default uses the mathematical axial-vector transformation for a literal mirror of the supplied `+X,+Y` field solution.

This is implemented explicitly and tested.

Before final production use, the reflected vector signs should still be compared with the TMS magnet-current/symmetry convention expected by the field model.

## 2. Source coordinates

The current source files are treated as already using the TMS physical axes:

```text
x -> X
y -> Y
z -> Z
```

The plotting orientation must never be confused with a numerical axis permutation.

## 3. Source coordinate units

The supplied `.fld` coordinates are interpreted as metres by default.

If a future source dataset is already in millimetres:

```bash
--source-length-unit mm
```

must be used.

## 4. Interpolation is a representation choice

The original finite-element mesh is not identical to the generated regular grid.

The output is a sampled/interpolated representation controlled by:

```text
grid spacing
interpolation method
neighbor count
IDW power
support radius
```

Different settings should be validated before physics comparison.

## 5. Current source naming is TMS-specific

The discovery code currently recognizes the existing:

```text
Plate1/Plate2
15/40/80 mm
```

file families.

The **output grid spacing is arbitrary**, but a completely different source naming/layout convention may require adapting the source-discovery logic.

---

# Troubleshooting

## `ModuleNotFoundError: No module named 'tms_mapper'`

Make sure you are in the repository root and using the virtual environment:

```bash
pwd
which python
source .venv/bin/activate
python -m pytest -q
```

Prefer:

```bash
python -m pytest
```

over plain `pytest` if system Python and virtual-environment Python are different.

## Scatter plot is completely zero

Check the selected slice.

A zero XY slice can be an air/inter-plate region.

Try a Z coordinate inside steel.

Also remember:

```bash
--slice 3.32 --units m
```

is equivalent to:

```bash
--slice 3320 --units mm
```

## Requested slice appears at the end of the detector

The plotter chooses the nearest available grid plane.

A common mistake is mixing units:

```text
--slice 3320 --units m
```

means 3320 metres, not 3320 mm, so the nearest available coordinate will be the detector edge.

## Quiver colors look lighter than scatter colors

Quiver plots contain thin arrows separated by white background, while scatter plots contain solid markers.

Use:

```bash
--stride 2
```

for a less crowded vector plot, and fix identical color scales with:

```bash
--vmin 0 --vmax 2.2
```

when comparing figures.

## Grid is too large

Run a dry run:

```bash
python build_mapper.py ... --dry-run -v
```

Increase the spacing or explicitly increase `--max-points`.

## Need to verify uniform spacing

Run:

```bash
python validate_mapper.py Mapper.txt
```

A valid output prints the declared uniform spacing and observed shape.

---

# Recommended workflow

For a new mapper build:

```text
1. Keep Field_maps/ unchanged.
2. Activate the Python virtual environment.
3. Run the tests.
4. Choose ΔX, ΔY, ΔZ.
5. Run build_mapper.py with --dry-run.
6. Inspect output bounds and point count.
7. Generate Mapper.txt.
8. Validate Mapper.txt.
9. Make an XY scatter plot inside steel.
10. Make a matching XY quiver plot.
11. Compare quarter and full geometry.
12. Inspect XZ/YZ longitudinal behavior.
13. Make a 3D diagnostic plot.
14. Verify field reflection/sign convention.
15. Attach Mapper.txt to the TMS ArbBField volume in GDML.
16. Run a small edep-sim smoke test before large production.
```

A typical sequence is:

```bash
source .venv/bin/activate

python -m pytest -q

python build_mapper.py \
  --grid-x-size 100 \
  --grid-y-size 100 \
  --grid-z-size 10 \
  --dry-run \
  -v

python build_mapper.py \
  --grid-x-size 100 \
  --grid-y-size 100 \
  --grid-z-size 10 \
  --geometry full \
  -o Mapper.txt \
  -v

python validate_mapper.py Mapper.txt

python plot_mapper.py Mapper.txt \
  --kind scatter \
  --plane xy \
  --slice 3.32 \
  --component mag \
  --geometry full \
  --units m \
  --vmin 0 \
  --vmax 2.2 \
  -o plots/full_xy_scatter.png

python plot_mapper.py Mapper.txt \
  --kind quiver \
  --plane xy \
  --slice 3.32 \
  --component mag \
  --geometry full \
  --units m \
  --stride 2 \
  --vmin 0 \
  --vmax 2.2 \
  -o plots/full_xy_quiver.png

python plot_mapper.py Mapper.txt \
  --kind 3d \
  --geometry full \
  --component mag \
  --units m \
  --stride 5 \
  --min-field 0.05 \
  --vmin 0 \
  --vmax 2.2 \
  -o plots/full_3d.png
```

---

# Summary

This repository converts the supplied **quarter-geometry, unstructured finite-element TMS magnetic-field solution** into a **full-detector, arbitrary-resolution, uniformly spaced Cartesian field map** suitable for `edep-sim`.

The essential contract is:

```text
INPUT
  supplied .fld field samples
  z_coord.csv
  user-selected ΔX, ΔY, ΔZ

PROCESS
  parse
  place along Z
  combine Plate1 + Plate2
  reflect +X,+Y -> full XY geometry
  interpolate
  zero unsupported/gap regions
  write regular grid

OUTPUT
  Mapper.txt
  X Y Z Bx By Bz |B|
  coordinates in mm
  field in tesla
  Z-fastest ordering
```

Use `validate_mapper.py` before simulation and `plot_mapper.py` to inspect the generated field before treating a map as production-ready.
