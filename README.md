# TMS Magnetic Field Mapper — Full Geometry on an Arbitrary Uniform Grid

This repository provides a mapper pipeline for the new TMS magnetic-field geometry supplied as
`Field_maps/*.fld` plus `Field_maps/z_coord.csv`.

## What this code does

The supplied finite-element solutions describe the **+x,+y quarter of TMS**.
The builder:

1. Reads every `Plate1_*` and `Plate2_*` `.fld` file.
2. Pairs Plate1 and Plate2 at each physical z position using `z_coord.csv`.
3. Converts the source coordinates to millimetres (`.fld` is treated as metres by default).
4. Detects whether the `.fld` z values are already globally positioned; otherwise it
   places each slab using `z_coord.csv`.
5. Reflects the +x,+y source quarter into all four transverse quadrants.
6. Resamples onto an arbitrary user-requested **uniform Cartesian grid**.
7. Writes a plain-text `Mapper.txt` in the format consumed by
   `ClarkMcGrew/edep-sim` `ArbBField`.
8. Provides a validator plus a plotting CLI for 2D, 3D, and quiver plots.

The raw files in `Field_maps/` are never modified.

## TMS coordinate convention

The mapper uses the TMS global detector axes everywhere — in the generated
`Mapper.txt`, in magnetic-field components, and in plotting:

```text
+Z : along the neutrino-beam direction
+Y : upward, opposite gravity
+X : transverse; +X is to the left of +Z
```

Marco's `.fld` columns are interpreted directly as
`x y z Bx By Bz -> X Y Z Bx By Bz`; no numerical axis permutation is applied
to the field map.  This is intentional: the source ranges themselves identify
the long detector/beam axis as Z.  The plotting utility only changes the
**display orientation** so Y is vertical and Z is visibly the beam direction.

For 2D plots:

```text
XY : transverse view (looking along Z/beam); X horizontal, Y vertical
XZ : top view; Z horizontal (beam), X transverse
YZ : side view; Z horizontal (beam), Y vertical (up)
```

---

## Install

From the repository:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## First run: inspect without generating

```bash
python build_mapper.py \
  --grid-x-size 100 \
  --grid-y-size 100 \
  --grid-z-size 10 \
  --dry-run -v
```

The three grid sizes are **arbitrary positive spacings in millimetres**.  The
requested spacing is authoritative; bounds are expanded outward to land exactly
on that lattice.

## Build the full TMS

```bash
python build_mapper.py \
  --grid-x-size 100 \
  --grid-y-size 100 \
  --grid-z-size 10 \
  --geometry full \
  --output Mapper.txt \
  -v
```

A higher-resolution example:

```bash
python build_mapper.py \
  --grid-x-size 50 \
  --grid-y-size 50 \
  --grid-z-size 5 \
  --output Mapper_50_50_5.txt \
  -v
```

The builder reports `NX`, `NY`, `NZ`, total points, bounds, and an approximate
text-file size before generation.  It refuses grids larger than 50 million
points unless `--max-points` is changed or `--force` is supplied.

---

## edep-sim output contract

The first non-comment line is

```text
offset_x offset_y offset_z GRID_X_SIZE GRID_Y_SIZE GRID_Z_SIZE
```

All coordinates are in **mm** and magnetic field components are in **tesla**.

Every following data row is

```text
X Y Z Bx By Bz Bmag
```

where `Z` is the beam coordinate and `Y` is positive opposite gravity.

Rows are written in the ordering required by edep-sim:

```text
Z changes fastest
Y changes next
X changes slowest
```

The file can be attached in GDML with

```xml
<auxiliary auxtype="ArbBField" auxvalue="/path/to/Mapper.txt"/>
```

## Reflection from the +x,+y source quarter

Magnetic field is an **axial vector**.  For a source value
`Bq=(Bxq,Byq,Bzq)` evaluated at `(|x|,|y|,z)`, the default full-TMS reflection is

```text
sx = +1 for x >= 0, -1 for x < 0
sy = +1 for y >= 0, -1 for y < 0

Bx = sy      * Bxq
By = sx      * Byq
Bz = sx * sy * Bzq
```

Equivalently:

```text
+x,+y : ( Bx,  By,  Bz)
-x,+y : ( Bx, -By, -Bz)
+x,-y : (-Bx,  By, -Bz)
-x,-y : (-Bx, -By,  Bz)
```

The CLI exposes `--reflection polar` and `--reflection none` only for
diagnostics/comparison; `axial` is the physical default for `B`.

---

## Source interpolation and geometry support

The `.fld` points are unstructured finite-element samples, while edep-sim needs a
regular grid.  The default source-to-grid interpolation is inverse-distance
weighting (IDW) with 8 nearest source points:

```bash
--interpolation idw --k-neighbors 8 --idw-power 2
```

Nearest-neighbor mapping is also available:

```bash
--interpolation nearest
```

By default `--support-radius-mm auto` estimates a support radius from the source
quarter XY mesh and writes zero outside that support.  This reduces accidental
field extrapolation across large holes or regions not represented by the source
volume.  For comparison studies:

```bash
--support-radius-mm none
```

or set an explicit distance:

```bash
--support-radius-mm 200
```

---

## z placement

`z_coord.csv` is always used to identify physical plate positions.  The builder
also checks whether each `.fld` already contains the expected global z
displacement.

Default:

```bash
--z-placement auto
```

Force raw source z:

```bash
--z-placement source
```

Force placement by the CSV offsets:

```bash
--z-placement zcoord
```

The base z is inferred from the first physical plate rather than hard-coded.

---

## Validate Mapper.txt before edep-sim

```bash
python validate_mapper.py Mapper.txt
```

The validator checks:

- six-value edep-sim header;
- positive `dx,dy,dz`;
- every coordinate lies on the declared uniform lattice;
- strict `z-fastest, y-next, x-slowest` row ordering;
- complete Cartesian grid (`rows = NX*NY*NZ`);
- `Bmag = sqrt(Bx^2 + By^2 + Bz^2)`.

---

# Plot CLI

All plots are made directly from `Mapper.txt`.

## 2D field slice

XY slice:

```bash
python plot_mapper.py Mapper.txt \
  --kind 2d \
  --plane xy \
  --slice -2000 \
  --component mag \
  --geometry full \
  -o plots/xy_zm2000.png
```

XZ slice of `By`:

```bash
python plot_mapper.py Mapper.txt \
  --kind 2d \
  --plane xz \
  --slice 0 \
  --component by \
  -o plots/xz_by.png
```

Valid components:

```text
mag  bx  by  bz
```

Valid planes:

```text
xy  xz  yz
```

If the requested slice is not exactly a grid coordinate, the nearest grid plane
is used.

## Quarter-only plot from a full Mapper.txt

```bash
python plot_mapper.py Mapper.txt \
  --kind 2d \
  --plane xy \
  --geometry quarter \
  --component mag
```

`--geometry quarter` filters to `x>=0, y>=0`.

## 3D field plot

```bash
python plot_mapper.py Mapper.txt \
  --kind 3d \
  --geometry full \
  --component mag \
  --stride 5 \
  --min-field 0.05 \
  -o plots/tms_3d.png
```

`--stride` and `--max-plot-points` keep rendering manageable.

## Quiver plot

```bash
python plot_mapper.py Mapper.txt \
  --kind quiver \
  --plane xy \
  --slice -2000 \
  --geometry full \
  --stride 2 \
  -o plots/quiver_xy.png
```

For each plane, the arrows use the corresponding in-plane field components:

```text
xy -> (Bx, By)
xz -> (Bx, Bz)
yz -> (By, Bz)
```

## Plot a full geometry even if Mapper.txt stores only the source quarter

If a mapper file has only non-negative x and y:

```bash
python plot_mapper.py QuarterMapper.txt \
  --kind 3d \
  --geometry full
```

the plotting utility reflects it into four quadrants using the same axial-vector
parity.  This is for visualization; the production mapper should normally be
generated directly with `--geometry full`.

---

## Run tests

```bash
pytest -q
```

The initial tests cover `.fld` parsing, coordinate-unit conversion, axial-vector
quadrant reflection, and edep-sim grid ordering.

---

## Recommended first production sequence

```bash
# 1. Inspect source and output dimensions
python build_mapper.py \
  --grid-x-size 100 \
  --grid-y-size 100 \
  --grid-z-size 10 \
  --dry-run -v

# 2. Generate full mapper
python build_mapper.py \
  --grid-x-size 100 \
  --grid-y-size 100 \
  --grid-z-size 10 \
  --output Mapper.txt -v

# 3. Validate exact edep-sim structure
python validate_mapper.py Mapper.txt

# 4. Inspect field
python plot_mapper.py Mapper.txt \
  --kind 2d --plane xy --component mag --geometry full \
  -o plots/xy.png

python plot_mapper.py Mapper.txt \
  --kind quiver --plane xy --geometry full --stride 2 \
  -o plots/quiver_xy.png

python plot_mapper.py Mapper.txt \
  --kind 3d --geometry full --stride 5 \
  -o plots/3d.png
```

## One physics validation still required

The code uses the mathematically correct axial-vector transformation for a
literal mirror of the +x,+y magnetic solution.  Before declaring the new map
production-ready, compare the reflected signs against the TMS magnet-current /
symmetry convention used to produce Marco's quarter model.  That check is kept
explicit instead of burying sign assumptions in the interpolator.
