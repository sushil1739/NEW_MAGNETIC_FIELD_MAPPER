import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

MAPPER_FILE = "Mapper.txt"
OUTDIR = Path("plots/TMS_central_band_weighted_B_flat_regions")
OUTDIR.mkdir(parents=True, exist_ok=True)

N_R1, N_R2, N_R3 = 34, 22, 24
T_R1, T_R2, T_R3 = 15.0, 40.0, 80.0

R1_ZMIN, R1_ZMAX = -3300.0, -2300.0
R2_ZMIN, R2_ZMAX = -1800.0, 200.0
R3_ZMIN, R3_ZMAX = 300.0, 3000.0

# Read Mapper.txt: skip comments and first numeric grid-metadata line.
rows = []
first_numeric = True
with open(MAPPER_FILE) as f:
    for line in f:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if first_numeric:
            first_numeric = False
            continue
        p = s.split()
        if len(p) < 6:
            continue
        try:
            rows.append([float(v) for v in p])
        except ValueError:
            pass

data = np.asarray(rows, dtype=float)
x, y, z = data[:,0], data[:,1], data[:,2]
Bx, By, Bz = data[:,3], data[:,4], data[:,5]
Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)

central = ((np.abs(x) >= 10) & (np.abs(x) <= 3730) & (np.abs(y) <= 1780))
zc, bc = z[central], Bmag[central]

# Per-z central-band profile.
zp, bp = [], []
for zz in np.unique(zc):
    vals = bc[zc == zz]
    if np.any(vals > 0):
        zp.append(zz)
        bp.append(np.mean(vals))
zp, bp = np.asarray(zp), np.asarray(bp)

# Group contiguous z slices into physical plates.
dz = np.diff(zp)
grid_dz = np.min(dz[dz > 0])
breaks = np.where(dz > 1.5 * grid_dz)[0]
groups, start = [], 0
for idx in breaks:
    groups.append(np.arange(start, idx+1))
    start = idx+1
groups.append(np.arange(start, len(zp)))

plates = []
for i, inds in enumerate(groups):
    plates.append({
        "plate": i+1,
        "z": float(np.mean(zp[inds])),
        "B": float(np.mean(bp[inds]))
    })
plates.sort(key=lambda r: r["z"])

for i, r in enumerate(plates):
    if i < N_R1:
        r["region"], r["t"] = "R1", T_R1
    elif i < N_R1 + N_R2:
        r["region"], r["t"] = "R2", T_R2
    else:
        r["region"], r["t"] = "R3", T_R3

def selected(r):
    if r["region"] == "R1":
        return R1_ZMIN < r["z"] < R1_ZMAX
    if r["region"] == "R2":
        return R2_ZMIN < r["z"] < R2_ZMAX
    return R3_ZMIN < r["z"] < R3_ZMAX

flat = [r for r in plates if selected(r)]
r1 = [r for r in flat if r["region"] == "R1"]
r2 = [r for r in flat if r["region"] == "R2"]
r3 = [r for r in flat if r["region"] == "R3"]

def weighted(rows):
    return sum(r["B"]*r["t"] for r in rows) / sum(r["t"] for r in rows)

B1, B2, B3, BF = weighted(r1), weighted(r2), weighted(r3), weighted(flat)

print("\nCENTRAL-BAND FLAT-REGION RESULTS")
print("--------------------------------")
print(f"R1 Thin         : {len(r1)} plates, B_rep = {B1:.6f} T")
print(f"R2 Thick        : {len(r2)} plates, B_rep = {B2:.6f} T")
print(f"R3 Double-thick : {len(r3)} plates, B_rep = {B3:.6f} T")
print(f"Full TMS flat-region thickness-weighted B_rep = {BF:.6f} T")

# Save CSV.
csv = OUTDIR / "flat_region_plate_results.csv"
with open(csv, "w") as f:
    f.write("plate,region,z_mm,thickness_mm,mean_B_T,selected\n")
    for r in plates:
        f.write(f'{r["plate"]},{r["region"]},{r["z"]:.3f},{r["t"]:.1f},{r["B"]:.8f},{int(selected(r))}\n')

# Plot baseline profile + highlighted flat regions.
fig, ax = plt.subplots(figsize=(13,6))
ax.plot([r["z"] for r in plates], [r["B"] for r in plates],
        "o-", markersize=4, linewidth=1, alpha=0.25, label="All central-band plate means")

for rows_, marker, label in [
    (r1, "o", f"R1 flat: B_rep = {B1:.3f} T"),
    (r2, "s", f"R2 flat: B_rep = {B2:.3f} T"),
    (r3, "^", f"R3 flat: B_rep = {B3:.3f} T")
]:
    ax.plot([r["z"] for r in rows_], [r["B"] for r in rows_],
            marker=marker, linewidth=1.5, label=label)

for b, lo, hi in [(B1,R1_ZMIN,R1_ZMAX),(B2,R2_ZMIN,R2_ZMAX),(B3,R3_ZMIN,R3_ZMAX)]:
    ax.hlines(b, lo, hi, linestyles="--")

ax.axhline(BF, linestyle="-.", linewidth=2,
           label=f"Full TMS flat weighted B_rep = {BF:.3f} T")
ax.set_xlabel("TMS local Z [mm]")
ax.set_ylabel("Per-plate mean |B| [T]")
ax.set_title("Central-band representative |B| using flat longitudinal regions")
ax.grid(alpha=0.25)
ax.legend(fontsize=9)
fig.tight_layout()

png = OUTDIR / "01_central_band_flat_regions_B_vs_z.png"
fig.savefig(png, dpi=300, bbox_inches="tight")
plt.close(fig)

summary = OUTDIR / "flat_region_summary.txt"
with open(summary, "w") as f:
    f.write(f"R1: {len(r1)} plates, B_rep = {B1:.6f} T\n")
    f.write(f"R2: {len(r2)} plates, B_rep = {B2:.6f} T\n")
    f.write(f"R3: {len(r3)} plates, B_rep = {B3:.6f} T\n")
    f.write(f"Full TMS flat-region thickness-weighted B_rep = {BF:.6f} T\n")

print("\nSaved:")
print(png)
print(csv)
print(summary)
