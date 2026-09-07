#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np

from tms_mapper.edep import iter_mapper_rows


parser = argparse.ArgumentParser()
parser.add_argument("mapper")
parser.add_argument("--z", type=float, default=0.0, help="Z slice in mm")
parser.add_argument("-o", "--output", default="scatter_xy.png")
args = parser.parse_args()

x = []
y = []

for row in iter_mapper_rows(args.mapper):
    X, Y, Z = row[:3]

    if np.isclose(Z, args.z, atol=1e-6):
        x.append(X)
        y.append(Y)

plt.figure(figsize=(9, 7))

plt.scatter(
    x,
    y,
    s=8
)

plt.xlabel("X [mm]")
plt.ylabel("Y [mm]")
plt.title(f"TMS XY grid coordinates at Z = {args.z:g} mm")

plt.axis("equal")
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig(args.output, dpi=200)

print(f"Saved {args.output}")
print(f"Points plotted: {len(x)}")

