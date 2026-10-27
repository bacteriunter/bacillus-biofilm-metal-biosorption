#!/usr/bin/env python3
"""
Module 40 - Figure 4

Conceptual framework linking aquatic phytoremediation,
Bacillus-associated extracellular matrix systems,
and metal-related genomic markers.

This is a conceptual illustration, not an experimental
result or a quantitative model.

Outputs: PNG, PDF, SVG.
"""

from pathlib import Path
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import (
    FancyBboxPatch, FancyArrowPatch,
    Circle, Ellipse, PathPatch
)
from matplotlib.path import Path as MplPath

OUT = Path("results/figures/module40")
OUT.mkdir(parents=True, exist_ok=True)

WINE = "#572238"
ROSE = "#D997A5"
PALE = "#F5DCE1"
DARK = "#333333"
GREEN = "#648773"
WATER = "#E8F2F1"
GRAY = "#6B6B6B"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

fig, ax = plt.subplots(figsize=(15, 8.5))
fig.patch.set_facecolor("white")
ax.set_xlim(0, 15)
ax.set_ylim(0.9, 8.5)
ax.axis("off")


def box(x, y, w, h, title, lines,
        edge=WINE, fill="white"):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.14,rounding_size=0.12",
        linewidth=1.5,
        edgecolor=edge,
        facecolor=fill,
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2, y + h - 0.30,
        title,
        ha="center", va="center",
        fontsize=11, fontweight="bold",
        color=edge,
    )
    ax.text(
        x + w / 2, y + h / 2 - 0.12,
        lines,
        ha="center", va="center",
        fontsize=9.2, color=DARK,
        linespacing=1.5,
    )


def arrow(x1, y1, x2, y2,
          color=WINE, dashed=False):
    patch = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle="-|>",
        mutation_scale=14,
        linewidth=1.6,
        linestyle="--" if dashed else "-",
        color=color,
        connectionstyle="arc3,rad=0",
    )
    ax.add_patch(patch)


# A. Aquatic environment
ax.add_patch(FancyBboxPatch(
    (0.45, 2.1), 5.0, 5.6,
    boxstyle="round,pad=0.12",
    facecolor=WATER,
    edgecolor="#B5CECA",
    linewidth=1.2,
))

ax.text(
    2.95, 7.35,
    "A. Contaminated water and aquatic plant",
    ha="center", fontsize=11,
    fontweight="bold", color=DARK,
)

# Water line
ax.plot(
    [0.75, 5.15], [5.65, 5.65],
    color="#83B7B1", linewidth=2.0,
)

# Stylized aquatic plant
ax.plot(
    [2.8, 2.8], [3.75, 6.55],
    color=GREEN, linewidth=4,
    solid_capstyle="round",
)

for x, y, angle in [
    (2.35, 6.0, 35),
    (3.25, 6.1, -35),
    (2.30, 6.55, 25),
    (3.30, 6.65, -25),
]:
    ax.add_patch(Ellipse(
        (x, y), 1.0, 0.35,
        angle=angle,
        facecolor="#A8C7B2",
        edgecolor=GREEN,
        linewidth=1,
    ))

# Roots
for dx in [-0.65, -0.30, 0, 0.30, 0.65]:
    vertices = [
        (2.8, 3.75),
        (2.8 + dx * 0.5, 3.05),
        (2.8 + dx, 2.65),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE3,
        MplPath.CURVE3,
    ]
    ax.add_patch(PathPatch(
        MplPath(vertices, codes),
        edgecolor=GREEN,
        facecolor="none",
        linewidth=2,
    ))


# Root-associated bacteria and extracellular matrix.
# Conceptual illustration, not an observed phenotype.
from matplotlib.patches import Arc

root_bacteria = [
    (2.49, 3.33, -35),
    (3.04, 3.28, 35),
    (2.20, 2.92, -40),
    (3.37, 2.91, 40),
    (2.75, 3.06, 75),
]

for bx, by, angle in root_bacteria:
    ax.add_patch(Ellipse(
        (bx, by), 0.30, 0.13,
        angle=angle,
        facecolor=WINE,
        edgecolor="#FFFFFF",
        linewidth=0.6,
        zorder=6,
    ))

# Translucent extracellular matrix envelope.
ax.add_patch(Ellipse(
    (2.79, 3.15),
    2.10, 1.35,
    facecolor=ROSE,
    edgecolor=WINE,
    linewidth=0.9,
    alpha=0.12,
    zorder=2,
))

ax.text(
    4.15, 2.78,
    "Biofilm / EPS",
    ha="center",
    va="center",
    fontsize=8.5,
    color=WINE,
)

# Dissolved metal symbols
for x, y in [
    (1.1, 4.8), (4.5, 4.8),
    (1.25, 3.45), (4.35, 3.15),
    (1.65, 2.65), (4.65, 4.0),
]:
    ax.add_patch(Circle(
        (x, y), 0.13,
        facecolor=ROSE,
        edgecolor=WINE,
        linewidth=1,
    ))

ax.text(
    2.95, 2.27,
    "Dissolved metal species",
    ha="center", fontsize=9,
    color=WINE,
)

# B. Plant-associated Bacillus
box(
    6.0, 5.65, 3.4, 1.55,
    r"B. Root-associated $\it{Bacillus}$",
    "Root colonization\nBiofilm formation\nExtracellular matrix",
    edge=WINE, fill="#FBF5F7",
)

box(
    6.0, 3.3, 3.4, 1.55,
    "Potential metal interactions",
    "Extracellular binding\nSurface retention\nMetal tolerance",
    edge=WINE, fill="#FBF5F7",
)

arrow(4.05, 3.45, 5.95, 6.15, color=GREEN)
arrow(7.7, 5.6, 7.7, 4.9)

# C. Genomic evidence

box(
    10.1, 5.65, 3.8, 1.55,
    "Extracellular matrix markers",
    "EPS-associated systems\nMatrix proteins\nPoly-gamma-glutamate\nCell-wall polymers",
    edge=WINE, fill=PALE,
)

box(
    10.1, 3.3, 3.8, 1.55,
    "Metal-related markers",
    "Metal homeostasis\nTransport and regulation\nStress-associated functions",
    edge=WINE, fill=PALE,
)

# Dashed arrows: inferred conceptual links
arrow(
    10.0, 6.4, 9.45, 6.4,
    color=GRAY, dashed=True,
)

arrow(
    10.0, 4.05, 9.45, 4.05,
    color=GRAY, dashed=True,
)

# Legend
ax.plot(
    [5.1, 5.6], [1.15, 1.15],
    color=WINE, linewidth=1.8,
)
ax.text(
    5.7, 1.15,
    "Biological process",
    va="center", fontsize=9,
)

ax.plot(
    [8.2, 8.7], [1.15, 1.15],
    color=GRAY, linewidth=1.8,
    linestyle="--",
)
ax.text(
    8.8, 1.15,
    "Proposed genomic connection",
    va="center", fontsize=9,
)

for ext in ["png", "pdf", "svg"]:
    fig.savefig(
        OUT / f"figure4_phytoremediation_framework.{ext}",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

plt.close(fig)

print("PASS: Module 40")
print(f"Output: {OUT}")
