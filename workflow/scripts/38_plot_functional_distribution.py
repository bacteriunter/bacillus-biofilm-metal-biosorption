#!/usr/bin/env python3
"""
Module 38 - Figure 1
Species-level functional marker prevalence.

Heatmap:
    123 species x 53 markers
    Six functional categories
    Prevalence scale: 0-100%

No statistical or genomic analyses are repeated.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle


INPUT = Path(
    "results/functional/module29/species_gene_prevalence.tsv"
)

OUTPUT = Path("results/figures/module38")

GROUPS = [
    ("biofilm_regulation", "Biofilm regulation"),
    ("eps_matrix", "EPS matrix"),
    ("protein_matrix", "Protein matrix"),
    ("poly_gamma_glutamate", "Poly-gamma-glutamate"),
    ("cell_wall_polymers", "Cell-wall polymers"),
    (
        "metal_homeostasis_not_biosorption",
        "Metal homeostasis",
    ),
]

COLORS = [
    "#FFFFFF",
    "#F5DCE1",
    "#D997A5",
    "#A34B64",
    "#572238",
]


def load_matrix():

    df = pd.read_csv(INPUT, sep="\t")
    df.columns = df.columns.str.strip()

    species_candidates = [
        "species",
        "species_name",
        "organism",
        "organism_name",
        "taxon",
    ]

    species_col = next(
        (
            c for c in species_candidates
            if c in df.columns
        ),
        None,
    )

    if species_col is None:
        raise ValueError(
            f"Species column not found: {list(df.columns)}"
        )

    marker_columns = [
        c for c in df.columns
        if "::" in c
    ]

    if marker_columns:

        matrix = df.set_index(
            species_col
        )[marker_columns]

    else:

        group_col = next(
            (
                c for c in [
                    "group",
                    "functional_group",
                    "category",
                    "system",
                ]
                if c in df.columns
            ),
            None,
        )

        gene_col = next(
            (
                c for c in [
                    "gene",
                    "marker",
                    "gene_name",
                    "marker_name",
                ]
                if c in df.columns
            ),
            None,
        )

        value_col = next(
            (
                c for c in [
                    "prevalence",
                    "gene_prevalence",
                    "prevalence_fraction",
                    "frequency",
                    "proportion",
                    "fraction",
                ]
                if c in df.columns
            ),
            None,
        )

        if gene_col is None or value_col is None:
            raise ValueError(
                f"Unrecognized format: {list(df.columns)}"
            )

        temp = df.copy()

        if group_col is not None:
            temp["marker_id"] = (
                temp[group_col].astype(str)
                + "::"
                + temp[gene_col].astype(str)
            )
        else:
            temp["marker_id"] = (
                temp[gene_col].astype(str)
            )

        if temp.duplicated(
            [species_col, "marker_id"]
        ).any():
            raise ValueError(
                "Duplicated species-marker combinations"
            )

        matrix = temp.pivot(
            index=species_col,
            columns="marker_id",
            values=value_col,
        )

    matrix = matrix.apply(
        pd.to_numeric,
        errors="raise",
    )

    if matrix.shape != (123, 53):
        raise ValueError(
            f"Expected (123, 53), found {matrix.shape}"
        )

    if matrix.isna().any().any():
        raise ValueError(
            "Missing prevalence values"
        )

    minimum = float(matrix.min().min())
    maximum = float(matrix.max().max())

    if minimum < 0:
        raise ValueError("Negative prevalence values")

    if maximum <= 1:
        matrix *= 100.0
    elif maximum > 100:
        raise ValueError(
            "Prevalence outside 0-100%"
        )

    return matrix


def main():

    OUTPUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    matrix = load_matrix()

    # Preserve the original heatmap organization:
    # markers grouped by functional category;
    # species ordered by mean prevalence.

    marker_order = []

    for group, label in GROUPS:
        marker_order.extend(
            sorted(
                c for c in matrix.columns
                if c.startswith(group + "::")
            )
        )

    if len(marker_order) != 53:
        raise ValueError(
            "Marker groups do not contain 53 markers"
        )

    matrix = matrix[marker_order]

    species_order = (
        matrix.mean(axis=1)
        .sort_values(
            ascending=False,
            kind="stable",
        )
        .index
    )

    matrix = matrix.loc[species_order]

    cmap = LinearSegmentedColormap.from_list(
        "rose_wine",
        COLORS,
        N=256,
    )

    norm = Normalize(
        vmin=0,
        vmax=100,
    )

    # Tall layout for 123 species.
    # Dedicated header axis prevents overlap
    # between functional category labels.

    fig = plt.figure(
        figsize=(16, 24),
        facecolor="white",
    )

    gs = GridSpec(
        2,
        2,
        figure=fig,
        width_ratios=[35, 1],
        height_ratios=[2.5, 20],
        wspace=0.07,
        hspace=0.015,
    )

    ax_header = fig.add_subplot(
        gs[0, 0]
    )

    ax = fig.add_subplot(
        gs[1, 0]
    )

    cax = fig.add_subplot(
        gs[1, 1]
    )

    values = matrix.to_numpy(
        dtype=float
    )

    image = ax.imshow(
        values,
        cmap=cmap,
        norm=norm,
        aspect="auto",
        interpolation="nearest",
        origin="upper",
    )

    ax.set_xticks(
        np.arange(len(marker_order))
    )

    ax.set_xticklabels(
        [
            marker.split("::", 1)[1]
            for marker in marker_order
        ],
        rotation=90,
        fontsize=8,
        fontstyle="italic",
    )

    ax.set_yticks(
        np.arange(len(matrix.index))
    )

    ax.set_yticklabels(
        matrix.index,
        fontsize=6.1,
        fontstyle="italic",
    )

    ax.tick_params(
        axis="both",
        length=0,
        pad=3,
    )

    ax.set_xlabel(
        "Functional markers",
        fontsize=11,
        labelpad=12,
    )


    # Header bands are aligned exactly
    # with the functional marker columns.

    ax_header.set_xlim(
        -0.5,
        len(marker_order) - 0.5,
    )

    ax_header.set_ylim(0, 1)

    ax_header.axis("off")

    start = 0

    for group, label in GROUPS:

        count = sum(
            marker.startswith(group + "::")
            for marker in marker_order
        )

        center = start + (count - 1) / 2

        ax_header.add_patch(
            Rectangle(
                (start - 0.5, 0.02),
                count,
                0.96,
                facecolor="#F2F0F1",
                edgecolor="white",
                linewidth=1.5,
            )
        )

        ax_header.text(
            center,
            0.50,
            label,
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
            rotation=90,
            clip_on=False,
        )

        start += count

        if start < len(marker_order):
            ax.axvline(
                start - 0.5,
                color="#FFFFFF",
                linewidth=1.5,
            )

    colorbar = fig.colorbar(
        image,
        cax=cax,
    )

    colorbar.set_ticks(
        [0, 25, 50, 75, 100]
    )

    colorbar.set_label(
        "Marker prevalence within species (%)",
        fontsize=10,
        labelpad=10,
    )

    colorbar.ax.tick_params(
        labelsize=9
    )

    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color("#555555")

    fig.subplots_adjust(
        left=0.22,
        right=0.94,
        top=0.975,
        bottom=0.065,
    )

    for extension in [
        "png",
        "pdf",
        "svg",
    ]:

        destination = OUTPUT / (
            "figure1_functional_distribution."
            + extension
        )

        fig.savefig(
            destination,
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
        )

    plt.close(fig)

    print("PASS: Module 38")
    print(f"Species: {matrix.shape[0]}")
    print(f"Markers: {matrix.shape[1]}")
    print("Palette: white / rose / wine")
    print(f"Output: {OUTPUT}")


if __name__ == "__main__":
    main()
