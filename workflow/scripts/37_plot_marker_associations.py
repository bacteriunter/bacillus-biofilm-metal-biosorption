#!/usr/bin/env python3
"""
Module 37 - Figure 3

Visualization of extracellular matrix and metal
homeostasis marker associations.

Input:
    results/statistics/module36/marker_associations.tsv

Outputs:
    results/figures/module37/figure_marker_associations.png
    results/figures/module37/figure_marker_associations.pdf
    results/figures/module37/figure_marker_associations.svg

Visualization only. No statistical tests are repeated.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import (
    LinearSegmentedColormap,
    TwoSlopeNorm,
)
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D


INPUT = Path(
    "results/statistics/module36/marker_associations.tsv"
)

OUTPUT = Path("results/figures/module37")

GROUPS = [
    ("eps_matrix", "EPS matrix"),
    ("protein_matrix", "Protein matrix"),
    ("poly_gamma_glutamate", "Poly-gamma-glutamate"),
    ("cell_wall_polymers", "Cell-wall polymers"),
]

# Negative / zero / positive log2 odds ratio
COLORS = [
    "#572238",
    "#D997A5",
    "#FFFFFF",
    "#E7A6B2",
    "#9D294F",
]


def main():

    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    OUTPUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(INPUT, sep="\t")

    required = {
        "marker_matrix",
        "marker_metal",
        "odds_ratio_global",
        "q_fisher_bh",
        "p_conditional_species",
        "q_conditional_bh",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing columns: {sorted(missing)}"
        )

    df["group"] = df["marker_matrix"].str.split(
        "::", n=1
    ).str[0]

    matrix_order = []

    for group, label in GROUPS:

        markers = sorted(
            df.loc[
                df["group"] == group,
                "marker_matrix",
            ].unique()
        )

        matrix_order.extend(markers)

    metal_order = sorted(
        df["marker_metal"].unique()
    )

    if len(matrix_order) != 30:
        raise ValueError(
            f"Expected 30 matrix markers; "
            f"found {len(matrix_order)}"
        )

    if len(metal_order) != 14:
        raise ValueError(
            f"Expected 14 metal markers; "
            f"found {len(metal_order)}"
        )

    if len(df) != 420:
        raise ValueError(
            f"Expected 420 pairs; found {len(df)}"
        )

    if df.duplicated(
        ["marker_matrix", "marker_metal"]
    ).any():
        raise ValueError(
            "Duplicated marker pairs"
        )

    def pivot(column):

        return df.pivot(
            index="marker_matrix",
            columns="marker_metal",
            values=column,
        ).reindex(
            index=matrix_order,
            columns=metal_order,
        )

    odds = pivot("odds_ratio_global")
    q_global = pivot("q_fisher_bh")
    p_conditional = pivot("p_conditional_species")

    raw = odds.to_numpy(dtype=float)

    finite_positive = raw[
        np.isfinite(raw) & (raw > 0)
    ]

    if finite_positive.size == 0:
        raise ValueError(
            "No positive finite odds ratios"
        )

    finite_log = np.log2(finite_positive)

    limit = max(
        2.0,
        float(
            np.percentile(
                np.abs(finite_log),
                95,
            )
        ),
    )

    log_or = np.full(
        raw.shape,
        np.nan,
        dtype=float,
    )

    positive = np.isfinite(raw) & (raw > 0)
    zero = raw == 0
    infinite = np.isposinf(raw)

    log_or[positive] = np.log2(raw[positive])
    log_or[zero] = -limit
    log_or[infinite] = limit

    log_or = np.clip(
        log_or,
        -limit,
        limit,
    )

    masked = np.ma.masked_invalid(log_or)

    cmap = LinearSegmentedColormap.from_list(
        "wine_diverging",
        COLORS,
        N=256,
    )

    cmap.set_bad("#E3E3E3")

    norm = TwoSlopeNorm(
        vmin=-limit,
        vcenter=0,
        vmax=limit,
    )

    # Three independent axes:
    # heatmap / group labels / colorbar.
    # Group labels cannot overlap the colorbar.

    fig = plt.figure(
        figsize=(12, 13),
        facecolor="white",
    )

    gs = GridSpec(
        1,
        3,
        figure=fig,
        width_ratios=[16, 2.5, 0.8],
        wspace=0.12,
    )

    ax = fig.add_subplot(gs[0, 0])
    ax_groups = fig.add_subplot(gs[0, 1])
    cax = fig.add_subplot(gs[0, 2])

    image = ax.imshow(
        masked,
        cmap=cmap,
        norm=norm,
        aspect="auto",
        interpolation="nearest",
    )

    ax.set_xticks(
        np.arange(len(metal_order))
    )

    ax.set_xticklabels(
        [
            marker.split("::", 1)[1]
            for marker in metal_order
        ],
        rotation=45,
        ha="left",
        fontsize=9,
        fontstyle="italic",
    )

    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")

    ax.set_xlabel(
        "Metal homeostasis-associated markers",
        fontsize=11,
        labelpad=15,
    )

    ax.set_yticks(
        np.arange(len(matrix_order))
    )

    ax.set_yticklabels(
        [
            marker.split("::", 1)[1]
            for marker in matrix_order
        ],
        fontsize=9,
        fontstyle="italic",
    )

    ax.set_ylabel(
        "Extracellular matrix-associated markers",
        fontsize=11,
    )

    ax.tick_params(
        axis="both",
        length=0,
    )

    q_values = q_global.to_numpy(dtype=float)

    conditional_values = (
        p_conditional.to_numpy(dtype=float)
    )

    for i in range(len(matrix_order)):
        for j in range(len(metal_order)):

            if (
                np.isfinite(q_values[i, j])
                and q_values[i, j] < 0.05
            ):

                ax.plot(
                    j,
                    i,
                    marker="o",
                    markersize=3,
                    color="black",
                    linestyle="none",
                )

            if np.isfinite(
                conditional_values[i, j]
            ):

                ax.add_patch(
                    Rectangle(
                        (j - 0.45, i - 0.45),
                        0.9,
                        0.9,
                        fill=False,
                        edgecolor="black",
                        linewidth=0.65,
                    )
                )

    # Dedicated group-label axis.
    ax_groups.set_xlim(0, 1)
    ax_groups.set_ylim(
        len(matrix_order) - 0.5,
        -0.5,
    )

    ax_groups.axis("off")

    start = 0

    for group, label in GROUPS:

        count = sum(
            marker.startswith(group + "::")
            for marker in matrix_order
        )

        center = start + (count - 1) / 2

        # Short labels for narrow functional groups.
        display_labels = {
            "eps_matrix": "EPS matrix",
            "protein_matrix": "Protein matrix",
            "poly_gamma_glutamate": "Poly-γ-glutamate",
            "cell_wall_polymers": "Cell-wall polymers",
        }

        ax_groups.text(
            0.5,
            center,
            display_labels[group],
            ha="center",
            va="center",
            fontsize=7.5,
            fontweight="bold",
            rotation=90,
            rotation_mode="anchor",
        )

        start += count

        if start < len(matrix_order):

            ax.axhline(
                start - 0.5,
                color="#FFFFFF",
                linewidth=1.5,
            )

    colorbar = fig.colorbar(
        image,
        cax=cax,
    )

    colorbar.set_label(
        r"Global $\log_2$ odds ratio",
        fontsize=11,
        labelpad=10,
    )

    colorbar.ax.tick_params(
        labelsize=9,
    )

    legend_elements = [
        Line2D(
            [0], [0],
            marker="o",
            color="none",
            markerfacecolor="black",
            markeredgecolor="black",
            markersize=5,
            label="Global Fisher test: BH q < 0.05",
        ),
        Patch(
            facecolor="none",
            edgecolor="black",
            label="Species-conditioned test evaluable",
        ),
        Patch(
            facecolor="#E3E3E3",
            edgecolor="none",
            label="Global odds ratio not estimable",
        ),
    ]

    fig.legend(
        handles=legend_elements,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.025),
        ncol=2,
        frameon=False,
        fontsize=9,
    )

    fig.subplots_adjust(
        left=0.14,
        right=0.91,
        top=0.84,
        bottom=0.13,
    )

    for extension in ["png", "pdf", "svg"]:

        fig.savefig(
            OUTPUT
            / f"figure_marker_associations.{extension}",
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
        )

    plt.close(fig)

    n_global = int(
        (df["q_fisher_bh"] < 0.05).sum()
    )

    n_conditional = int(
        (df["q_conditional_bh"] < 0.05).sum()
    )

    n_evaluable = int(
        df["p_conditional_species"].notna().sum()
    )

    print("PASS: Module 37")
    print(f"Marker pairs: {len(df)}")
    print(f"Global significant: {n_global}")
    print(
        "Species-conditioned significant: "
        f"{n_conditional}/{n_evaluable}"
    )
    print(f"Output directory: {OUTPUT}")


if __name__ == "__main__":
    main()
