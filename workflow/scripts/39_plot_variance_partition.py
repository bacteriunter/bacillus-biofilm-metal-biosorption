#!/usr/bin/env python3
"""
Module 39 - Figure 2

Descriptive partition of binary marker-detection variance
into between-species and within-species components.

Input:
    results/functional/module32/gene_variance_partition.tsv

Outputs:
    results/figures/module39/figure2_variance_partition.png
    results/figures/module39/figure2_variance_partition.pdf
    results/figures/module39/figure2_variance_partition.svg

This script performs visualization only.
No variance decomposition or statistical tests are repeated.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


INPUT = Path(
    "results/functional/module32/gene_variance_partition.tsv"
)

OUTPUT = Path(
    "results/figures/module39"
)

GROUPS = [
    ("biofilm_regulation", "A. Biofilm regulation"),
    ("eps_matrix", "B. EPS matrix"),
    ("protein_matrix", "C. Protein matrix"),
    (
        "poly_gamma_glutamate",
        "D. Poly-gamma-glutamate",
    ),
    ("cell_wall_polymers", "E. Cell-wall polymers"),
    (
        "metal_homeostasis_not_biosorption",
        "F. Metal homeostasis",
    ),
]

COLOR_BETWEEN = "#572238"
COLOR_WITHIN = "#D997A5"
COLOR_INVARIANT = "#D9D9D9"


def load_data():

    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    df = pd.read_csv(
        INPUT,
        sep="\t",
    )

    required = [
        "label",
        "genomes",
        "total_binary_variance",
        "between_species_component",
        "within_species_component",
        "between_fraction_of_total",
        "within_fraction_of_total",
    ]

    missing = [
        column
        for column in required
        if column not in df.columns
    ]

    if missing:
        raise ValueError(
            f"Missing columns: {missing}"
        )

    if len(df) != 53:
        raise ValueError(
            f"Expected 53 markers, found {len(df)}"
        )

    if df["label"].duplicated().any():
        raise ValueError(
            "Duplicated functional markers"
        )

    if not (df["genomes"] == 638).all():
        raise ValueError(
            "Unexpected genome counts"
        )

    df["group"] = df["label"].str.split(
        "::", n=1
    ).str[0]

    df["gene"] = df["label"].str.split(
        "::", n=1
    ).str[1]

    numeric_columns = [
        "total_binary_variance",
        "between_species_component",
        "within_species_component",
        "between_fraction_of_total",
        "within_fraction_of_total",
    ]

    for column in numeric_columns:
        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        )

    if df[
        [
            "total_binary_variance",
            "between_species_component",
            "within_species_component",
        ]
    ].isna().any().any():
        raise ValueError(
            "Missing variance components"
        )

    if (
        df["total_binary_variance"] < -1e-10
    ).any():
        raise ValueError(
            "Negative total variance"
        )

    # Invariant markers have undefined variance fractions.
    df["invariant"] = np.isclose(
        df["total_binary_variance"],
        0.0,
        atol=1e-12,
    )

    variable = df.loc[~df["invariant"]]

    if variable[
        [
            "between_fraction_of_total",
            "within_fraction_of_total",
        ]
    ].isna().any().any():
        raise ValueError(
            "Missing variance fractions for variable markers"
        )

    fractions_sum = (
        variable["between_fraction_of_total"]
        + variable["within_fraction_of_total"]
    )

    if not np.allclose(
        fractions_sum,
        1.0,
        atol=1e-5,
    ):
        raise ValueError(
            "Variance fractions do not sum to 1"
        )

    if (
        (
            variable[
                [
                    "between_fraction_of_total",
                    "within_fraction_of_total",
                ]
            ] < -1e-8
        ).any().any()
    ):
        raise ValueError(
            "Negative variance fractions"
        )

    return df


def main():

    OUTPUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    df = load_data()

    fig, axes = plt.subplots(
        3,
        2,
        figsize=(13, 14),
    )

    for ax, (group, title) in zip(
        axes.flat,
        GROUPS,
    ):

        subset = (
            df[df["group"] == group]
            .sort_values("gene")
            .reset_index(drop=True)
        )

        if subset.empty:
            raise ValueError(
                f"No markers for group: {group}"
            )

        y = np.arange(len(subset))

        between = (
            subset["between_fraction_of_total"]
            .fillna(0)
            .to_numpy()
            * 100
        )

        within = (
            subset["within_fraction_of_total"]
            .fillna(0)
            .to_numpy()
            * 100
        )

        invariant = subset["invariant"].to_numpy()

        # Plot only variable markers as stacked bars.
        ax.barh(
            y[~invariant],
            between[~invariant],
            height=0.65,
            color=COLOR_BETWEEN,
            edgecolor="white",
            linewidth=0.5,
        )

        ax.barh(
            y[~invariant],
            within[~invariant],
            left=between[~invariant],
            height=0.65,
            color=COLOR_WITHIN,
            edgecolor="white",
            linewidth=0.5,
        )

        # Invariant markers: visually distinct and
        # explicitly labelled, without inventing fractions.
        for i in np.flatnonzero(invariant):

            ax.barh(
                i,
                100,
                height=0.65,
                color=COLOR_INVARIANT,
                edgecolor="white",
                linewidth=0.5,
            )

            ax.text(
                50,
                i,
                "Invariant",
                ha="center",
                va="center",
                fontsize=8,
                color="#555555",
            )

        ax.set_yticks(y)

        ax.set_yticklabels(
            subset["gene"],
            fontsize=9,
            fontstyle="italic",
        )

        ax.invert_yaxis()

        ax.set_xlim(0, 100)

        ax.set_xticks(
            [0, 25, 50, 75, 100]
        )

        ax.set_xlabel(
            "Fraction of total binary variance (%)",
            fontsize=9,
        )

        ax.set_title(
            title,
            loc="left",
            fontsize=11,
            fontweight="bold",
            pad=10,
        )

        ax.grid(
            axis="x",
            alpha=0.15,
            linewidth=0.6,
        )

        ax.set_axisbelow(True)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.tick_params(
            axis="both",
            length=0,
        )

    legend = [
        Patch(
            facecolor=COLOR_BETWEEN,
            label="Between-species component",
        ),
        Patch(
            facecolor=COLOR_WITHIN,
            label="Within-species component",
        ),
        Patch(
            facecolor=COLOR_INVARIANT,
            label="Invariant marker",
        ),
    ]

    fig.legend(
        handles=legend,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.025),
        ncol=3,
        frameon=False,
        fontsize=10,
    )

    fig.suptitle(
        "Partition of functional-marker variability "
        r"across $\it{Bacillus}$ genomes",
        fontsize=14,
        y=0.995,
    )

    fig.tight_layout(
        rect=[0.03, 0.065, 0.98, 0.97],
        h_pad=2.5,
        w_pad=2.0,
    )

    for extension in [
        "png",
        "pdf",
        "svg",
    ]:

        destination = OUTPUT / (
            "figure2_variance_partition."
            + extension
        )

        fig.savefig(
            destination,
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
        )

    plt.close(fig)

    n_invariant = int(
        df["invariant"].sum()
    )

    print("PASS: Module 39")
    print(f"Markers: {len(df)}")
    print(f"Variable markers: {len(df) - n_invariant}")
    print(f"Invariant markers: {n_invariant}")
    print(f"Output directory: {OUTPUT}")

    for extension in ["png", "pdf", "svg"]:
        print(
            OUTPUT
            / f"figure2_variance_partition.{extension}"
        )


if __name__ == "__main__":
    main()
