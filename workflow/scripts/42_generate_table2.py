#!/usr/bin/env python3
"""
Module 42 - Table 2

Detection prevalence and variance partitioning
of functional markers in Bacillus.

Input:
    results/functional/module32/gene_variance_partition.tsv

Outputs:
    results/tables/module42/table2_functional_variability.tsv
    results/tables/module42/table2_functional_variability.csv
    results/tables/module42/table2_marker_details.tsv

Descriptive summary of existing results only.
No statistical tests or genomic analyses.
"""

from pathlib import Path
import numpy as np
import pandas as pd


INPUT = Path(
    "results/functional/module32/"
    "gene_variance_partition.tsv"
)

OUT = Path("results/tables/module42")

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


def median_or_na(values):
    values = pd.to_numeric(
        values, errors="coerce"
    ).dropna()

    if values.empty:
        return "NA"

    return f"{values.median():.2f}"


def main():

    if not INPUT.exists():
        raise FileNotFoundError(INPUT)

    OUT.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(INPUT, sep="\t")

    required = {
        "label",
        "genomes",
        "genomes_with_hit",
        "prevalence",
        "total_binary_variance",
        "between_fraction_of_total",
        "within_fraction_of_total",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing columns: {sorted(missing)}"
        )

    if len(df) != 53:
        raise ValueError(
            f"Expected 53 markers; found {len(df)}"
        )

    if df["label"].duplicated().any():
        raise ValueError(
            "Duplicated marker labels"
        )

    if not (df["genomes"] == 638).all():
        raise ValueError(
            "Unexpected genome count"
        )

    df["group"] = df["label"].str.split(
        "::", n=1
    ).str[0]

    df["gene"] = df["label"].str.split(
        "::", n=1
    ).str[1]

    df["variable"] = (
        df["total_binary_variance"] > 0
    )

    df["marker_prevalence_percent"] = (
        df["prevalence"] * 100
    )

    df["between_percent"] = (
        df["between_fraction_of_total"] * 100
    )

    df["within_percent"] = (
        df["within_fraction_of_total"] * 100
    )

    expected_groups = {
        group for group, _ in GROUPS
    }

    actual_groups = set(df["group"])

    if actual_groups != expected_groups:
        raise ValueError(
            "Functional groups differ from expected: "
            f"{sorted(actual_groups ^ expected_groups)}"
        )

    rows = []

    for group, display_name in GROUPS:

        subset = df.loc[
            df["group"] == group
        ].copy()

        variable = subset.loc[
            subset["variable"]
        ]

        invariant = subset.loc[
            ~subset["variable"]
        ]

        if (
            len(variable) + len(invariant)
            != len(subset)
        ):
            raise ValueError(
                f"Inconsistent counts: {group}"
            )

        # For variable markers, both fractions
        # must be defined and sum to 100%.
        if not variable.empty:

            fractions = (
                variable["between_fraction_of_total"]
                + variable["within_fraction_of_total"]
            )

            if not np.allclose(
                fractions.to_numpy(dtype=float),
                1.0,
                atol=1e-6,
            ):
                raise ValueError(
                    f"Variance fractions inconsistent: {group}"
                )

        rows.append({
            "Functional group": display_name,
            "Markers evaluated (n)": len(subset),
            "Invariant markers (n)": len(invariant),
            "Variable markers (n)": len(variable),
            "Median marker prevalence (%)":
                median_or_na(
                    subset["marker_prevalence_percent"]
                ),
            "Between-species variance (%)":
                median_or_na(
                    variable["between_percent"]
                ),
            "Within-species variance (%)":
                median_or_na(
                    variable["within_percent"]
                ),
        })

    table = pd.DataFrame(rows)

    # Overall row: marker-level medians,
    # not averages of group medians.
    variable_all = df.loc[df["variable"]]

    overall = pd.DataFrame([{
        "Functional group": "Overall",
        "Markers evaluated (n)": len(df),
        "Invariant markers (n)":
            int((~df["variable"]).sum()),
        "Variable markers (n)":
            int(df["variable"].sum()),
        "Median marker prevalence (%)":
            median_or_na(
                df["marker_prevalence_percent"]
            ),
        "Between-species variance (%)":
            median_or_na(
                variable_all["between_percent"]
            ),
        "Within-species variance (%)":
            median_or_na(
                variable_all["within_percent"]
            ),
    }])

    table = pd.concat(
        [table, overall],
        ignore_index=True,
    )

    table.to_csv(
        OUT / "table2_functional_variability.tsv",
        sep="\t",
        index=False,
    )

    table.to_csv(
        OUT / "table2_functional_variability.csv",
        index=False,
    )

    details = df[[
        "group",
        "gene",
        "label",
        "genomes",
        "genomes_with_hit",
        "prevalence",
        "total_binary_variance",
        "between_fraction_of_total",
        "within_fraction_of_total",
    ]].copy()

    details.to_csv(
        OUT / "table2_marker_details.tsv",
        sep="\t",
        index=False,
    )

    print("\nTABLE 2")
    print(table.to_string(index=False))

    print("\nPASS: Module 42")
    print(f"Output: {OUT}")


if __name__ == "__main__":
    main()
