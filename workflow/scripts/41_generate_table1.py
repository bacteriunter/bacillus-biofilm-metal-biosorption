#!/usr/bin/env python3
"""
Module 41 - Table 1

Characteristics of the Bacillus genomic dataset
and functional marker screening.

Uses existing outputs only.
No new genomic or statistical analyses.
"""

from pathlib import Path
import json

import numpy as np
import pandas as pd


COHORT = Path(
    "data/processed/functional_cohort_638.tsv"
)

SPECIES = Path(
    "results/functional/module29/species_summary.tsv"
)

SUMMARY = Path(
    "results/functional/module28/summary.json"
)

OUT = Path("results/tables/module41")
OUT.mkdir(parents=True, exist_ok=True)


def fmt(value, decimals=2):
    return f"{value:.{decimals}f}"


def main():

    cohort = pd.read_csv(COHORT, sep="\t")
    species = pd.read_csv(SPECIES, sep="\t")

    summary = json.loads(
        SUMMARY.read_text(encoding="utf-8")
    )

    if len(cohort) != 638:
        raise ValueError("Unexpected cohort size")

    if len(species) != 123:
        raise ValueError("Unexpected species count")

    if summary["genomes"] != len(cohort):
        raise ValueError(
            "Cohort and module28 genome counts differ"
        )

    if cohort["assembly_accession"].duplicated().any():
        raise ValueError("Duplicated assemblies")

    if cohort["species"].nunique() != len(species):
        raise ValueError(
            "Species counts are inconsistent"
        )

    if species["genomes"].sum() != len(cohort):
        raise ValueError(
            "Species-level genome totals differ"
        )

    completeness = pd.to_numeric(
        cohort["completeness"],
        errors="raise"
    )

    contamination = pd.to_numeric(
        cohort["contamination"],
        errors="raise"
    )

    genomes_per_species = species["genomes"]

    rows = []

    def add(section, characteristic, value):
        rows.append({
            "Section": section,
            "Characteristic": characteristic,
            "Value": value,
        })

    add(
        "Genomic dataset",
        "Genomes analyzed",
        len(cohort)
    )

    add(
        "Genomic dataset",
        "Species represented",
        cohort["species"].nunique()
    )

    add(
        "Genomic dataset",
        "Genomes per species, median (range)",
        (
            f"{fmt(genomes_per_species.median(), 1)} "
            f"({genomes_per_species.min()}–"
            f"{genomes_per_species.max()})"
        )
    )

    for level, count in (
        cohort["assembly_level"]
        .value_counts()
        .sort_index()
        .items()
    ):

        percentage = 100 * count / len(cohort)

        add(
            "Assembly characteristics",
            f"Assembly level: {level}",
            f"{count} ({fmt(percentage, 1)}%)"
        )

    add(
        "Assembly quality",
        "Completeness, median (IQR), %",
        (
            f"{fmt(completeness.median())} "
            f"({fmt(completeness.quantile(0.25))}–"
            f"{fmt(completeness.quantile(0.75))})"
        )
    )

    add(
        "Assembly quality",
        "Contamination, median (IQR), %",
        (
            f"{fmt(contamination.median())} "
            f"({fmt(contamination.quantile(0.25))}–"
            f"{fmt(contamination.quantile(0.75))})"
        )
    )

    add(
        "Functional screening",
        "Reference protein sequences",
        summary["reference_sequences"]
    )

    add(
        "Functional screening",
        "Functional marker labels",
        summary["gene_group_labels"]
    )

    add(
        "Functional screening",
        "Functional groups",
        6
    )

    add(
        "Functional screening",
        "Genomes with at least one filtered hit",
        summary["genomes_with_filtered_hits"]
    )

    add(
        "Functional screening",
        "Maximum E-value",
        summary["filter"]["max_evalue"]
    )

    add(
        "Functional screening",
        "Minimum query coverage",
        summary["filter"]["min_query_coverage"]
    )

    table = pd.DataFrame(rows)

    table.to_csv(
        OUT / "table1_genomic_dataset.tsv",
        sep="\t",
        index=False
    )

    table.to_csv(
        OUT / "table1_genomic_dataset.csv",
        index=False
    )

    audit = {
        "cohort_genomes": len(cohort),
        "species": cohort["species"].nunique(),
        "assembly_level_counts": (
            cohort["assembly_level"]
            .value_counts()
            .to_dict()
        ),
        "completeness": completeness.describe().to_dict(),
        "contamination": contamination.describe().to_dict(),
        "genomes_per_species": (
            genomes_per_species.describe().to_dict()
        ),
        "functional_summary": summary,
    }

    with open(
        OUT / "table1_audit.json",
        "w",
        encoding="utf-8"
    ) as handle:

        json.dump(
            audit,
            handle,
            indent=2,
            ensure_ascii=False
        )

    print("\nTABLE 1")
    print(table.to_string(index=False))

    print("\nPASS: Module 41")
    print(f"Output: {OUT}")


if __name__ == "__main__":
    main()
