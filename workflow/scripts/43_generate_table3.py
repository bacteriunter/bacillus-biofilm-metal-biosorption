#!/usr/bin/env python3
"""
Module 43 - Table 3

Published phenotypic evidence and genomic marker
detection in selected Bacillus-related strains.

Uses existing module28 outputs.
No additional homology searches or statistical tests.

Outputs:
  results/tables/module43/table3_phenotype_genome_evidence.tsv
  results/tables/module43/table3_phenotype_genome_evidence.csv
"""

from pathlib import Path
import pandas as pd


OUT = Path("results/tables/module43")

MAIN = Path(
    "results/functional/module28/gene_presence.tsv"
)

EXTERNAL = Path(
    "results/functional/external_phenotypes/"
    "module28/gene_presence.tsv"
)

GROUPS = {
    "biofilm_regulation": "Biofilm regulation",
    "eps_matrix": "EPS matrix",
    "protein_matrix": "Protein matrix",
    "poly_gamma_glutamate": "Poly-gamma-glutamate",
    "cell_wall_polymers": "Cell-wall polymers",
    "metal_homeostasis_not_biosorption":
        "Metal homeostasis",
}

# Previously established strain-to-assembly links.
# G5 is deliberately excluded from genomic comparisons.
RECORDS = [
    {
        "Strain": "AB211",
        "Species": "Priestia aryabhattai",
        "Assembly": "GCF_001858395.1",
        "Source": "external",
        "Phenotype": "Biofilm formation",
        "Experimental evidence":
            "Crystal violet OD570; 48, 72 and 96 h; M9",
        "Metal": "Not evaluated in cited assay",
        "Reference":
            "Bhattacharyya et al. (2017)",
        "DOI": "10.3389/fmicb.2017.00411",
    },
    {
        "Strain": "AR156",
        "Species": "Bacillus cereus",
        "Assembly": "GCA_001880305.1",
        "Source": "external",
        "Phenotype": "Biofilm formation",
        "Experimental evidence":
            "Crystal violet OD595; 48 h; "
            "LBGM with root-exudate treatments",
        "Metal": "Not evaluated in cited assay",
        "Reference": "Wang et al. (2019)",
        "DOI": "10.3389/fmicb.2019.00098",
    },
    {
        "Strain": "C10-4",
        "Species": "Bacillus altitudinis",
        "Assembly": "GCF_028462265.1",
        "Source": "main",
        "Phenotype": "Cd(II) biosorption",
        "Experimental evidence":
            "Reported capacity: 3.31 mg/g "
            "fresh biomass",
        "Metal": "Cd(II)",
        "Reference":
            "Published biosorption study",
        "DOI": "Not recorded",
    },
    {
        "Strain": "G5",
        "Species": "Bacillus cereus",
        "Assembly": "Not verified",
        "Source": "none",
        "Phenotype": "Biofilm formation",
        "Experimental evidence":
            "Crystal violet OD570; 20 h; "
            "5% LB; mono-/coculture",
        "Metal": "Not evaluated in cited assay",
        "Reference": "Cheng et al. (2014)",
        "DOI": "10.1371/journal.pone.0094220",
    },
]


def load_matrix(path):
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(path, sep="\t")

    if "genome" not in df.columns:
        raise ValueError(
            f"Missing genome column: {path}"
        )

    if df["genome"].duplicated().any():
        raise ValueError(
            f"Duplicated genome identifiers: {path}"
        )

    marker_columns = [
        col for col in df.columns
        if "::" in col
    ]

    if len(marker_columns) != 53:
        raise ValueError(
            f"Expected 53 markers in {path}; "
            f"found {len(marker_columns)}"
        )

    for col in marker_columns:
        values = pd.to_numeric(
            df[col], errors="raise"
        )

        if not values.isin([0, 1]).all():
            raise ValueError(
                f"Nonbinary marker values: {col}"
            )

        df[col] = values.astype(int)

    return df.set_index("genome")


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    main_df = load_matrix(MAIN)
    ext_df = load_matrix(EXTERNAL)

    if len(main_df) != 638:
        raise ValueError(
            f"Expected 638 main genomes; "
            f"found {len(main_df)}"
        )

    rows = []

    for record in RECORDS:
        row = {
            key: value
            for key, value in record.items()
            if key != "Source"
        }

        source = record["Source"]
        accession = record["Assembly"]

        if source == "main":
            matrix = main_df
        elif source == "external":
            matrix = ext_df
        else:
            matrix = None

        if matrix is not None:
            if accession not in matrix.index:
                raise ValueError(
                    f"Assembly not found: {accession}"
                )

            values = matrix.loc[accession]

            for group, name in GROUPS.items():
                columns = [
                    col for col in matrix.columns
                    if col.startswith(group + "::")
                ]

                row[name + " (detected/evaluated)"] = (
                    f"{int(values[columns].sum())}"
                    f"/{len(columns)}"
                )
        else:
            for name in GROUPS.values():
                row[name + " (detected/evaluated)"] = (
                    "Not available"
                )

        row["Genome-phenotype linkage"] = (
            "Exact strain assembly"
            if matrix is not None
            else "No verified genome"
        )

        rows.append(row)

    table = pd.DataFrame(rows)

    table.to_csv(
        OUT / "table3_phenotype_genome_evidence.tsv",
        sep="\t",
        index=False,
    )

    table.to_csv(
        OUT / "table3_phenotype_genome_evidence.csv",
        index=False,
    )

    print("\nTABLE 3")
    print(table.to_string(index=False))

    print("\nPASS: Module 43")
    print(f"Output: {OUT}")


if __name__ == "__main__":
    main()
