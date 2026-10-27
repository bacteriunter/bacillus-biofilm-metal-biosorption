#!/usr/bin/env python3
"""
Module 44 - Table 4

Global and species-conditioned marker associations.

Uses existing module36 results only.
No new statistical tests.
"""

from pathlib import Path
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment


INPUT = Path(
    "results/statistics/module36/"
    "marker_associations.tsv"
)

OUT = Path("results/tables/module44")
OUT.mkdir(parents=True, exist_ok=True)


def main():

    df = pd.read_csv(INPUT, sep="\t")

    required = [
        "marker_matrix",
        "marker_metal",
        "odds_ratio_global",
        "q_fisher_bh",
        "p_conditional_species",
        "q_conditional_bh",
    ]

    missing = set(required) - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing columns: {sorted(missing)}"
        )

    if len(df) != 420:
        raise ValueError(
            f"Expected 420 pairs, found {len(df)}"
        )

    global_q = pd.to_numeric(
        df["q_fisher_bh"], errors="coerce"
    )

    conditional_p = pd.to_numeric(
        df["p_conditional_species"],
        errors="coerce"
    )

    conditional_q = pd.to_numeric(
        df["q_conditional_bh"],
        errors="coerce"
    )

    n_global = int((global_q < 0.05).sum())
    n_conditional = int(conditional_p.notna().sum())
    n_conditional_sig = int(
        (conditional_q < 0.05).sum()
    )

    rows = [
        (
            "Marker pairs evaluated",
            len(df),
        ),
        (
            "Pairs with available global adjusted p-values",
            int(global_q.notna().sum()),
        ),
        (
            "Globally significant associations "
            "(Fisher test, BH q < 0.05)",
            n_global,
        ),
        (
            "Pairs with evaluable species-conditioned tests",
            n_conditional,
        ),
        (
            "Significant species-conditioned associations "
            "(BH q < 0.05)",
            n_conditional_sig,
        ),
    ]

    table = pd.DataFrame(
        rows,
        columns=["Characteristic", "Value"]
    )

    table.to_csv(
        OUT / "table4_marker_associations.tsv",
        sep="\t",
        index=False,
    )

    table.to_csv(
        OUT / "table4_marker_associations.csv",
        index=False,
    )

    wb = Workbook()
    ws = wb.active
    ws.title = "Table 4"

    ws.append(list(table.columns))

    for record in table.itertuples(
        index=False, name=None
    ):
        ws.append(list(record))

    for cell in ws[1]:
        cell.font = Font(bold=True)
        cell.alignment = Alignment(
            horizontal="center"
        )

    ws.column_dimensions["A"].width = 72
    ws.column_dimensions["B"].width = 18

    wb.save(
        OUT / "table4_marker_associations.xlsx"
    )

    print("\nTABLE 4")
    print(table.to_string(index=False))

    print("\nPASS: Module 44")
    print(f"Output: {OUT}")


if __name__ == "__main__":
    main()
