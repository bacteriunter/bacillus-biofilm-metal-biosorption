#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests


GROUPS_MATRIX = [
    "eps_matrix",
    "protein_matrix",
    "poly_gamma_glutamate",
    "cell_wall_polymers",
]
GROUP_METAL = "metal_homeostasis_not_biosorption"


def load_metadata(path, genomes):
    meta = pd.read_csv(path, sep="\t", dtype=str)
    meta.columns = meta.columns.str.strip()

    genome_candidates = [
        "genome", "assembly_accession", "accession",
        "assembly", "genome_id", "assembly_id",
    ]
    species_candidates = [
        "species", "species_name", "organism_name",
        "organism", "taxon", "taxon_name",
    ]

    genome_col = next(
        (c for c in genome_candidates if c in meta.columns),
        None,
    )
    species_col = next(
        (c for c in species_candidates if c in meta.columns),
        None,
    )

    if genome_col is None or species_col is None:
        raise ValueError(
            "No se identificaron columnas de genoma/especie. "
            f"Columnas disponibles: {list(meta.columns)}"
        )

    meta = meta[[genome_col, species_col]].copy()
    meta.columns = ["genome", "species"]

    if meta.isna().any().any():
        raise ValueError("Existen identificadores o especies vacíos.")

    if meta["genome"].duplicated().any():
        raise ValueError(
            "Hay identificadores de genoma duplicados en metadatos."
        )

    meta = meta.set_index("genome")

    missing = genomes.difference(meta.index)
    if len(missing):
        raise ValueError(
            f"Faltan {len(missing)} genomas en metadatos. "
            f"Ejemplos: {list(missing[:5])}"
        )

    return meta.loc[genomes, "species"]


def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan)
    valid = np.isfinite(values)

    if valid.any():
        result[valid] = multipletests(
            values[valid], method="fdr_bh"
        )[1]

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matrix",
        default="results/functional/module28/gene_presence.tsv",
    )
    parser.add_argument(
        "--metadata",
        default="data/processed/functional_cohort_638.tsv",
    )
    parser.add_argument(
        "--output",
        default="results/statistics/module36",
    )
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20261026)
    args = parser.parse_args()

    if args.permutations < 99:
        raise ValueError("Se requieren al menos 99 simulaciones.")

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    matrix = pd.read_csv(args.matrix, sep="\t")
    matrix.columns = matrix.columns.str.strip()

    if "genome" not in matrix.columns:
        raise ValueError("La matriz no contiene la columna genome.")

    if matrix["genome"].duplicated().any():
        raise ValueError("La matriz contiene genomas duplicados.")

    matrix = matrix.set_index("genome")

    if matrix.shape[0] != 638:
        raise ValueError(
            f"Se esperaban 638 genomas; encontrados {matrix.shape[0]}."
        )

    selected = [
        c for c in matrix.columns
        if c.split("::", 1)[0] in
        GROUPS_MATRIX + [GROUP_METAL]
    ]

    matrix = matrix[selected].apply(
        pd.to_numeric, errors="raise"
    )

    if not matrix.isin([0, 1]).all().all():
        raise ValueError("La matriz debe ser binaria (0/1).")

    species = load_metadata(args.metadata, matrix.index)

    if species.nunique() != 123:
        raise ValueError(
            f"Se esperaban 123 especies; "
            f"encontradas {species.nunique()}."
        )

    matrix_markers = [
        c for c in selected
        if c.split("::", 1)[0] in GROUPS_MATRIX
    ]
    metal_markers = [
        c for c in selected
        if c.split("::", 1)[0] == GROUP_METAL
    ]

    if not matrix_markers or not metal_markers:
        raise ValueError("Faltan grupos funcionales requeridos.")

    # Índices por especie: cada simulación conserva
    # las frecuencias marginales dentro de cada especie.
    groups = [
        np.flatnonzero(species.to_numpy() == name)
        for name in species.unique()
    ]

    data = matrix.to_numpy(dtype=np.uint8)
    col_index = {name: i for i, name in enumerate(selected)}
    rng = np.random.default_rng(args.seed)

    records = []

    for marker_a in matrix_markers:
        a = data[:, col_index[marker_a]]

        for marker_b in metal_markers:
            b = data[:, col_index[marker_b]]

            n11 = int(np.sum((a == 1) & (b == 1)))
            n10 = int(np.sum((a == 1) & (b == 0)))
            n01 = int(np.sum((a == 0) & (b == 1)))
            n00 = int(np.sum((a == 0) & (b == 0)))

            table = [[n11, n10], [n01, n00]]

            global_testable = (
                (n11 + n10 > 0)
                and (n01 + n00 > 0)
                and (n11 + n01 > 0)
                and (n10 + n00 > 0)
            )

            if global_testable:
                odds, p_fisher = fisher_exact(
                    table, alternative="two-sided"
                )
            else:
                odds, p_fisher = np.nan, np.nan

            expected = 0.0
            simulated = np.zeros(
                args.permutations, dtype=np.int32
            )
            informative_species = 0

            for idx in groups:
                n = len(idx)
                na = int(a[idx].sum())
                nb = int(b[idx].sum())

                expected += na * nb / n

                # Solo las especies con variación en
                # ambos marcadores aportan aleatoriedad.
                if 0 < na < n and 0 < nb < n:
                    informative_species += 1
                    simulated += hypergeom.rvs(
                        M=n,
                        n=na,
                        N=nb,
                        size=args.permutations,
                        random_state=rng,
                    )
                else:
                    simulated += min(na, nb)

            if informative_species:
                deviation = abs(n11 - expected)
                simulated_deviation = abs(
                    simulated - expected
                )

                # Prueba bilateral basada en la distancia
                # respecto a la esperanza condicionada.
                extreme = int(np.count_nonzero(
                    simulated_deviation >= deviation - 1e-12
                ))
                p_conditional = (
                    extreme + 1
                ) / (args.permutations + 1)
            else:
                p_conditional = np.nan

            records.append({
                "marker_matrix": marker_a,
                "marker_metal": marker_b,
                "n11": n11,
                "n10": n10,
                "n01": n01,
                "n00": n00,
                "odds_ratio_global": odds,
                "p_fisher": p_fisher,
                "expected_codetection_species": expected,
                "observed_minus_expected": n11 - expected,
                "informative_species": informative_species,
                "p_conditional_species": p_conditional,
            })

    results = pd.DataFrame(records)

    results["q_fisher_bh"] = bh_adjust(
        results["p_fisher"]
    )
    results["q_conditional_bh"] = bh_adjust(
        results["p_conditional_species"]
    )

    results["significant_global"] = (
        results["q_fisher_bh"] < 0.05
    )
    results["significant_conditional"] = (
        results["q_conditional_bh"] < 0.05
    )

    results = results.sort_values(
        ["q_conditional_bh", "p_conditional_species"],
        na_position="last",
    )

    results.to_csv(
        out / "marker_associations.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary = {
        "genomes": int(matrix.shape[0]),
        "species": int(species.nunique()),
        "matrix_markers": len(matrix_markers),
        "metal_markers": len(metal_markers),
        "tested_pairs": len(results),
        "globally_testable_pairs": int(
            results["p_fisher"].notna().sum()
        ),
        "conditionally_testable_pairs": int(
            results["p_conditional_species"].notna().sum()
        ),
        "significant_global_fdr05": int(
            results["significant_global"].sum()
        ),
        "significant_conditional_fdr05": int(
            results["significant_conditional"].sum()
        ),
        "permutations": args.permutations,
        "seed": args.seed,
        "method": (
            "Fisher exact test and species-conditioned "
            "hypergeometric Monte Carlo test"
        ),
        "multiple_testing": (
            "Benjamini-Hochberg, separately by test family"
        ),
    }

    (out / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print(json.dumps(summary, indent=2))
    print(f"\nResultados: {out.resolve()}")


if __name__ == "__main__":
    main()
