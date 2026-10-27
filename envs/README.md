# Computational environments

## Scientific analysis environment

The scientific Python libraries recorded in
`verified_versions.tsv` were inspected using the Python
3.10.18 interpreter associated with the `docking` environment.

This file also records command-line bioinformatics tools
available through the system PATH during inspection.

## Workflow management environment

The `bmetal_workflow` Conda environment contains:

- Python 3.10.22
- Snakemake 7.32.4

Versions are recorded in `snakemake_versions.tsv`.

## Environment specification

`environment.yml` is a proposed dependency specification.
It has not yet been validated through a clean installation
and should not be treated as an exact environment lockfile.

## Reproducibility

Scientific dependencies and workflow-management dependencies
are documented separately to avoid conflating environments.

Exact software versions should be preserved when reproducing
the analyses. Differences in versions may affect results.
