# From bacterial biofilms to plant-assisted metal remediation

## Assessing genomic predictors of extracellular matrix formation and metal biosorption in Bacillus

This repository contains the computational workflow, analytical scripts, selected derived datasets, and manuscript figures and tables associated with a comparative genomic investigation of extracellular matrix formation and metal-related functions in Bacillus and related taxa.

## Study overview

The study evaluates the distribution, variability, and co-detection of genomic markers associated with biofilm regulation, extracellular polymeric substances, proteinaceous matrix components, poly-gamma-glutamate production, cell-wall polymers, and metal homeostasis.

A curated cohort of 638 genomes representing 123 species was analyzed using a functional reference panel comprising 53 marker labels.

The study further examines whether marker associations remain supported after accounting for species structure and assesses the availability of experimentally documented biofilm and metal-biosorption phenotypes.

## Main findings

- Functional marker distributions exhibited substantial species-associated structure.
- Of 53 evaluated markers, 41 varied across the genomic cohort and 12 were invariant.
- Global marker co-detection analyses identified 87 statistically significant associations after multiple-testing correction.
- No associations remained significant among the species-conditioned comparisons evaluated.
- Available phenotype–genome evidence was insufficient for a standardized validation of genomic predictors of metal biosorption.
- Genomic homology detection should not be interpreted as direct evidence of gene function, phenotype, or remediation performance.

Six reference targets (epsA, epsB, pgsE, capE, zntA, and zitB) could not be incorporated into the evaluated reference panel and must not be interpreted as biologically absent.

## Repository structure

- `workflow/scripts/`: analytical scripts and workflow modules.
- `Snakefile`: Snakemake workflow for selected analytical and manuscript-generation stages.
- `envs/`: environment specification and verified software versions.
- `config/`: project configuration.
- `resources/`: metadata, reference checksums, and supporting resources.
- `data/processed/`: selected derived datasets.
- `data/reference/`: functional reference panel.
- `results/functional/`: selected functional matrices and summaries.
- `results/statistics/`: selected statistical outputs.
- `results/figures/`: manuscript figures.
- `results/tables/`: manuscript tables.

Raw genome assemblies, downloaded full texts, intermediate datasets, and large alignment collections are retained locally and excluded from the lightweight public repository.

## Reproducibility

The project was developed using Python and command-line bioinformatics tools, including MMseqs2, skani, NCBI Datasets, and BLAST+.

Software versions are documented in `envs/verified_versions.tsv`, and a proposed environment specification is provided in `envs/environment.yml`.

The current `Snakefile` connects selected genomic validation, functional analysis, statistical analysis, and manuscript-generation stages. It is not yet a complete end-to-end implementation of every literature retrieval and genome processing module.

Existing outputs can be inspected without executing computational analyses using:

```bash
snakemake --cores 1 --dry-run
```

Full reproduction of stages requiring external resources depends on access to the corresponding source datasets and reference files.

## Manuscript outputs

| Manuscript item | Generating module |
|---|---|
| Figure 1 | Module 38 |
| Figure 2 | Module 39 |
| Figure 3 | Module 37 |
| Figure 4 | Module 40 |
| Table 1 | Module 41 |
| Table 2 | Module 42 |
| Table 3 | Module 44 |
| Table 4 | Module 43 |

The manuscript numbering is authoritative even where internal filenames retain historical numbering.

## Scientific limitations

Homology-based marker detection identifies candidate homologs rather than experimentally confirmed orthologs or functional phenotypes.

Global genomic associations may reflect species composition. The species-conditioned analysis did not provide statistically significant evidence supporting the evaluated associations independently of species structure.

Experimental phenotype observations were heterogeneous and insufficient for quantitative validation of genomic biosorption predictors.

The proposed relationship between bacterial extracellular matrices and plant-assisted metal remediation is a conceptual framework rather than a demonstration of remediation performance.

## License

Original project code, workflow configurations, and original documentation are distributed under the MIT License. See `LICENSE`.

Third-party software, sequence data, publications, and externally sourced reference materials remain subject to their respective licenses and terms of use. See `NOTICE.md`.

## Citation

Citation information will be provided in `CITATION.cff` and updated when the associated manuscript and Zenodo release receive permanent identifiers.

## Reproduction modes and data availability

Two distinct levels of reproducibility are supported by the
repository structure.

### 1. Reproducing manuscript outputs

The repository includes selected processed matrices,
statistical results, manuscript figures, and tables.

The corresponding analysis scripts are preserved under
`workflow/scripts/`.

The Snakemake workflow provides the dependency structure
for the principal analyses and manuscript outputs.

The workflow has been checked using:

```bash
snakemake --cores 1 --dry-run
```

This check confirms that the declared targets are currently
up to date in the original analysis environment. It does
not constitute independent end-to-end validation.

### 2. Reconstructing the complete genomic analysis

Full reconstruction additionally requires the original
genome assemblies, intermediate ANI results, and functional
reference sequences.

Large genome files and selected intermediate outputs are
excluded from the public repository.

The workflow includes modules for within-species ANI,
redundancy assessment, representative selection,
functional homolog searches, matrix construction,
statistical analyses, and manuscript outputs.

The complete workflow has not yet been validated on a
fresh installation with independently recovered inputs.

### Software environments

Environment documentation is available in:

- `envs/environment.yml`
- `envs/verified_versions.tsv`
- `envs/snakemake_versions.tsv`
- `envs/README.md`

The environment specification is provisional and should
not be interpreted as a validated installation lockfile.

### Configuration

Scientific parameters are documented in:

`config/config.yaml`

Some Snakefile rules still contain explicit parameters.
The configuration file is currently a reference record,
not a fully integrated runtime configuration.

### Data distribution

The repository distribution inventory is recorded in:

`resources/repository_data_manifest.tsv`

The inventory distinguishes files permitted by
`.gitignore` from files excluded from distribution.

Files permitted by `.gitignore` are not automatically
cleared for public redistribution.

### Functional reference sequences

Reference-panel provenance and limitations are described in:

`resources/references/README.md`

Third-party protein sequences are not automatically covered
by the MIT license applied to original project code.

### Scientific scope

Sequence similarity is interpreted as evidence of candidate
homology, not confirmed orthology or experimentally
demonstrated biological function.

The analyses do not establish validated genomic prediction
of metal biosorption phenotypes.
