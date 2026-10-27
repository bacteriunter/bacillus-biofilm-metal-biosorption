# Functional reference panel

## Overview

The functional reference panel used in this study contains
99 protein reference records.

Files:

- data/reference/functional_panel/functional_reference_panel_merged.faa
- data/reference/functional_panel/functional_reference_panel_merged.tsv

All 99 records in the merged catalog identify
`reference_catalog.tsv` as their source catalog.

The catalog records protein accessions, gene symbols,
functional groups, reference organisms, sequence lengths,
and sequence-level SHA-256 checksums.

## Reference retrieval

Original reference files:

- reference_catalog.tsv
- reference_proteins.faa
- reference_retrieval.json

Supplementary retrieval records:

- reference_supplement_26b.faa
- reference_supplement_26b.tsv
- reference_supplement_26b_provenance.json

The supplementary retrieval log documents searches for
epsA, epsB, pgsB, pgsC, pgsE, capE, zntA, and zitB.

The log records accepted results for pgsB and pgsC.

The absence of an accepted reference for a target does
not establish biological absence of that gene.

## Interpretation

Reference-based sequence similarity identifies candidate
homologs. It does not independently establish orthology,
biochemical activity, or experimentally validated function.

## Reproducibility

Checksums for the merged reference files are recorded in:

resources/references/functional_panel_checksums.json

## Licensing and redistribution

The repository MIT license applies to original project
code and documentation.

Third-party protein sequences retain their respective
source attribution and applicable terms.

Redistribution permissions for the reference sequences
must be established before public release.

## UniProtKB attribution

The functional reference protein sequences and associated
source annotations were retrieved from UniProtKB using
the UniProt REST API.

Source: UniProt Consortium
Website: https://www.uniprot.org/
License: Creative Commons Attribution 4.0 International
License URL: https://creativecommons.org/licenses/by/4.0/
UniProt license information: https://www.uniprot.org/help/license

The reference panel is a project-specific selection and
reformatting of UniProtKB records.

Original protein accessions and source URLs are retained
in the reference catalog. Sequence SHA-256 checksums are
provided for provenance verification.

UniProt-derived material remains subject to its applicable
source license. The MIT license of this repository applies
to original project code and documentation and does not
replace the UniProtKB license.

UniProt does not endorse this project or its scientific
interpretations.
