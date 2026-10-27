#!/usr/bin/env python3
"""bmetal: reproducible screening worksheet, not automated study selection.

All PubMed records are retained. No keyword-based exclusion or automatic
eligibility decision is made. Manual title/abstract and full-text decisions
must be recorded with reasons; genome matching is a separate stage.
"""
import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

FIELDS = ["pmid", "doi", "title", "year", "journal", "search_groups",
          "abstract_available", "taxon_level", "phenotype_domain",
          "title_abstract_decision", "title_abstract_reason",
          "full_text_decision", "full_text_reason", "strain_identifier",
          "experimental_measurement", "measurement_unit", "metal",
          "biosorbent_type", "genome_accession", "genome_match_verified",
          "reviewer", "review_date", "notes"]

PROTOCOL = """# bmetal — screening protocol v1.0 (prospective worksheet)

## Scope
Original experimental studies on Bacillus biofilm formation, extracellular
matrix or metal biosorption; plant–metal studies are retained for contextual
assessment. Species-only identifications may be retained as bibliographic
evidence. An identified isolate without a resolved species can also qualify.

## Two distinct evidence sets
- Bibliographic evidence: eligible experimental findings at species or isolate
  level, with phenotype, methods and quantitative result checked in full text.
- Genotype–phenotype set: only a phenotype measured in an identifiable isolate
  with a verified matching genome; species-level genome substitution prohibited.

## Title/abstract stage
Read every record independently of exploratory prioritization scores.
Use include_for_full_text, exclude, or uncertain. Record the explicit reason
for exclusion; absent information defaults to uncertain, not exclude.
Reviews may be excluded from the experimental evidence set but retained as
background references. Duplicate publications must be checked separately.

## Full-text stage
Record include, exclude, or unresolved. Verify original experimental work,
organism identity, phenotype and measurement, biosorbent type (intact cells,
biofilm, isolated EPS, dead biomass, etc.), assay conditions, and whether
results are extractable. For exclusion record a specific verifiable reason.
Do not equate metal tolerance, metal reduction/precipitation or plant growth
with adsorption/biosorption; classify them as distinct outcomes.

## Methodological safeguards
Automated keywords may assist navigation only, not decide inclusion.
Maintain a full-text retrieval log and a decision log. Record reviewer identity,
date, and disagreements; ideally use independent double screening with
adjudication. PubMed alone is not an exhaustive environmental microbiology
search; document additional databases and exact search dates before claiming
comprehensive coverage. Never infer a genome match from species alone.

## Status
Protocol for human screening. No final inclusion/exclusion decisions have
been made by this script.
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="data/raw/pubmed_records.tsv")
    ap.add_argument("--output", default="data/interim/pubmed_eligibility_review.tsv")
    ap.add_argument("--summary", default="results/reports/pubmed_eligibility_setup.json")
    ap.add_argument("--protocol", default="docs/screening_protocol.md")
    args = ap.parse_args()
    source = Path(args.input)
    if not source.is_file():
        ap.error(f"Input not found: {source}")
    with source.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"pmid", "doi", "title", "year", "journal", "search_groups", "abstract"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            ap.error(f"Missing input columns: {sorted(missing)}")
        records = list(reader)
    ids = [r["pmid"].strip() for r in records]
    if any(not p for p in ids) or len(ids) != len(set(ids)):
        ap.error("Missing or repeated PMID; audit input before proceeding")
    output = Path(args.output)
    if output.exists():
        ap.error(f"Output already exists: {output}. Preserve human decisions; choose another --output.")
    protocol = Path(args.protocol)
    if protocol.exists():
        ap.error(f"Protocol already exists: {protocol}; do not overwrite it.")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        for r in records:
            row = {k: r.get(k, "") for k in FIELDS}
            row["abstract_available"] = "yes" if r["abstract"].strip() else "no"
            row["title_abstract_decision"] = "pending_manual_review"
            row["full_text_decision"] = "not_assessed"
            writer.writerow(row)
    protocol.parent.mkdir(parents=True, exist_ok=True)
    protocol.write_text(PROTOCOL, encoding="utf-8")
    summary = {"generated_utc": datetime.now(timezone.utc).isoformat(),
               "source": str(source), "records": len(records),
               "pending_manual_review": len(records),
               "no_abstract": sum(not r["abstract"].strip() for r in records),
               "missing_doi": sum(not r["doi"].strip() for r in records),
               "automatic_exclusions": 0,
               "note": "No eligibility decisions made automatically."}
    report = Path(args.summary)
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Records retained: {len(records)}; automatic exclusions: 0")
    print(f"Review worksheet: {output}")
    print(f"Protocol: {protocol}")
    print(f"Report: {report}")

if __name__ == "__main__":
    main()
