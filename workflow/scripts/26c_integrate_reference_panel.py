#!/usr/bin/env python3
import csv
import hashlib
import json
from pathlib import Path

root = Path("data/reference/functional_panel")

sources = [
    (root / "reference_catalog.tsv",
     root / "reference_proteins.faa"),
    (root / "reference_supplement_26b.tsv",
     root / "reference_supplement_26b.faa"),
]

def read_fasta(path):
    records = {}
    header = None
    seq = []

    def save():
        if header is not None:
            accession = header.split("|")[0].split()[0]
            sequence = "".join(seq).upper()
            if not sequence:
                raise ValueError(f"Secuencia vacía: {accession}")
            if accession in records:
                raise ValueError(f"Accesión duplicada en FASTA: {accession}")
            records[accession] = sequence

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        save()

    return records

combined = []
accessions = {}
sequences = {}
duplicates = []

for catalog_path, fasta_path in sources:
    fasta = read_fasta(fasta_path)

    with catalog_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    used = set()

    for row in rows:
        accession = (
            row.get("uniprot_accession")
            or row.get("protein_accession")
            or ""
        ).strip()

        if accession not in fasta:
            raise ValueError(f"FASTA faltante: {accession}")

        sequence = fasta[accession]
        digest = hashlib.sha256(sequence.encode()).hexdigest()

        expected = row.get("sequence_sha256", "")
        if expected and expected != digest:
            raise ValueError(f"SHA-256 inconsistente: {accession}")

        used.add(accession)

        if accession in accessions:
            if accessions[accession] != digest:
                raise ValueError(f"Accesión con secuencias distintas: {accession}")
            duplicates.append([accession, "accession"])
            continue

        if digest in sequences:
            duplicates.append([accession, "identical_sequence"])
            continue

        accessions[accession] = digest
        sequences[digest] = accession

        combined.append({
            "accession": accession,
            "gene": row.get("query_gene") or row.get("gene") or "",
            "group": row.get("group") or row.get("family") or "",
            "organism": row.get("organism") or row.get("reference_organism") or "",
            "source_catalog": catalog_path.name,
            "length_aa": len(sequence),
            "sequence_sha256": digest,
            "sequence": sequence,
        })

    extra = set(fasta) - used
    if extra:
        raise ValueError(f"FASTA sin catálogo en {fasta_path}: {sorted(extra)}")

combined.sort(key=lambda r: (r["group"], r["gene"], r["accession"]))

out_fasta = root / "functional_reference_panel_merged.faa"
out_catalog = root / "functional_reference_panel_merged.tsv"

with out_fasta.open("w") as handle:
    for row in combined:
        handle.write(f">{row['accession']}|{row['group']}|{row['gene']}\n")
        seq = row["sequence"]
        for i in range(0, len(seq), 80):
            handle.write(seq[i:i+80] + "\n")

columns = [
    "accession", "gene", "group", "organism",
    "source_catalog", "length_aa", "sequence_sha256"
]

with out_catalog.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
    writer.writeheader()
    for row in combined:
        writer.writerow({key: row[key] for key in columns})

summary = {
    "status": "PASS",
    "input_references": 131,
    "unique_references": len(combined),
    "duplicates_removed": len(duplicates),
    "duplicate_details": duplicates,
    "fasta": str(out_fasta),
    "catalog": str(out_catalog),
    "fasta_sha256": hashlib.sha256(out_fasta.read_bytes()).hexdigest(),
    "catalog_sha256": hashlib.sha256(out_catalog.read_bytes()).hexdigest(),
}

(root / "functional_reference_panel_merged_summary.json").write_text(
    json.dumps(summary, indent=2) + "\n"
)

print(json.dumps(summary, indent=2))
