#!/usr/bin/env python

"""
Flatten NCBI datasets download structure into a single directory with one genome
per species TaxID, named {taxid}.fna so mashtree tip labels match analysis scripts.
"""

import argparse
import json
import os
import shutil
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--download-dir", required=True, help="NCBI datasets download directory")
parser.add_argument("--output-dir", required=True, help="Output directory for flattened genomes")
args = parser.parse_args()

download_dir = os.path.abspath(args.download_dir.rstrip("/")) + "/"
output_dir = os.path.abspath(args.output_dir.rstrip("/")) + "/"
os.makedirs(output_dir, exist_ok=True)

# Collect all report files and data dirs across batch subdirectories
report_paths = []
data_dirs = []

# Check for batch subdirectories (batch_0/, batch_1/, ...)
for entry in sorted(Path(download_dir).iterdir()):
    if entry.is_dir() and entry.name.startswith("batch_"):
        report = entry / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
        data = entry / "ncbi_dataset" / "data"
        if report.is_file():
            report_paths.append(report)
        if data.is_dir():
            data_dirs.append(data)

# Also check for non-batched layout (single download)
single_report = Path(download_dir) / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
single_data = Path(download_dir) / "ncbi_dataset" / "data"
if single_report.is_file():
    report_paths.append(single_report)
if single_data.is_dir():
    data_dirs.append(single_data)

accession_to_taxid = {}
taxid_to_accession = {}

for report_path in report_paths:
    with open(str(report_path)) as f:
        for line in f:
            record = json.loads(line.strip())
            accession = record.get("accession", "")
            tax_id = str(record.get("organism", {}).get("taxId", ""))
            if accession and tax_id:
                if tax_id not in taxid_to_accession:
                    taxid_to_accession[tax_id] = accession
                    accession_to_taxid[accession] = tax_id

copied = 0
skipped_no_taxid = 0
skipped_duplicate = 0
seen_taxids = set()

for data_dir_path in data_dirs:
    for accession_dir in sorted(data_dir_path.iterdir()):
        if not accession_dir.is_dir():
            continue
        accession = accession_dir.name

        fna_files = list(accession_dir.glob("*_genomic.fna")) + list(accession_dir.glob("*.fna"))
        if not fna_files:
            continue

        fna_file = fna_files[0]

        tax_id = accession_to_taxid.get(accession, "")
        if not tax_id:
            skipped_no_taxid += 1
            continue

        if tax_id in seen_taxids:
            skipped_duplicate += 1
            continue

        seen_taxids.add(tax_id)
        dest = output_dir + tax_id + ".fna"
        shutil.copy2(str(fna_file), dest)
        copied += 1

print("Flattened genomes: {} copied, {} skipped (no TaxID), {} skipped (duplicate TaxID)".format(
    copied, skipped_no_taxid, skipped_duplicate))
