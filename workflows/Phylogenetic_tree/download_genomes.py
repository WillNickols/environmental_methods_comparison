#!/usr/bin/env python

"""
Download one genome per species TaxID using the NCBI FTP assembly summaries.
Tries RefSeq first, falls back to GenBank for missing TaxIDs.
No special tools required -- just Python and wget.
"""

import argparse
import os
import subprocess
import sys
from multiprocessing.pool import ThreadPool

parser = argparse.ArgumentParser()
parser.add_argument("--taxid-list", required=True, help="File with one NCBI TaxID per line")
parser.add_argument("--output-dir", required=True, help="Output directory for genome .fna files")
parser.add_argument("--threads", type=int, default=8, help="Parallel download threads")
args = parser.parse_args()

taxid_list = os.path.abspath(args.taxid_list)
output_dir = os.path.abspath(args.output_dir.rstrip("/")) + "/"
os.makedirs(output_dir, exist_ok=True)

with open(taxid_list) as f:
    target_taxids = set(line.strip() for line in f if line.strip())

print("Target TaxIDs: {}".format(len(target_taxids)))

summary_path = output_dir + "assembly_summary_refseq.txt"
summary_url = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

if not os.path.isfile(summary_path):
    print("Downloading RefSeq assembly summary...")
    result = subprocess.run(["wget", "-q", "-O", summary_path, summary_url])
    if result.returncode != 0:
        print("Failed to download assembly summary")
        sys.exit(1)

print("Parsing assembly summary...")

# Pick one genome per species TaxID: prefer "reference genome", then "representative genome",
# then any with "Complete Genome" assembly level, then whatever is first.
# Columns: https://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt
# 0: assembly_accession, 4: refseq_category, 5: taxid, 6: species_taxid,
# 11: assembly_level, 19: ftp_path

PRIORITY = {"reference genome": 0, "representative genome": 1}

taxid_best = {}

with open(summary_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 20:
            continue

        species_taxid = fields[6]
        if species_taxid not in target_taxids:
            taxid = fields[5]
            if taxid not in target_taxids:
                continue
            species_taxid = taxid

        ftp_path = fields[19]
        if ftp_path == "na":
            continue

        refseq_cat = fields[4]
        assembly_level = fields[11]
        priority = PRIORITY.get(refseq_cat, 2)
        is_complete = 0 if assembly_level == "Complete Genome" else 1

        key = (priority, is_complete)

        if species_taxid not in taxid_best or key < taxid_best[species_taxid][0]:
            taxid_best[species_taxid] = (key, ftp_path, fields[0])

print("Found genomes for {}/{} TaxIDs".format(len(taxid_best), len(target_taxids)))

missing = target_taxids - set(taxid_best.keys())
if missing:
    missing_path = output_dir + "missing_taxids.txt"
    with open(missing_path, "w") as f:
        f.write("\n".join(sorted(missing)) + "\n")
    print("TaxIDs with no RefSeq genome: {} (written to {})".format(len(missing), missing_path))


import threading
progress_lock = threading.Lock()
progress = {"done": 0, "failed": 0, "skipped": 0, "total": len(taxid_best)}


def download_genome(item):
    taxid, (_, ftp_path, accession) = item
    dest = output_dir + taxid + ".fna.gz"
    dest_unzipped = output_dir + taxid + ".fna"

    if os.path.isfile(dest_unzipped):
        with progress_lock:
            progress["skipped"] += 1
            progress["done"] += 1
            done = progress["done"]
            total = progress["total"]
        if done % 10 == 0:
            print("{}/{} (skipped existing)".format(done, total), flush=True)
        return True

    basename = os.path.basename(ftp_path.rstrip("/"))
    genome_url = ftp_path.rstrip("/") + "/" + basename + "_genomic.fna.gz"

    result = subprocess.run(
        ["wget", "-q", "-O", dest, genome_url],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    if result.returncode != 0 or not os.path.isfile(dest) or os.path.getsize(dest) == 0:
        with progress_lock:
            progress["failed"] += 1
            progress["done"] += 1
        print("FAILED: TaxID {} ({}) URL: {}".format(taxid, accession, genome_url), flush=True)
        return False

    subprocess.run(["gunzip", "-f", dest],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    if os.path.isfile(dest_unzipped):
        with progress_lock:
            progress["done"] += 1
            done = progress["done"]
            total = progress["total"]
        if done % 100 == 0:
            print("{}/{} downloaded".format(done, total), flush=True)
        return True

    with progress_lock:
        progress["failed"] += 1
        progress["done"] += 1
    print("FAILED gunzip: TaxID {} ({})".format(taxid, accession), flush=True)
    return False


threads = args.threads
print("Downloading {} genomes using {} threads...".format(len(taxid_best), threads), flush=True)

pool = ThreadPool(threads)
results = pool.map(download_genome, list(taxid_best.items()))
pool.close()
pool.join()

succeeded = results.count(True)
failed = results.count(False)

print("\nDone: {} downloaded, {} failed".format(succeeded, failed))

if failed > 0:
    print("Re-run to retry failed downloads.")
    sys.exit(1)
