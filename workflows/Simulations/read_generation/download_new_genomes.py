import sys
import pandas as pd
import Bio
import re
from Bio import Entrez
import urllib.request
import os
from os.path import exists
from ftplib import FTP
import tarfile
import shutil
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--out-folder", help="output folder")
parser.add_argument("--tax-dump-folder", help="folder for NCBI tax dump")
parser.add_argument("--start-date", help="start date in YYYY/MM/DD")
parser.add_argument("--end-date", help="end date in YYYY/MM/DD")
parser.add_argument("--current-date", help="current_date from datetime")
parser.add_argument("--n", help="number of species to download", default="all")
parser.add_argument("--email", help="email for Entrez")
parser.add_argument("--api-key", help="api key for Entrez")
args = parser.parse_args()

current_date = args.current_date

output = args.out_folder
output = os.path.abspath(output.rstrip("/")) + "/"
if not os.path.isdir(output):
    os.makedirs(output)
if not os.path.isdir(output + "genomes/"):
    os.makedirs(output + "genomes/")

taxdump = args.tax_dump_folder
taxdump = os.path.abspath(taxdump.rstrip("/")) + "/"
if not os.path.isdir(taxdump):
    os.makedirs(taxdump)

n = args.n

Entrez.email = args.email
Entrez.api_key = args.api_key
Bio.Entrez.max_tries = 20
Bio.Entrez.sleep_between_tries = 5

handle = Entrez.esearch("taxonomy", args.start_date + ":" + args.end_date + "[edat] AND Bacteria[SubTree] AND species[Rank] NOT unclassified[prop] NOT uncultured[prop] AND (\"above species level\"[prop] OR specified[prop])", retmax = 1000)
record = Entrez.read(handle)
handle.close()
taxIDs = record["IdList"]

urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz", taxdump + "ncbi-taxonomy_" + current_date + ".tar.gz")

if not os.path.isdir("/" + taxdump.strip("/") + "_tmp/NCBI/"):
    os.makedirs("/" + taxdump.strip("/") + "_tmp/NCBI/")
file = tarfile.open(taxdump + "ncbi-taxonomy_" + current_date + ".tar.gz")
file.extractall("/" + taxdump.strip("/") + "_tmp/NCBI/")
file.close()

os.remove(taxdump + "ncbi-taxonomy_" + current_date + ".tar.gz")

with tarfile.open(taxdump + "ncbi-taxonomy_" + current_date + ".tar.gz", "w:gz") as tar:
    tar.add("/" + taxdump.strip("/") + "_tmp/", arcname=os.path.basename("/" + taxdump.strip("/") + "_tmp/"))

shutil.rmtree("/" + taxdump.strip("/") + "_tmp/")

column_names = ["genome_ID", "NCBI_ID", "filepath", "novelty_category", "strain_ID"]
df = pd.DataFrame(columns = column_names)

if n == "all":
    n = len(taxIDs)

if len(taxIDs) < int(n):
    raise ValueError('Too many taxa requested')

for i in range(n):
    handle = Entrez.efetch(db = "taxonomy", id = str(taxIDs[i]), retmode="xml", rettype="gb")
    lines = handle.readlines()
    lines = [line.rstrip() for line in lines]
    str_match = [s for s in lines if "ScientificName" in s]
    if len(str_match) > 0:
        print(re.search('<ScientificName>(.*)</ScientificName>', str_match[0]).group(1))
        name = re.search('<ScientificName>(.*)</ScientificName>', str_match[0]).group(1)
    else:
        raise ValueError('Taxonomic name not present in NCBI lookup')

    ftp = FTP('ftp.ncbi.nih.gov')
    ftp.login()

    try:
        ftp.cwd('genomes/refseq/bacteria/' + name.replace(" ", "_") + '/latest_assembly_versions')
        strains = ftp.nlst()
    except:
        print("No strains found")
        strains = []

    for j in range(len(strains)):
        try:
            if not exists(output + "genomes/" + strains[j] + "_genomic.fna.gz"):
                print("Retrieving from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/" + name.replace(" ", "_") + "/latest_assembly_versions/" + strains[j] + "/" + strains[j] + "_genomic.fna.gz")
                urllib.request.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/" + name.replace(" ", "_") + "/latest_assembly_versions/" + strains[j] + "/" + strains[j] + "_genomic.fna.gz", output + "genomes/" + strains[j] + "_genomic.fna.gz")
        except:
            strains.remove(strains[j])
            j -= 1

    if len(strains) > 0:
        df = pd.concat([df,pd.DataFrame.from_dict({'genome_ID': ["NEW_" + name.replace(" ", "_")] * len(strains), 'NCBI_ID': [taxIDs[i]] * len(strains), 'filepath': [output + "genomes/" + strain + "_genomic.fa" for strain in strains], 'novelty_category': ["known_species"] * len(strains), 'strain_ID': strains})], ignore_index = True)

df.to_csv(output + "new_genomes_metadata.tsv", index = False, sep = "\t")

#https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=statistics&unclassified=hide&uncultured=hide&unspecified=hide
