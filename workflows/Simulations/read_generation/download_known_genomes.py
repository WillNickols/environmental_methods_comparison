import sys
import pandas as pd
import Bio
import re
from Bio import Entrez
import urllib.request
import os
from os.path import exists
from ftplib import FTP
import gzip
import shutil
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--reference-profile", help="reference profile tsv")
parser.add_argument("--out-folder", help="output folder")
parser.add_argument("--email", help="email for Entrez")
parser.add_argument("--api-key", help="api key for Entrez")
args = parser.parse_args()

input = args.reference_profile
output = args.out_folder

output = os.path.abspath(output.rstrip("/")) + "/"
if not os.path.isdir(output):
    os.makedirs(output)
if not os.path.isdir(output + "genomes/"):
    os.makedirs(output + "genomes/")

df = pd.read_csv(input, sep = '\t')

Entrez.email = args.email
Entrez.api_key = args.api_key
Bio.Entrez.max_tries = 20
Bio.Entrez.sleep_between_tries = 5

column_names = ["genome_ID", "NCBI_ID", "filepath", "novelty_category",  "strain_ID"]
output_df = pd.DataFrame(columns = column_names)

ncbi_branches = ["bacteria", "archaea", "viral", "plasmid", "fungi", "protozoa", "invertebrate", "plant", "vertebrate_mammalian", "vertebrate_other"]

names = []
for tax_id in df['NCBI_ID']:
    try:
        handle = Entrez.efetch(db = "taxonomy", id = str(tax_id), retmode="xml", rettype="gb")
    except:
        handle = Entrez.efetch(db = "taxonomy", id = str(tax_id), retmode="xml", rettype="gb")
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

    branch_index = 0
    while True:
        if branch_index == len(ncbi_branches):
            raise ValueError("Genome not found")
        branch = ncbi_branches[branch_index]
        try:
            ftp.cwd('genomes/refseq/' + branch + '/' + name.replace(" ", "_") + '/')
            dir_list = ftp.nlst()
            ftp.close()
            if 'reference' in dir_list:
                ftp = FTP('ftp.ncbi.nih.gov')
                ftp.login()
                ftp.cwd('genomes/refseq/' + branch + '/' + name.replace(" ", "_") + '/reference/')
                strains = ftp.nlst()
                ftp.close()
                genome_type = 'reference'
            elif 'representative' in dir_list:
                ftp = FTP('ftp.ncbi.nih.gov')
                ftp.login()
                ftp.cwd('genomes/refseq/' + branch + '/' + name.replace(" ", "_") + '/representative/')
                strains = ftp.nlst()
                ftp.close()
                genome_type = 'representative'
            break
        except:
            branch_index+=1

        print("No genome found")
        strains = []


    for j in range(len(strains)):
        try:
            if not exists(output + "genomes/" + strains[j] + "_genomic.fna"):
                print("Retrieving from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/" + branch + "/" + name.replace(" ", "_") + "/" + genome_type + "/" + strains[j] + "/" + strains[j] + "_genomic.fna")
                urllib.request.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/" + branch + "/" + name.replace(" ", "_") + "/" + genome_type + "/" + strains[j] + "/" + strains[j] + "_genomic.fna.gz", output + "genomes/" + strains[j] + "_genomic.fna.gz")
        except:
            print("Download failed on " + strains[j])
            strains.remove(strains[j])
            j -= 1

    output_df = pd.concat([output_df,pd.DataFrame.from_dict({'genome_ID': [name.replace(" ", "_")] * len(strains), 'NCBI_ID': [tax_id] * len(strains), 'filepath': [output + "genomes/" + strain + "_genomic.fa" for strain in strains], 'novelty_category': ["known_species"] * len(strains), 'strain_ID': strains})], ignore_index = True)

output_df.to_csv(output + "general_reference_metadata.tsv", index = False, sep = "\t")

# /n/holystore01/LABS/huttenhower_lab/Users/wnickols/CAMISIM/tools/assembly_summary_complete_genomes.txt
