import pandas as pd
import io
import Bio
import re
from Bio import Entrez
import json
from os.path import exists
import argparse
import random
from ftplib import FTP
import os
import time

parser = argparse.ArgumentParser()
parser.add_argument("--known-species", help="list from makeProfiles.R")
parser.add_argument("--sgb-filepath", help="path to sgbs")
parser.add_argument("--sgb-info", help="SGB_INFO.tsv file")
parser.add_argument("--gtdbtk-relab", help="gtdbtk_relab.tsv from GTDBTK on SGBs")
parser.add_argument("--profile", help="output profile name")
args = parser.parse_args()

Entrez.email = os.environ["EMAIL"]
Entrez.api_key = os.environ["API_KEY"]
Entrez.max_tries = 20
Entrez.sleep_between_tries = 5

profile = pd.DataFrame(columns=['genome_ID', 'NCBI_ID', 'filepath', 'novelty_category'])

known_species = pd.read_csv(args.known_species, sep='\t', header=None)
known_species.columns = ['species']

ncbi_branches = ["bacteria", "archaea", "viral", "plasmid", "fungi", "protozoa", "invertebrate", "plant", "vertebrate_mammalian", "vertebrate_other"]

for _, row in known_species.iterrows():
    try:
        handle = Entrez.efetch(db = "taxonomy", id = str(row['species']), retmode="xml", rettype="gb")
    except:
        time.sleep(5)
        handle = Entrez.efetch(db = "taxonomy", id = str(row['species']), retmode="xml", rettype="gb")

    lines = handle.readlines()
    lines = [line.rstrip() for line in lines]
    str_match = [s for s in lines if "ScientificName" in s]
    if len(str_match) > 0:
        species_name = re.search('<ScientificName>(.*)</ScientificName>', str_match[0]).group(1)
    else:
        raise ValueError("Taxon not found in NCBI")

    accessible = False

    ftp = FTP('ftp.ncbi.nih.gov')
    ftp.login()

    print(species_name)

    branch_index = 0
    while branch_index < len(ncbi_branches):
        branch = ncbi_branches[branch_index]
        try:
            ftp.cwd('genomes/refseq/' + branch + '/' + species_name.replace(" ", "_") + '/')
            dir_list = ftp.nlst()
            ftp.close()
            if 'reference' in dir_list or 'representative' in dir_list:
                accessible = True
            break
        except:
            branch_index+=1

    if accessible:
        addition = pd.DataFrame({'genome_ID': [str(species_name)], 'NCBI_ID': [int(row['species'])], 'filepath': [''], 'novelty_category': ['known_species']})
        profile = pd.concat([profile, addition], ignore_index = True)

sgb_info = pd.read_csv(args.sgb_info, sep='\t')
gtdbtk_relab = pd.read_csv(args.gtdbtk_relab, sep='\t')

sgb_info['for_sort'] = - (sgb_info['completeness'] - 5 * sgb_info['contamination'])
sgb_info = sgb_info.sort_values(by=['for_sort'])

sgb_info = sgb_info.merge(gtdbtk_relab, 'left', left_on='sgb', right_on='mag')
sgb_list = sgb_info['sgb'].tolist()
taxon_list = sgb_info['taxon'].tolist()

taxon_abbrevs = {'k':'kingdom', 'p': 'phylum', 'c':'class', 'o':'order', 'f':'family', 'g':'genus', 's':'species'}

last_known_index = len(profile.index) - 1
i = 1
for sgb, taxon in zip(sgb_list, taxon_list):
    if taxon == 'UNKNOWN':
        taxID = 2
        genome_ID = 'BIN_UNKNOWN_' + str(i)
        novelty_category = 'known_kingdom'
        i += 1
        addition = pd.DataFrame({'genome_ID': [genome_ID], 'NCBI_ID': [int(taxID)], 'filepath': [args.sgb_filepath + sgb + ".fa"], 'novelty_category': [novelty_category]})
    else:
        while True:
            taxID = None
            novelty_category = ""

            if taxon.count('|') == 0:
                if "Archaea" in taxon:
                    taxID = 2157
                    genome_ID = 'BIN_UNKNOWN_' + str(i)
                    novelty_category = 'known_kingdom'
                    i += 1
                else:
                    taxID = 2
                    genome_ID = 'BIN_UNKNOWN_' + str(i)
                    novelty_category = 'known_kingdom'
                    i += 1
                addition = pd.DataFrame({'genome_ID': [genome_ID], 'NCBI_ID': [int(taxID)], 'filepath': [args.sgb_filepath + sgb + ".fa"], 'novelty_category': [novelty_category]})
                break

            pieces = taxon.rsplit("|", 1)
            taxon = pieces[0]
            full_name = pieces[1]
            just_name = full_name.split("__", 1)[1]

            try:
                handle = Entrez.esearch("taxonomy", just_name.replace("_", " "))
                record = Entrez.read(handle)
                handle.close()
                taxID = record["IdList"][0]

                try:
                    handle2 = Entrez.efetch(db = "taxonomy", id = taxID, retmode="xml", rettype="gb")
                except:
                    time.sleep(5)
                    handle2 = Entrez.efetch(db = "taxonomy", id = taxID, retmode="xml", rettype="gb")
                    
                lines = handle2.readlines()
                lines = [line.rstrip() for line in lines]
                str_match = [s for s in lines if "Rank" in s]
                return_level = re.search('<Rank>(.*)</Rank>', str_match[0]).group(1)
                if return_level == taxon_abbrevs[full_name[0]]:
                    novelty_category = 'known_' + taxon_abbrevs[full_name[0]]
                if novelty_category == "":
                    taxID = None
                    raise ValueError("No novelty category")
            except:
                pass

            if taxID is not None:
                if novelty_category == "known_species":
                    addition = pd.DataFrame({'genome_ID': ["BIN_" + just_name.replace("_", " ")], 'NCBI_ID': [int(taxID)], 'filepath': [args.sgb_filepath + sgb + ".fa"], 'novelty_category': [novelty_category]})
                else:
                    addition = pd.DataFrame({'genome_ID': ["BIN_UNKNOWN_" + str(i) + "_" + just_name.replace("_", " ")], 'NCBI_ID': [int(taxID)], 'filepath': [args.sgb_filepath + sgb + ".fa"], 'novelty_category': [novelty_category]})
                    i += 1
                break
    if int(addition['NCBI_ID'][0]) in profile['NCBI_ID'].tolist():
        indexes_of_interest = [i for i, x in enumerate(profile['NCBI_ID'].tolist()) if x == int(addition['NCBI_ID'][0])]

        # Keep up to 10 genera of new SGBs
        if len(indexes_of_interest) > 10 and addition['novelty_category'][0] == 'known_genus':
            indexes_of_interest_filepaths = [profile.loc[index_of_interest]['filepath'].split("/")[-1].split(".fa")[0] for index_of_interest in indexes_of_interest]

            scores = [sgb_info.loc[sgb_info['sgb'] == index_of_interest_filepath]['for_sort'].tolist()[0] for index_of_interest_filepath in indexes_of_interest_filepaths]
            worst = max(range(len(scores)), key=scores.__getitem__)
            if (sgb_info.loc[sgb_info['sgb'] == sgb]['for_sort']).tolist()[0] < scores[worst]:
                profile.loc[indexes_of_interest_filepaths[worst]] = addition.loc[0]
                print("SGB of same genus replaced")

        # Replace reference if SGB is > 90% complete and <5% contaminated
        elif len(indexes_of_interest) == 1 and profile.loc[indexes_of_interest[0]]['filepath'].split("/")[-1].split(".fa")[0] == "" and sgb_info.loc[sgb_info['sgb'] == sgb]['completeness'].tolist()[0] > 90 and sgb_info.loc[sgb_info['sgb'] == sgb]['contamination'].tolist()[0] < 5:
            profile.loc[indexes_of_interest[0]] = addition.loc[0]
            print("Replacing reference genome with SGB")

        # Otherwise add the SGB as is
        else:
            profile = pd.concat([profile, addition]).reset_index(drop=True)

    elif novelty_category == "known_species":
        index_to_insert = random.randint(0, last_known_index)
        profile = pd.concat([profile.iloc[:index_to_insert], addition, profile.iloc[index_to_insert:]]).reset_index(drop=True)
        last_known_index += 1
    else:
        profile = pd.concat([profile, addition]).reset_index(drop=True)


profile.to_csv(args.profile, index = False, sep = "\t")
