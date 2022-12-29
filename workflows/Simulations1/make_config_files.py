import os
import pandas as pd
import sys
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--profile", help="input profile")
parser.add_argument("--inputs-folder", help="input folder")
parser.add_argument("--configs-folder", help="configs folder")
parser.add_argument("--strain-n", help="strain number")
parser.add_argument("--local-jobs", help="number of local jobs", type=int)
parser.add_argument("--samples", help="number of samples to generate", type=int)
parser.add_argument("--tmp-dir", help="temporary directory to be deleted")
parser.add_argument("--profiles-dir", help="directory for the output profiles")
parser.add_argument("--camisim-path", help="directory for CAMISIM")
parser.add_argument("--sample-size", help="size of samples to generate", type=float)
parser.add_argument("--current-date", help="current date")
parser.add_argument("--gc-len-file", help="GC and genome length file")
parser.add_argument("--threads", help="GC and genome length file")
args = parser.parse_args()

profile = pd.read_csv(args.profile, sep='\t')
inputs_folder = args.inputs_folder
config_files = args.configs_folder
strain_n = args.strain_n
threads = args.threads
samples = args.samples
tmp = args.tmp_dir
profiles = args.profiles_dir
camisim_path = args.camisim_path
sample_size = args.sample_size
current_date = args.current_date

def to_int(in_list):
    return [int(item) for item in in_list]

# Create metadata.tsv
metadata = pd.DataFrame.from_dict({'genome_ID': profile[["BIN_" in item for item in profile['genome_ID']]]['genome_ID'].tolist(), 'OTU': to_int(profile[["BIN_" in item for item in profile['genome_ID']]]['NCBI_ID'].tolist()), 'NCBI_ID': to_int(profile[["BIN_" in item for item in profile['genome_ID']]]['NCBI_ID'].tolist()), 'novelty_category': profile[["BIN_" in item for item in profile['genome_ID']]]['novelty_category'].tolist()})
metadata = metadata.drop_duplicates('genome_ID')

reference_metadata = pd.read_csv(inputs_folder + "reference_metadata.tsv", sep='\t')

if str(strain_n) == str(1):
    reference_metadata = reference_metadata.drop_duplicates('genome_ID')
    genome_IDs = [gen  + "." + strain for gen, strain in zip(reference_metadata['genome_ID'].tolist(), reference_metadata['strain_ID'].tolist())]
    metadata = pd.concat([metadata,pd.DataFrame.from_dict({'genome_ID': genome_IDs, 'OTU': to_int(reference_metadata['NCBI_ID'].tolist()), 'NCBI_ID': to_int(reference_metadata['NCBI_ID'].tolist()), 'novelty_category': ["known_species"] * reference_metadata.shape[0]})], ignore_index = True)
elif str(strain_n) == "all":
    genome_IDs = [gen  + "." + strain for gen, strain in zip(reference_metadata['genome_ID'].tolist(), reference_metadata['strain_ID'].tolist())]
    metadata = pd.concat([metadata,pd.DataFrame.from_dict({'genome_ID': genome_IDs, 'OTU': to_int(reference_metadata['NCBI_ID'].tolist()), 'NCBI_ID': to_int(reference_metadata['NCBI_ID'].tolist()), 'novelty_category': ["known_species"] * reference_metadata.shape[0]})], ignore_index = True)
else:
    raise ValueError("strain-n must be 1 or all")

new_genomes_metadata = pd.read_csv(inputs_folder + "new_genomes_metadata.tsv", sep='\t')

if str(strain_n) == str(1):
    new_genomes_metadata = new_genomes_metadata.drop_duplicates('genome_ID')
    genome_IDs = ["UNKNOWN_" + gen  + "." + strain for gen, strain in zip(new_genomes_metadata['genome_ID'].tolist(), new_genomes_metadata['strain_ID'].tolist())]
    metadata = pd.concat([metadata,pd.DataFrame.from_dict({'genome_ID': genome_IDs, 'OTU': to_int(new_genomes_metadata['NCBI_ID'].tolist()), 'NCBI_ID': to_int(new_genomes_metadata['NCBI_ID'].tolist()), 'novelty_category': new_genomes_metadata['novelty_category'].tolist()})], ignore_index = True)
elif str(strain_n) == "all":
    genome_IDs = ["UNKNOWN_" + gen  + "." + strain for gen, strain in zip(new_genomes_metadata['genome_ID'].tolist(), new_genomes_metadata['strain_ID'].tolist())]
    metadata = pd.concat([metadata,pd.DataFrame.from_dict({'genome_ID': genome_IDs, 'OTU': to_int(new_genomes_metadata['NCBI_ID'].tolist()), 'NCBI_ID': to_int(new_genomes_metadata['NCBI_ID'].tolist()), 'novelty_category': new_genomes_metadata['novelty_category'].tolist()})], ignore_index = True)

# Create genome_to_id.tsv
genome_to_id = pd.DataFrame.from_dict({'genome_ID': profile[["BIN_" in item for item in profile['genome_ID']]]['genome_ID'].tolist(), 'filepath': [inputs_folder + "genomes/" + file.split("/")[-1] for file in profile[["BIN_" in item for item in profile['genome_ID']]]['filepath'].tolist()]})
genome_to_id = genome_to_id.drop_duplicates("genome_ID")

reference_metadata = pd.read_csv(inputs_folder + "reference_metadata.tsv", sep='\t')

if str(strain_n) == str(1):
    reference_metadata = reference_metadata.drop_duplicates('genome_ID')
    genome_IDs = [gen  + "." + strain for gen, strain in zip(reference_metadata['genome_ID'].tolist(), reference_metadata['strain_ID'].tolist())]
    genome_to_id = pd.concat([genome_to_id,pd.DataFrame.from_dict({'genome_ID': genome_IDs, 'filepath': reference_metadata['filepath']})], ignore_index = True)
elif str(strain_n) == "all":
    genome_IDs = [gen  + "." + strain for gen, strain in zip(reference_metadata['genome_ID'].tolist(), reference_metadata['strain_ID'].tolist())]
    genome_to_id = pd.concat([genome_to_id,pd.DataFrame.from_dict({'genome_ID': genome_IDs, 'filepath': reference_metadata['filepath']})], ignore_index = True)

new_genomes_metadata = pd.read_csv(inputs_folder + "new_genomes_metadata.tsv", sep='\t')

if str(strain_n) == str(1):
    new_genomes_metadata = new_genomes_metadata.drop_duplicates('genome_ID')
    genome_IDs = ["UNKNOWN_" + gen  + "." + strain for gen, strain in zip(new_genomes_metadata['genome_ID'].tolist(), new_genomes_metadata['strain_ID'].tolist())]
    genome_to_id = pd.concat([genome_to_id,pd.DataFrame.from_dict({'genome_ID': genome_IDs, 'filepath': new_genomes_metadata['filepath']})], ignore_index = True)
elif str(strain_n) == "all":
    genome_IDs = ["UNKNOWN_" + gen  + "." + strain for gen, strain in zip(new_genomes_metadata['genome_ID'].tolist(), new_genomes_metadata['strain_ID'].tolist())]
    genome_to_id = pd.concat([genome_to_id,pd.DataFrame.from_dict({'genome_ID': genome_IDs, 'filepath': new_genomes_metadata['filepath']})], ignore_index = True)

metadata = metadata.astype({"NCBI_ID": int, "OTU": int})
metadata.to_csv(config_files + "metadata.tsv", index = False, sep = "\t")
genome_to_id.to_csv(config_files + "genome_to_id.tsv", index = False, sep = "\t", header=False)

lines = ["[Main]",
"seed=" + str(config_files[-2]),
"phase=0",
"max_processors=" + str(threads),
"dataset_id=RL",
"output_directory=" + profiles,
"temp_directory=" + tmp,
"gsa=False",
"pooled_gsa=False",
"anonymous=True",
"compress=1",
"",
"[ReadSimulator]",
"readsim=" + camisim_path + "tools/art_illumina-2.3.6/art_illumina",
"error_profiles=" + camisim_path + "tools/art_illumina-2.3.6/profiles",
"samtools=" + camisim_path + "tools/samtools-1.3/samtools",
"profile=mbarc",
"type=art",
"fragments_size_mean=270",
"fragment_size_standard_deviation=27",
"",
"[CommunityDesign]",
"size=" + str(sample_size),
"distribution_file_paths=" + ",".join([config_files + "abundances.tsv"] * int(samples)),
"ncbi_taxdump=" + inputs_folder + "/NCBI_tax/ncbi-taxonomy_" + current_date + ".tar.gz",
"strain_simulation_template=" + camisim_path + "scripts/StrainSimulationWrapper/sgEvolver/simulation_dir",
"number_of_samples=" + str(samples),
"",
"[community0]",
"metadata=" + config_files + "metadata.tsv",
"id_to_genome_file=" + config_files + "genome_to_id.tsv",
"id_to_gff_file=",
"genomes_total=" + str(metadata.shape[0]),
"num_real_genomes=" + str(metadata.shape[0]),
"max_strains_per_otu=" + str(1 if (strain_n == 1) else 999),
"ratio=1",
"mode=differential",
"log_mu=1",
"log_sigma=2",
"gauss_mu=1",
"gauss_sigma=1",
"view=False"]

with open(config_files + "config.ini", 'w') as config:
    config.write('\n'.join(lines))
