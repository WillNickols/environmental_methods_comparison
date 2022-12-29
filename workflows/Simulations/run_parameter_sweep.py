# Takes in the profile without the abundances but sorted from most abundant to least abundant

from anadama2 import Workflow
import numpy as np
from numpy.random import default_rng
import os
import pandas as pd
import subprocess
from pathlib import Path
from datetime import date
import math

workflow = Workflow(version="0.1", description="Simulation workflow")
workflow.add_argument("config-file", desc="simulation configuration file")
workflow.add_argument("run-sim", desc="Whether to run simulators", action="store_true")
args = workflow.parse_args()

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

try:
	camisim_path = os.environ["CAMISIM_PATH"]
except:
	raise ValueError("CAMISIM_PATH not provided or set")

in_profile = args.input
output = args.output

with open(args.config_file) as f:
    lines = f.readlines()

configs = dict()
for line in lines:
    line = line.strip()
    splits = line.split(": ", 1)
    configs[splits[0]] = splits[1]

strain_n = 1

genomes_dir = output + "genomes_dir/"
os.makedirs(genomes_dir, exist_ok=True)

new_genomes_dir = genomes_dir + "new/"
os.makedirs(new_genomes_dir, exist_ok=True)

known_genomes_dir = genomes_dir + "known/"
os.makedirs(known_genomes_dir, exist_ok=True)

current_date = date.today().strftime("%Y%m%d")

def download_genomes():
    workflow.add_task("python " + this_folder + "download_new_genomes.py --out-folder " + new_genomes_dir + " --tax-dump-folder " + new_genomes_dir + "NCBI_tax/ --start-date " + configs['start_date'] + " --end-date " + configs['end_date'] + " --current-date " + current_date + " --n all --email " + os.environ["EMAIL"] + " --api-key " + os.environ["API_KEY"],
        targets=[new_genomes_dir + "NCBI_tax/ncbi-taxonomy_" + current_date + ".tar.gz", new_genomes_dir + "new_genomes_metadata.tsv"])

    command = "python " + this_folder + "gunzip_files.py -i " + new_genomes_dir + "genomes/ --local-jobs 1 -o " + new_genomes_dir + "genomes/ && touch " + new_genomes_dir + "genomes/gunzip.done"
    workflow.add_task(command, targets=new_genomes_dir + "genomes/gunzip.done", depends=[new_genomes_dir + "new_genomes_metadata.tsv"])

    profile = pd.read_csv(in_profile, sep='\t')
    profile[["UNKNOWN" not in item and "BIN_" not in item for item in profile['genome_ID']]].to_csv(known_genomes_dir + "general_reference_profile.tsv", index = False, sep = "\t")

    workflow.add_task("python " + this_folder + "download_known_genomes.py --reference-profile " + known_genomes_dir + "general_reference_profile.tsv --out-folder " + known_genomes_dir + " --email " + os.environ["EMAIL"] + " --api-key " + os.environ["API_KEY"],
        targets=known_genomes_dir + "general_reference_metadata.tsv", depends=known_genomes_dir + "general_reference_profile.tsv")

    command = "python " + this_folder + "gunzip_files.py -i " + known_genomes_dir + "genomes/ --local-jobs 1 -o " + known_genomes_dir + "genomes/ && touch " + known_genomes_dir + "genomes/gunzip.done"
    depends = [known_genomes_dir + "general_reference_metadata.tsv"]
    workflow.add_task(command, targets=known_genomes_dir + "genomes/gunzip.done", depends=depends)

download_genomes()

def make_profile(n, sample_size, k, sigma, unknown_prop, mutation_rate, sample_num):
    out_path = output + "original/" + "profile.n." + str(n) + ".size." + str(sample_size) + ".k." + str(k) + ".sigma." + str(sigma) + ".up." + str(unknown_prop) + ".mut_rate." + str(mutation_rate) + "_sample_" + str(sample_num) + "/"
    if not os.path.isfile(out_path + "out/profiles/" + "profile_" + str(sample_num) + "/taxonomic_profile_" + str(sample_num) +  ".txt"):
        np.random.seed(int(int(sample_num) + 3 * int(n) + 5 * float(sample_size) + 11 * int(k) + 13 * float(sigma) + 17 * float(unknown_prop) + 19 * float(mutation_rate)))
        profile = pd.read_csv(in_profile, sep='\t')
        if strain_n == 1:
            known_profile = profile[~profile['genome_ID'].str.contains("UNKNOWN")]
            known_profile = known_profile.drop_duplicates("genome_ID")
            unknown_profile = profile[profile['genome_ID'].str.contains("UNKNOWN")]
            profile = pd.concat([known_profile, unknown_profile]).reset_index(drop=True)

        num_unknown = len(profile[profile['genome_ID'].str.contains("UNKNOWN")]["genome_ID"])
        num_known = len(profile.index) - num_unknown

        # Split the profile into two known halves randomly and then reinterleave them for increased sample variability
        # This prevents the species from always being in the same abundance order while uusally keeping the more abundant
        # species towards the top of the list
        rng = default_rng()
        indexes = np.sort(rng.choice(num_known, size=int(num_known * 0.5), replace=False))
        split_1 = profile.iloc[indexes].reset_index(drop=True)
        split_2 = profile.drop(indexes).reset_index(drop=True)
        profile = pd.concat([split_1,split_2]).sort_index().reset_index(drop=True)

        profile = pd.concat([profile[:int(float((1 - float(unknown_prop)) * int(float(n))))], profile[(len(profile.index)-num_unknown):]]).reset_index(drop=True)
        abundances = np.sort(np.random.lognormal(-3, float(sigma), int(float((1 - float(unknown_prop)) * int(float(n))))))[::-1]
        abundances = (1 - float(unknown_prop)) * abundances/np.sum(abundances)
        abundances = np.concatenate((abundances, np.zeros(len(profile.index) - len(abundances))))
        profile.insert(loc=3, column='abundance', value=abundances)

        addition = pd.DataFrame({'genome_ID': ["UNKNOWN"], 'NCBI_ID': [""], 'filepath': [""], 'abundance': [unknown_prop], 'novelty_category': [""]})
        profile = pd.concat([profile, addition])

        if not os.path.isdir(out_path):
            os.makedirs(out_path)

        if not os.path.isdir(out_path + "profile/"):
            os.makedirs(out_path + "profile/")

        if not os.path.isdir(out_path + "out/"):
            os.makedirs(out_path + "out/")

        tmp = out_path + "tmp/"
        if not os.path.isdir(tmp):
            os.makedirs(tmp)

        inputs_dir = out_path + "out/inputs/"
        if not os.path.isdir(inputs_dir):
            os.makedirs(inputs_dir)

        profiles = out_path + "out/profiles/"
        if not os.path.isdir(profiles):
            os.makedirs(profiles)

        config_dir = out_path + "out/config/"
        if not os.path.isdir(config_dir):
            os.makedirs(config_dir)

        # Write profile and validate the inputs
        profile.to_csv(out_path + "profile/" + "profile.tsv", index = False, sep = "\t")
        if not profile.columns.tolist() == ["genome_ID", "NCBI_ID", "filepath", "abundance", "novelty_category"]:
            raise ValueError("Wrong column names")
        if not abs(sum(profile['abundance']) - 1) < 0.001:
            raise ValueError("Abundances don't sum to 1")
        if any([item == "" or item is None or pd.isnull(item) for item in profile.loc[["BIN_" in item for item in profile['genome_ID']], 'filepath'].tolist()]):
            raise ValueError("Assembled bins need file paths")
        if any([item == "" or item is None or pd.isnull(item) for item in profile.loc[["UNKNOWN" not in item for item in profile['genome_ID']], 'novelty_category'].tolist()]):
            raise ValueError("Non-UNKNOWN entries need a novelty category")
        if any([item == "" or item is None or pd.isnull(item) for item in profile['abundance'].tolist()]):
            raise ValueError("All entries need an abundance")
        if any([item == "" or item is None or pd.isnull(item) for item in profile.loc[["UNKNOWN" not in item for item in profile['genome_ID']], 'NCBI_ID'].tolist()]):
            raise ValueError("Non-UNKNOWN entries need an NCBI ID")

        workflow.add_task("cp " + new_genomes_dir + ". -r " + inputs_dir + ' && sed -i -- \'s/' + new_genomes_dir.replace("/", "\/") + '/' + inputs_dir.replace("/", "\/") + '/g\' ' + inputs_dir + "new_genomes_metadata.tsv",
            depends=[new_genomes_dir + "NCBI_tax/ncbi-taxonomy_" + current_date + ".tar.gz", new_genomes_dir + "new_genomes_metadata.tsv", new_genomes_dir + "genomes/gunzip.done"],
            targets=[inputs_dir + "NCBI_tax/ncbi-taxonomy_" + current_date + ".tar.gz", inputs_dir + "new_genomes_metadata.tsv"])

        profile[["UNKNOWN" not in item and "BIN_" not in item for item in profile['genome_ID']]].to_csv(inputs_dir + "reference_profile.tsv", index = False, sep = "\t")

        workflow.add_task("python " + this_folder + "move_known_genomes.py --general-reference-metadata " + known_genomes_dir + "general_reference_metadata.tsv" + " --reference-profile " + inputs_dir + "reference_profile.tsv --src " + known_genomes_dir + "genomes/" + " --dest " + inputs_dir + "genomes/",
            targets=inputs_dir + "reference_metadata.tsv", depends=[inputs_dir + "reference_profile.tsv", known_genomes_dir + "general_reference_metadata.tsv", known_genomes_dir + "genomes/gunzip.done"])

        files_to_move = profile[["BIN_" in item for item in profile['genome_ID']]]['filepath'].tolist()

        if not os.path.isdir(inputs_dir + "genomes/"):
            os.makedirs(inputs_dir + "genomes/")

        with open(inputs_dir + "files_to_copy.txt", 'w') as f:
            f.write("\n".join(files_to_move) + "\n")

        workflow.add_task("cat " + inputs_dir + "files_to_copy.txt | while read in; do cp \"$in\" " + inputs_dir + "genomes/ ; done && touch " + inputs_dir + "genomes/done",
                targets=inputs_dir + "genomes/done", depends=inputs_dir + "files_to_copy.txt")

        if mutation_rate > 0:
            workflow.add_task("python " + this_folder + "mutate.py --in-dir " + inputs_dir + "genomes/" + " --out-dir " + inputs_dir + "genomes/" + " --rate " + str(mutation_rate) + " --threads 1",
                targets=inputs_dir + "genomes/mutate.done", depends=[inputs_dir + "genomes/done", inputs_dir + "reference_metadata.tsv", inputs_dir + "new_genomes_metadata.tsv"])

        if not os.path.isfile(inputs_dir + "GC_and_lengths.tsv"):
            command = "python " + this_folder + "calculateGC.py --in-dir " + inputs_dir + "genomes/" + " --out-file " + inputs_dir + "GC_and_lengths.tsv" + " --threads 1"
            depends = [inputs_dir + "genomes/done", inputs_dir + "reference_metadata.tsv", inputs_dir + "new_genomes_metadata.tsv"]
            if mutation_rate > 0:
                depends.extend([inputs_dir + "genomes/mutate.done"])
            workflow.add_task(command, targets=inputs_dir + "GC_and_lengths.tsv", depends=depends)

        if not os.path.isdir(profiles + "profile_" + str(sample_num) + "/"):
            os.makedirs(profiles + "profile_" + str(sample_num) + "/")

        if not os.path.isdir(tmp + "tmp_" + str(sample_num) + "/"):
            os.makedirs(tmp + "tmp_" + str(sample_num) + "/")

        if not os.path.isdir(config_dir + "configs_" + str(sample_num) + "/"):
            os.makedirs(config_dir + "configs_" + str(sample_num) + "/")

        depends = [inputs_dir + "genomes/done", inputs_dir + "reference_profile.tsv", inputs_dir + "new_genomes_metadata.tsv", inputs_dir + "GC_and_lengths.tsv"]

        workflow.add_task("python " + this_folder + "make_config_files.py --profile " + out_path + "profile/" + "profile.tsv" + " --inputs-folder " + inputs_dir + " --configs-folder " + config_dir + "configs_" + str(sample_num) + "/ --strain-n " + str(strain_n) + " --threads 1 " +
            " --samples " + str(1) + " --tmp-dir " + tmp + "tmp_" + str(sample_num) + "/ --profiles-dir " + profiles + "profile_" + str(sample_num) + "/ --camisim-path " + os.environ["CAMISIM_PATH"] + " --sample-size " + str(sample_size) + " --current-date " + current_date + " --gc-len-file " + inputs_dir + "GC_and_lengths.tsv",
            targets=[config_dir + "configs_" + str(sample_num) + "/metadata.tsv", config_dir + "configs_" + str(sample_num) + "/genome_to_id.tsv", config_dir + "configs_" + str(sample_num) + "/config.ini"], depends=depends)

        workflow.add_task("Rscript " + this_folder + "get_genome_size_abundances.R --k " + str(k) + " --sigma " + str(sigma) + " --num_species " + str(n) + " --genome_to_id " + config_dir + "configs_" + str(sample_num) + "/genome_to_id.tsv" + " --profile " + out_path + "profile/" + "profile.tsv" + " --GC_and_lengths " + inputs_dir + "GC_and_lengths.tsv" + " -o " + config_dir + "abundances.tsv",
            targets=[config_dir + "abundances.tsv"], depends=[config_dir + "configs_" + str(sample_num) + "/genome_to_id.tsv"])

        workflow.add_task("cp " + config_dir + "abundances.tsv" + " " + config_dir + "configs_" + str(sample_num) + "/abundances.tsv",
            targets=[config_dir + "configs_" + str(sample_num) + "/abundances.tsv"], depends=[config_dir + "abundances.tsv"])

        if args.run_sim:
            if sample_num == 0:
                if not os.path.isfile(profiles + "profile_" + str(sample_num) + "/taxonomic_profile_" + str(sample_num) +  ".txt"):
                    workflow.add_task_gridable("python " + os.environ["CAMISIM_PATH"] + "metagenomesimulation.py " + config_dir + "configs_" + str(sample_num) + "/config.ini && for filename in " + profiles + "profile_" + str(sample_num) + "/*_sample*; do mv $filename " + profiles + "profile_" + str(sample_num) + "/sample_" + str(sample_num) + "; done && rm -r " + profiles + "profile_" + str(sample_num) + "/source_genomes/ && rm -r " + inputs_dir + "genomes/",
                        targets=[profiles + "profile_" + str(sample_num) + "/taxonomic_profile_" + str(sample_num) +  ".txt"],
                        depends=[config_dir + "configs_" + str(sample_num) + "/abundances.tsv", config_dir + "configs_" + str(sample_num) + "/metadata.tsv", config_dir + "configs_" + str(sample_num) + "/genome_to_id.tsv", config_dir + "configs_" + str(sample_num) + "/config.ini"],
                        time=120 + math.ceil(float(sample_size)) * 10,
                		mem=20000 + math.ceil(float(sample_size)) * 3000,
                		cores=int(configs['threads_per_simulator']),
                		partition=args.grid_partition,
                        interpret_deps_and_targs=False)
            else:
                if not os.path.isfile(profiles + "profile_" + str(sample_num) + "/taxonomic_profile_" + str(sample_num) +  ".txt"):
                    workflow.add_task_gridable("python " + os.environ["CAMISIM_PATH"] + "metagenomesimulation.py " + config_dir + "configs_" + str(sample_num) + "/config.ini && for filename in " + profiles + "profile_" + str(sample_num) + "/*_sample*; do mv $filename " + profiles + "profile_" + str(sample_num) + "/sample_" + str(sample_num) + "; done && mv " + profiles + "profile_" + str(sample_num) + "/taxonomic_profile_0.txt " + profiles + "profile_" + str(sample_num) + "/taxonomic_profile_" + str(sample_num) +  ".txt && rm -r " + profiles + "profile_" + str(sample_num) + "/source_genomes/ && rm -r " + inputs_dir + "genomes/",
                        targets=[profiles + "profile_" + str(sample_num) + "/taxonomic_profile_" + str(sample_num) +  ".txt"],
                        depends=[config_dir + "configs_" + str(sample_num) + "/abundances.tsv", config_dir + "configs_" + str(sample_num) + "/metadata.tsv", config_dir + "configs_" + str(sample_num) + "/genome_to_id.tsv", config_dir + "configs_" + str(sample_num) + "/config.ini"],
                        time=120 + math.ceil(float(sample_size)) * 10,
                		mem=20000 + math.ceil(float(sample_size)) * 3000,
                		cores=int(configs['threads_per_simulator']),
                		partition=args.grid_partition,
                        interpret_deps_and_targs=False)

list_to_run = [(int(n), float(configs['sample_sizes'].split(";")[1]), int(configs['ks'].split(";")[1]), float(configs['sigmas'].split(";")[1]), float(configs['unknown_props'].split(";")[1]), float(configs['mutation_rates'].split(";")[1])) for n in configs['ns'].split(";")[0].split(", ")]
list_to_run.extend([(int(configs['ns'].split(";")[1]), float(sample_size), int(configs['ks'].split(";")[1]), float(configs['sigmas'].split(";")[1]), float(configs['unknown_props'].split(";")[1]), float(configs['mutation_rates'].split(";")[1])) for sample_size in configs['sample_sizes'].split(";")[0].split(", ")])
list_to_run.extend([(int(configs['ns'].split(";")[1]), float(configs['sample_sizes'].split(";")[1]), int(k), float(configs['sigmas'].split(";")[1]), float(configs['unknown_props'].split(";")[1]), float(configs['mutation_rates'].split(";")[1])) for k in configs['ks'].split(";")[0].split(", ")])
list_to_run.extend([(int(configs['ns'].split(";")[1]), float(configs['sample_sizes'].split(";")[1]), int(configs['ks'].split(";")[1]), float(sigma), float(configs['unknown_props'].split(";")[1]), float(configs['mutation_rates'].split(";")[1])) for sigma in configs['sigmas'].split(";")[0].split(", ")])
list_to_run.extend([(int(configs['ns'].split(";")[1]), float(configs['sample_sizes'].split(";")[1]), int(configs['ks'].split(";")[1]), float(configs['sigmas'].split(";")[1]), float(unknown_prop), float(configs['mutation_rates'].split(";")[1])) for unknown_prop in configs['unknown_props'].split(";")[0].split(", ")])
list_to_run.extend([(int(configs['ns'].split(";")[1]), float(configs['sample_sizes'].split(";")[1]), int(configs['ks'].split(";")[1]), float(configs['sigmas'].split(";")[1]), float(configs['unknown_props'].split(";")[1]), float(mutation_rate)) for mutation_rate in configs['mutation_rates'].split(";")[0].split(", ")])
list_to_run = list(set(list_to_run))

for item in list_to_run:
    n, sample_size, k, sigma, unknown_prop, mutation_rate = item
    for i in range(int(configs['samples'])):
        make_profile(n, sample_size, k, sigma, unknown_prop, mutation_rate, i)

workflow.go()

#
