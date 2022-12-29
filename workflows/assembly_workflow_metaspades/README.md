# Installation

The workflow can be installed with the following commands.  It is okay if the Checkm2 database install throws an error starting with `File "checkm2/bin/checkm2", line 244, in <module>...`.  The PhyloPhlAn run will intentionally fail in order to download the database but will take a while to download the database.
```
git clone https://github.com/WillNickols/assembly_workflow
cd assembly_workflow
conda env create -f assembly_environment.yml
conda activate biobakery_assembly
git clone --recursive https://github.com/chklovski/checkm2.git
checkm2/bin/checkm2 database --download --path databases/checkm/
export CHECKM_DATA_PATH=$(pwd)/databases/checkm/
mkdir tmp_to_delete && mkdir -p databases/phylophlan && touch tmp_to_delete/tmp.fa && phylophlan_metagenomic -d SGB.Jul20 --database_folder databases/phylophlan/ -i tmp_to_delete/; rm -r tmp_to_delete/
export PHYLOPHLAN_PATH=$(pwd)/databases/phylophlan/
```

Run the following commands to install the necessary R packages.
```
R
install.packages(c("docopt", "dplyr", "data.table", "stringr", "doParallel", "tidyr"))
q()
```

Once the conda environment is created and you are in the `assembly_workflow` directory, you can activate the environment with these commands:
```
conda activate biobakery_assembly
export CHECKM_DATA_PATH=$(pwd)/databases/checkm/
export PHYLOPHLAN_PATH=$(pwd)/databases/phylophlan/
export SPADES_BIN=/n/holystore01/LABS/huttenhower_lab/Users/wnickols/metaspades/SPAdes-3.15.4-Linux/bin/
```

# Example runs

```
python assembly_workflow.py \
  -i example/input/sample_1/ \
  -o example/output/sample_1/ \
  --abundance-type by_sample \
  --input-extension fastq.gz \
  --pair-identifier _R1 \
  --cores 8 \
  --local-jobs 12 \
  --time-per-gb 40 \
  --mem-per-gb 90
```

The output, `example/output/sample_1/final_profile_by_sample.tsv` is a MetaPhlAn-like output with the abundance of each taxonomic group identified in the sample.  The output shows that Pseudoalteromonas marina is present in the sample along with a species, `sgb_01` that represents a new species genome bin not within a Mash distance of 0.05 of any known SGB, 0.15 of any known GGB, or 0.3 of any known FGB.
