# metaSPADES assembly workflow

This folder contains code for the metaSPADES-based assembly workflow, which is based on Aaron Walsh's pets assembly workflow and Pasolli et al. 2019.  The workflow consists of the following main steps:
- Create contigs using metaSPADES (for consistency with MEGAHIT, the contigs rather than the scaffolds were used as the contigs in subsequent steps)
- Build a Bowtie 2 index from the contigs and align reads with the flags `--very-sensitive-local` and `--no-unal`
- Calculate the contig read depth mean and variance with samtools' `view` and `sort` and `jgi_summarize_bam_contig_depths`
- Bin the contigs with MetaBAT 2 with a minimum contig length of 1500
- If using the `by_sample` approach, align each sample's reads to its bins with Bowtie 2 and calculate the relative abundance of each bin (accounting for genome size) with samtools' `view`, `sort`, and `index` and Checkm's `coverage` and `profile` scripts
- If using the `by_dataset` approach, align each sample's reads to all bins in the dataset with Bowtie 2 and calculate the relative abundance of each bin (accounting for genome size) with samtools' `view`, `sort`, and `index` and Checkm's `coverage` and `profile` scripts
- Calculate the n50 of each bin and the completeness and contamination of each bin using Checkm 2
- Place each bin using PhyloPhlAn 3
- If the bin has a Mash distance less than 0.05 from an SGB, less than 0.15 from a GGB, or less than 0.3 from a FGB, assign it to that bin
- Otherwise, recluster all remaining bins into SGBs based on a 0.05 Mash threshold

Because of the memory and time requirements of metaSPADES, some samples failed to assemble, so those were excluded from the merged output by running `assembly_workflow_cleanup.py` after `assembly_workflow.py`.

# Versioning

- SPADES: 3.15.4
- Bowtie2: 2.5.0
- Samtools: 1.16.1
- MetaBAT: 2.15
- Checkm 2: 1.0.0
- PhyloPhlAn 3: 3.0.67
  - Database: SGB.Jul20
- Mash: 2.3

# Installation

The requirements for the workflow can be installed by cloning this GitHub directory and following the commands below.  The PhyloPhlAn run will intentionally fail in order to download the database.
```
cd assembly_workflow_metaSPADES
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

Additionally, metaSPADES should be downloaded according to the instructions on [its GitHub](https://github.com/ablab/spades#sec2), and the bin should be saved as an environmental variable:
```
wget http://cab.spbu.ru/files/release3.15.4/SPAdes-3.15.4-Linux.tar.gz
tar -xzf SPAdes-3.15.4-Linux.tar.gz
export SPADES_BIN=$(pwd)/SPAdes-3.15.4-Linux/bin/
```

# Running

Commands of the following form were run:
```
python assembly_workflow.py \
  -i kneaddata_cleaned \
  -o assembly_metaSPADES \
  --abundance-type by_sample \
  --input-extension fastq.gz \
  --paired paired
```