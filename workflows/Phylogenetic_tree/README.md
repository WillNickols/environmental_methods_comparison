# Phylogenetic tree workflow

This folder contains code for building a phylogenetic tree from reference genomes.

# Installation

Run `install_tree_dependencies.sh` with the `biobakery_assembly` conda environment activated. This installs mashtree, quicktree, and ncbi-datasets-cli from source.

# Running

1. Collect TaxIDs (local, in RStudio): run `analysis/scripts/Phylogenetic_tree_generation/collect_all_taxids.R`

2. Build tree (cluster):

```
python workflows/Phylogenetic_tree/build_phylogenetic_tree.py \
  --taxid-list analysis/databases/Phylogenetic_trees/all_species_taxids.txt \
  -o workflows/Phylogenetic_tree/phylo_tree_output/ \
  --grid-partition test \
  --grid-jobs 1 \
  --cores 16
```
