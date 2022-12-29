# Simulation generation workflow

This folder contains code for the simulating metagenomic reads from environmental communitites.  The read generation steps are performed by a modified version of CAMISIM in which the gold standard assembly steps are removed.  These steps were computationally expensive, and the binning itself was not evaluated here.  Everything else was left as in the installation.  The workflow is split into two parts with the following main steps:

- Part 1: Profile building
  - Input MAGs were clustered into SGBs with `mash_clusters.R` and assigned a taxonomy by `gtdbtk_add_tax_assignment.py`
  - The most abundant NCBI-available species for each environment as determined by averaging the results of all the taxonomic profilers on the relevant real data sets was computed by `makeProfiles.R`.  This was sorted by average abundance and will be referred to as the "known profile."
  - All the NCBI-available species determined in the previous step were added to the profile so that their genomes would be downloaded from RefSeq.  Then, for each SGB with a GTDB-Tk assignment, the lowest-level taxonomy recognized by NCBI was determined with an Entrez lookup.  If the SGB did not have a recognized species, it was marked as an unknown species of the lowest taxonomic group for which it had an NCBI-recognized name.  Up to 10 SGBs from each genus were allowed (chosen based on a ranking of the completeness minus five times the contamination).  If the SGB had a recognized species assignment, it was added to the profile as that species.  If that species was already present in the known profile and the SGB was more than 90% complete with less than 5% contamination, the RefSeq genome was replaced with the SGB under the assumption that the SGB would better exemplify a strain of the species found in the environment.
- Part 2: Read generation
  - Genomes of species added to NCBI after the most recent database update of any method were downloaded as intentionally unknown species
  - RefSeq "representative" or "reference" genomes were downloaded for each species in the known profile that had not been replaced by a high-quality SGB
  - For the known-abundance proportion of the sample, the sorted known profile was randomly split into two groups, interleaved, and then assigned abundances from a normalized log normal distribution.  Because the original known profile was sorted by abundance in real samples, this splitting and interleaving allows species abundance rankings to change between samples while still ensuring that the most abundant species in real samples are still given high abundances in the simulation.
  - If the configuration specified a non-zero mutation rate, the genomes had point mutations induced at the specified rate
  - Genome lengths were calculated for each RefSeq genome and SGB
  - Abundance profiles were generated as follows:
    - The known profile abundances computed above were left as is
    - A list of desired genome sizes was generated to fit the distribution specified by the k parameter
    - The known genomes were placed at slots in this list closest to their genome sizes
    - The unknown species were chosen from new NCBI species and SGBs without an NCBI-recognized species assignment so that their genome sizes best fit the remaining slots.  The proportion of total species in the profile accounted for by these additions was equal to the unknown proportion from the configuration file.
    - These species were then assigned abundances from a log normal distribution with the same parameters as before but now scaled to the unknown species abundance proportion
  - Configuration files for CAMISIM were created
  - The modified version of CAMISIM was run

From here, reads were cleaned and input into the workflow pipelines for each tool.

# Versioning

- GTDB-Tk: 2.1.0
  - Database: R207_v2
- Mash: 2.3
- CAMISIM: 1.1.0 (Modified)

# Installation

GTDB-Tk was installed [based on GTDB-Tk's installation instructions](https://ecogenomics.github.io/GTDBTk/installing/bioconda.html) using Conda.  CAMISIM was installed from [its GitHub](https://github.com/CAMI-challenge/CAMISIM/wiki/User-manual) according to its installation instructions.  In `metagenomesimulation.py`, lines 101-111, 116, and 119-122 were removed to prevent the creation of gold standard assemblies.  (These steps were computationally costly, and the gold standards were not used.)

# Running

For the animal gut simulations, input SGBs were generated from MAGs from the dog (PRJEB34360), cat(PRJNA758898), and wild gut (PRJEB42019) samples in the MEGAHIT workflow.  For the ocean simulations, input SGBs were generated only from MAGs created from Tara Polar (PRJEB9740) samples in the MEGAHIT workflow.  Because the forest soil (PRJEB12502) samples did not generate enough MAGs for the number of SGBs needed in the soil simulation, samples from these projects were also used: PRJEB36534, PRJNA597671, PRJNA823845, PRJNA630300, PRJEB38557, PRJEB41174, PRJEB31111.  The MEGAHIT workflow was used on all soil samples.  Configuration files are available in the `configs` folder.

Commands of the following form were run after all MAGs were copied into `mags/`:
```
mkdir mash
ls mags/* > mash/input_list.txt
mkdir mash/main
mash sketch -l mash/input_list.txt -o mash/main/sketches
mash dist -t mash/main/sketches.msh mash/main/sketches.msh > mash/mash.tsv
mkdir checkm2
Rscript profile_building/make_input_sgbs/mash_clusters.R --mash mash/mash.tsv \
  --checkm checkm2/quality_report.tsv \
  --mag_dir mags \
  --out_dir sgbs \

gtdbtk classify_wf --genome_dir sgbs/sgbs/ --out_dir gtdbtk/ -x .fa
python profile_building/make_input_sgbs/gtdbtk_add_tax_assignment.py \
  --bac gtdbtk/gtdbtk.bac120.summary.tsv \
  --arc gtdbtk/gtdbtk.ar53.summary.tsv \
  --output gtdbtk/gtdbtk_relab.tsv

python profile_building/build_profile.py \
  --known-species known_species_list.tsv \
  --sgb-filepath sgbs/sgbs/ \
  --sgb-info sgbs/sgbs/SGB_info.tsv \
  --gtdbtk-relab gtdbtk/gtdbtk_relab.tsv \
  --profile profile.tsv

python read_generation/run_parameter_sweep.py \
  -i profile.tsv \
  -o simulations/ \
  --config-file configs/soil_config.txt \
  --run-sim
```
