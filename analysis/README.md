# Analysis

This folder contains the workflow outputs, the scripts used to analyze them, and the resulting figures.  The `databases`, `figures`, and `metadata` folders contain their own READMEs.

- The `databases` folder contains scripts and data to determine what taxa each tool can assign to a sample.  Its `figures` folder contains upset plots showing the overlap in these databases.
- The `figures` folder contains the outputs of analysis scripts for both the real and simulated data.
- The `metadata` folder contains information about the reads from each dataset (individual and aggregated) and metadata from the North American forest soil dataset used for running PERMANOVAs.
- The `scripts` folder contains scripts to find the read numbers of the samples in each dataset (`get_read_number.R`), calculate the genome lengths for genomes used in the simulations (`simulation_genome_lengths.py`), analyze the real data and create figures (`RealAnalysis.R` and `helpers.R`), and analyze the simulated data and create figures (`SimulationAnalysis.R` and `helpers.R`).
- The `real_data_outputs` folder contains the results of the taxonomic assignment workflows for each dataset.
- The `simulation_outputs` folder contains the results of the taxonomic assignment workflows for each dataset.