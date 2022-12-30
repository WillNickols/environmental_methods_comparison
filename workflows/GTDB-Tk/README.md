# GTDB-Tk workflow

This folder contains code for the GTDB-Tk workflow, which is designed to be run after the `assembly_workflow_MEGAHIT` or `assembly_workflow_metaSPADES` workflow.

# Versioning

- GTDB-Tk: 2.1.0
  - Database: R207_v2

# Installation

GTDB-Tk was installed [based on GTDB-Tk's installation instructions](https://ecogenomics.github.io/GTDBTk/installing/bioconda.html) using Conda.

# Running

Commands of the following form were run:

```
python gtdbtk_workflow.py -i assembly_MEGAHIT/ \
  -o gtdbtk/ \
  --abundance-type by_sample
```