# Metaxa 2 workflow

This folder contains code for the Metaxa 2 workflow.

# Versioning:

- Metaxa 2: 2.2.3
  - Database: Default database with Metaxa version 2.2.3

# Installation

Metaxa 2 was installed according to [its manual instructions](https://microbiology.se/publ/metaxa2_users_guide_2.2.pdf).

# Running

Commands of the following form were run:

```
python metaxa2_workflow.py -i kneaddata_cleaned \
  -o metaxa2_output \
  --input-extension fastq.gz \
  --paired paired
```
