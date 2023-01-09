# Centrifuge

This folder contains code for the Centrifuge workflow.  The scripts `kreport2mpa.py` and `combine_mpa.py` were copied from Bracken into the `centrifuge` folder.

# Versioning:

- Centrifuge: 1.0.4
  - Database: Equivalent to the archaeal, bacterial, and viral database constructed on March 25, 2022.

# Installation

Centrifuge was installed according to [its manual instructions](https://github.com/DaehwanKimLab/centrifuge).  Because Centrifuge was added late in the evaluation process, its database was downloaded on January 5, 2023 but parsed to match what was available to Kraken 2 on March 25, 2022 (see the `database_building` folder).

# Running

Commands of the following form were run:

```
python centrifuge_workflow.py -i kneaddata_cleaned \
  -o centrifuge_output \
  --input-extension fastq.gz \
  --paired paired
```