# Kraken 2 / Bracken workflow

This folder contains code for the Kraken 2 / Bracken workflow.

# Versioning:

- Kraken 2: 2.1.2
  - Database: Default created from `kraken2-build --standard` on March 25, 2022
- Bracken: 2.6.2
  - Database: Created with default options in `bracken-build` on March 25, 2022

# Installation

Kraken 2 was installed according to [its manual instructions](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown) and the standard database was built on March 25, 2022 with `kraken2-build --standard`.  Bracken was installed with `install_bracken.sh` from [the Bracken GitHub](https://github.com/jenniferlu717/Bracken), and its database was built with default parameters on March 25, 2022.

# Running

Commands of the following form were run:

```
python kraken_workflow.py -i kneaddata_cleaned \
  -o kraken_output \
  --input-extension fastq.gz \
  --paired paired
```
