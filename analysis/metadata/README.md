# Metadata

The `forest_soil` folder contains two files related to metadata from Wilhelm et al. 2018: 
- `41396_2018_279_MOESM2_ESM.csv` contains a reformatted version of the Excel table in S1. 
- `41396_2018_279_MOESM4_ESM.csv` contains a reformatted version of the Excel table in S4.
- `filereport_read_run_PRJEB12502.tsv` contains the European Nucleotide Archive's summary data for project PRJEB12502.

The `human` folder contains the file downloaded with `https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv`.

The `sample_info` folder contains information about the samples used in the analysis.
- `read_info.tsv` contains read information on each sample individually.
- `Summary Read Information.pptx` contains summarized read information for each dataset.
- `biosamples_metadata.tsv` was created with `Rscript scripts/get_ena_file_reports.R -a sample_info/accessions.tsv`

Using the scripts in `scripts`, the `biosamples_metadata.tsv` file was split into a metadata file for each dataset. For the forest soil and human datasets, these were complemented with the metadata above.