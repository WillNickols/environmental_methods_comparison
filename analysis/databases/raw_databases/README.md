# Raw databases

This folder contains the raw database for each tool used.  Instructions to download the databases are as follows:

### MetaPhlAn 2
```
python get_database_MPA2.py -i mpa_v20_m200.pkl -o dbs/
```

### MetaPhlAn 3
```
python get_database_MPA3.py -i mpa_v30_CHOCOPhlAn_201901.pkl -o dbs/
```

### MetaPhlAn 4
```
python get_database_MPA4.py -i mpa_vJan21_CHOCOPhlAnSGB_202103.pkl -o dbs/
```

### mOTUs3
The mOTUs3 database was downloaded from [Zenodo](https://zenodo.org/record/5140350#.Y7tT5uzMKre).

### Metaxa 2
The Metaxa 2 database was downloaded from the `metaxa2_db/SSU/` folder of its installation.

### Kraken 2 / Bracken 2
The `seqid2taxid.map` file from Kraken was downloaded from its database folder.

### Centrifuge
The `seqid2taxid.map` file from Centrifuge was downloaded from its database folder.

### PhyloPhlAn 3
The `SGB.Jul20.txt` file was downloaded from the PhyloPhlAn 3 database folder.

### GTDB-Tk
The `gtdb_taxonomy.tsv` file was downloaded from the GTDB-Tk `release207_v2/taxonomy/` folder.