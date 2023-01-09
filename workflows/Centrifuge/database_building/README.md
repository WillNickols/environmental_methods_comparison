Because Centrifuge was added later in the evaluation, the database was made to match the Kraken 2 database by removing taxa added after March 25, 2022.  

```
git clone https://github.com/infphilo/centrifuge
cd centrifuge
make

centrifuge-download -o taxonomy taxonomy
centrifuge-download -o library -m -d "archaea,bacteria,viral" refseq > seqid2taxid.map
python remove_new_genomes.py
# Rename _tmp files once correct
cat library/*/*.fna > input-sequences.fna
centrifuge-build -p 24 --conversion-table seqid2taxid.map \
                 --taxonomy-tree taxonomy/nodes.dmp --name-table taxonomy/names.dmp \
                 input-sequences.fna abv
```