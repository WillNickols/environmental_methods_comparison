
```
/n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/build/bin/simka \
  -in /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/will_test/in.txt \
  -out /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/will_results/ \
  -out-tmp /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/will_tmp/ \
  -max-memory 20000 \
  -nb-cores 8

```

```
hutlab load centos7/python3/biobakery_workflows/3.0.0-beta-devel-dependsUpdate

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simka_workflow.py \
  --input /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/coastal_new/kneaddata/ \
  --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/coastal_new/simka/ \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/coastal/simka/ \
  --grid-partition 'shared' \
  --grid-jobs 96 \
  --cores 8 \
  --mem 30000 \
  --paired paired \
  --simka-path /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/build/bin/simka \
  --local-jobs 16

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simka_workflow.py \
  --input /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mine_new/kneaddata/ \
  --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mine_new/simka/ \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/mine/simka/ \
  --grid-partition 'shared' \
  --grid-jobs 96 \
  --cores 8 \
  --mem 30000 \
  --paired paired \
  --simka-path /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/build/bin/simka \
  --local-jobs 16

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simka_workflow.py \
  --input /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/gators_new/kneaddata/ \
  --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/gators_new/simka/ \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/gators/simka/ \
  --grid-partition 'shared' \
  --grid-jobs 96 \
  --cores 8 \
  --mem 30000 \
  --paired unpaired \
  --simka-path /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/build/bin/simka \
  --local-jobs 16

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simka_workflow.py \
  --input /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/delignification_new/kneaddata/ \
  --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/delignification_new/simka/ \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/delignification/simka/ \
  --grid-partition 'shared' \
  --grid-jobs 96 \
  --cores 8 \
  --mem 60000 \
  --paired paired \
  --simka-path /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/build/bin/simka \
  --local-jobs 16

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simka_workflow.py \
  --input /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/cats_out/kneaddata/ \
  --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/cats_out/simka/ \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/cats_out/simka/ \
  --grid-partition 'shared' \
  --grid-jobs 96 \
  --cores 8 \
  --mem 30000 \
  --paired paired \
  --simka-path /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/build/bin/simka \
  --local-jobs 16

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simka_workflow.py \
  --input /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/dogs_out/kneaddata/ \
  --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/dogs_out/simka/ \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/dogs_out/simka/ \
  --grid-partition 'shared' \
  --grid-jobs 96 \
  --cores 8 \
  --mem 30000 \
  --paired paired \
  --simka-path /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/build/bin/simka \
  --local-jobs 16

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simka_workflow.py \
  --input /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/human/kneaddata/ \
  --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/human/simka/ \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/human/simka/ \
  --grid-partition 'shared' \
  --grid-jobs 96 \
  --cores 8 \
  --mem 30000 \
  --paired paired \
  --simka-path /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/build/bin/simka \
  --local-jobs 16 \
  --no-unmatched y

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simka_workflow.py \
  --input /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/saltmarsh/kneaddata/ \
  --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/saltmarsh/simka/ \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/saltmarsh/simka/ \
  --grid-partition 'shared' \
  --grid-jobs 96 \
  --cores 8 \
  --mem 30000 \
  --paired paired \
  --simka-path /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/build/bin/simka \
  --local-jobs 16

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simka_workflow.py \
  --input /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/wild_gut/kneaddata2/ \
  --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/wild_gut/simka/ \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/wild_gut/simka/ \
  --grid-partition 'shared' \
  --grid-jobs 96 \
  --cores 8 \
  --mem 30000 \
  --paired paired \
  --simka-path /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/build/bin/simka \
  --local-jobs 16

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simka_workflow.py \
  --input /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/tara_polar_new/kneaddata/ \
  --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/tara_polar_new/simka/ \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/tara_polar_new/simka/ \
  --grid-partition 'shared' \
  --grid-jobs 96 \
  --cores 8 \
  --mem 30000 \
  --paired paired \
  --simka-path /n/holystore01/LABS/huttenhower_lab/Users/wnickols/simka/simka/build/bin/simka \
  --local-jobs 16
```

Re-running with full dataset comparison
```
conda activate biobakery_assembly
export CHECKM_DATA_PATH=$(pwd)/databases/checkm/
export PHYLOPHLAN_PATH=$(pwd)/databases/phylophlan/

python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/gators_new/kneaddata/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/gators_new/assembly/ \
  --abundance-type by_dataset --input-extension fastq.gz --paired unpaired \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/gators/assembly/ \
  --grid-partition 'shared' --grid-jobs 96 --cores 8 --time 10000 --mem 40000 \
  --local-jobs 12 \
  --skip-contigs y \
  --remove-intermediate-files y \
  --grid-options="--account=nguyen_lab"

python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mine_new/kneaddata/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mine_new/assembly/ \
  --abundance-type by_dataset --input-extension fastq.gz --paired paired \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/mine/assembly/ \
  --grid-partition 'shared' --grid-jobs 96 --cores 8 --time 10000 --mem 40000 \
  --local-jobs 12 \
  --skip-contigs y \
  --remove-intermediate-files y \
  --grid-options="--account=nguyen_lab"

python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/delignification_new/kneaddata/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/delignification_new/assembly/ \
  --abundance-type by_dataset --input-extension fastq.gz --paired paired \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/delignification/assembly/ \
  --grid-partition 'shared' --grid-jobs 96 --cores 8 --time 10000 --mem 40000 \
  --local-jobs 12 \
  --skip-contigs y \
  --remove-intermediate-files y \
  --grid-options="--account=nguyen_lab"

python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/tara_polar_new/kneaddata/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/tara_polar_new/assembly/ \
  --abundance-type by_dataset --input-extension fastq.gz --paired paired \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/tara_polar/assembly/ \
  --grid-partition 'shared' --grid-jobs 96 --cores 8 --time 10000 --mem 40000 \
  --local-jobs 12 \
  --skip-contigs y \
  --remove-intermediate-files y \
  --grid-options="--account=nguyen_lab"

python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/coastal_new/kneaddata/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/coastal_new/assembly/ \
  --abundance-type by_dataset --input-extension fastq.gz --paired paired \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/coastal/assembly/ \
  --grid-partition 'shared' --grid-jobs 96 --cores 8 --time 10000 --mem 40000 \
  --local-jobs 12 \
  --skip-contigs y \
  --remove-intermediate-files y \
  --grid-options="--account=nguyen_lab" \
  --dry-run

python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/cats_out/kneaddata/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/cats_out/assembly/ \
  --abundance-type by_dataset --input-extension fastq.gz --paired paired \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/cats_out/assembly/ \
  --grid-partition 'shared' --grid-jobs 96 --cores 8 --time 10000 --mem 40000 \
  --local-jobs 12 \
  --skip-contigs y \
  --remove-intermediate-files y \
  --grid-options="--account=nguyen_lab" \
  --dry-run

python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/dogs_out/kneaddata/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/dogs_out/assembly/ \
  --abundance-type by_dataset --input-extension fastq.gz --paired paired \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/dogs_out/assembly/ \
  --grid-partition 'shared' --grid-jobs 96 --cores 8 --time 10000 --mem 40000 \
  --local-jobs 12 \
  --skip-contigs y \
  --remove-intermediate-files y \
  --grid-options="--account=nguyen_lab" \
  --dry-run

python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/human/kneaddata/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/human/assembly/ \
  --abundance-type by_dataset --input-extension fastq.gz --paired paired \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/human/assembly/ \
  --grid-partition 'shared' --grid-jobs 96 --cores 8 --time 10000 --mem 40000 \
  --local-jobs 12 \
  --skip-contigs y \
  --remove-intermediate-files y \
  --grid-options="--account=nguyen_lab" \
  --dry-run

python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/saltmarsh/kneaddata/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/saltmarsh/assembly/ \
  --abundance-type by_dataset --input-extension fastq.gz --paired paired \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/saltmarsh/assembly/ \
  --grid-partition 'shared' --grid-jobs 96 --cores 8 --time 10000 --mem 40000 \
  --local-jobs 12 \
  --skip-contigs y \
  --remove-intermediate-files y \
  --grid-options="--account=nguyen_lab" \
  --dry-run

python assembly_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/wild_gut/kneaddata2/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/wild_gut/assembly/ \
  --abundance-type by_dataset --input-extension fastq.gz --paired paired \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/wild_gut/assembly/ \
  --grid-partition 'shared' --grid-jobs 96 --cores 8 --time 10000 --mem 40000 \
  --local-jobs 12 \
  --skip-contigs y \
  --remove-intermediate-files y \
  --grid-options="--account=nguyen_lab" \
  --dry-run
```


```
source /n/holystore01/LABS/huttenhower_lab/Users/wnickols/gtdbtk/env/bin/activate
hutlab load centos7/python3/biobakery_workflows/3.0.0-beta-devel
export GTDBTK_DATA_PATH=/n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/gtdbtk/release207_v2/
conda activate gtdbtk

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/gtdbtk_workflow/gtdbtk_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/coastal_new/assembly/bins/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/coastal_new/gtdbtk/ \
  --abundance-type by_dataset \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/coastal/gtdbtk/ \
  --cores 8 \
  --mags-sgb-dir /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/coastal_new/assembly/ \
  --local-jobs 10

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/gtdbtk_workflow/gtdbtk_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/gators_new/assembly/bins/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/gators_new/gtdbtk/ \
  --abundance-type by_dataset \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/gators/gtdbtk/ \
  --cores 8 \
  --mags-sgb-dir /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/gators_new/assembly/ \
  --local-jobs 10

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/gtdbtk_workflow/gtdbtk_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mine_new/assembly/bins/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mine_new/gtdbtk/ \
  --abundance-type by_dataset \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/mine/gtdbtk/ \
  --cores 8 \
  --mags-sgb-dir /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/mine_new/assembly/ \
  --local-jobs 10

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/gtdbtk_workflow/gtdbtk_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/delignification_new/assembly/bins/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/delignification_new/gtdbtk/ \
  --abundance-type by_dataset \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/delignification/gtdbtk/ \
  --cores 8 \
  --mags-sgb-dir /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/delignification_new/assembly/ \
  --local-jobs 10

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/gtdbtk_workflow/gtdbtk_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/tara_polar_new/assembly/bins/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/tara_polar_new/gtdbtk/ \
  --abundance-type by_dataset \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/tara_polar/gtdbtk/ \
  --cores 8 \
  --mags-sgb-dir /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/tara_polar_new/assembly/ \
  --local-jobs 10

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/gtdbtk_workflow/gtdbtk_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/cats_out/assembly/bins/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/cats_out/gtdbtk/ \
  --abundance-type by_dataset \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/cats_out/gtdbtk/ \
  --cores 8 \
  --mags-sgb-dir /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/cats_out/assembly/ \
  --local-jobs 10

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/gtdbtk_workflow/gtdbtk_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/dogs_out/assembly/bins/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/dogs_out/gtdbtk/ \
  --abundance-type by_dataset \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/dogs_out/gtdbtk/ \
  --cores 8 \
  --mags-sgb-dir /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/dogs_out/assembly/ \
  --local-jobs 10

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/gtdbtk_workflow/gtdbtk_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/human/assembly/bins/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/human/gtdbtk/ \
  --abundance-type by_dataset \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/human/gtdbtk/ \
  --cores 8 \
  --mags-sgb-dir /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/human/assembly/ \
  --local-jobs 10

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/gtdbtk_workflow/gtdbtk_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/saltmarsh/assembly/bins/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/saltmarsh/gtdbtk/ \
  --abundance-type by_dataset \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/saltmarsh/gtdbtk/ \
  --cores 8 \
  --mags-sgb-dir /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/saltmarsh/assembly/ \
  --local-jobs 10

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/gtdbtk_workflow/gtdbtk_workflow.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/wild_gut/assembly/bins/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/wild_gut/gtdbtk/ \
  --abundance-type by_dataset \
  --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/wild_gut/gtdbtk/ \
  --cores 8 \
  --mags-sgb-dir /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/wild_gut/assembly/ \
  --local-jobs 10
```

