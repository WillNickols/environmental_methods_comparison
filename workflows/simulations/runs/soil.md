```
conda activate simulate
hutlab load centos7/python3/biobakery_workflows/3.0.0-beta-devel-dependsUpdate

export EMAIL=willnickols@college.harvard.edu
export API_KEY=788159a1065634208a394b7635077c683c07
export CAMISIM_PATH=/n/holystore01/LABS/huttenhower_lab/Users/wnickols/CAMISIM/

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simulations/profile_building/build_profile.py \
  --known-species /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/known_species_list.tsv \
  --sgb-filepath /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/new_sgbs/sgbs/sgbs/ \
  --sgb-info /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/new_sgbs/sgbs/sgbs/SGB_info.tsv \
  --gtdbtk-relab /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/new_sgbs/gtdbtk/gtdbtk_relab.tsv \
  --profile /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/profile.tsv

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simulations/run_parameter_sweep.py \
  -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/profile.tsv \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/simulations/ \
  --config-file /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/configs/soil_config.txt \
  --local-jobs 8 \
  --grid-jobs 200 \
  --grid-partition 'shared' \
  --run-sim

python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simulations/copy_from_parameter_sweep.py --output /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/simulations/

conda deactivate
conda activate bbmap
python split_inputs_workflow.py -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/ocean/simulations/reorganized/reads/ \
  -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/ocean/simulations/reorganized/reads_split/ \
  --input-extension fastq.gz \
  --local-jobs 5
```