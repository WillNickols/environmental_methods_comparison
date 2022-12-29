from anadama2 import Workflow
import os
import pandas as pd
from datetime import date
import math

workflow = Workflow(version = "0.4")
args = workflow.parse_args()

workflow.add_task_gridable("python /n/holystore01/LABS/huttenhower_lab/Users/wnickols/code/simulations/run_parameter_sweep.py -i /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/profile.tsv -o /n/holylfs05/LABS/nguyen_lab/Everyone/wnickols/simulations/soil/simulations/ --start-date 2022/05/11 --end-date 2022/07/24 --strain-n 1 --camisim-path /n/holystore01/LABS/huttenhower_lab/Users/wnickols/CAMISIM/ --email willnickols@college.harvard.edu --api-key 788159a1065634208a394b7635077c683c07 --threads-per-simulator 2 --simulators-running 4 --samples 5 --ns 10 40 80 220 --n-default 80 --sample-sizes 0.075 0.15 1.5 7.5 45 --sample-size-default 7.5 --gcs 0.5 0.6 0.7 --gc-default unrestricted --ks 40 60 80 100 200 --k-default 100 --sigmas 0.5 1 1.5 2 2.5 --sigma-default 1 --unknown-props 0 0.2 0.4 0.6 0.8 1 --unknown-prop-default 0.4 --dry-run n --grid-jobs 500 --grid-partition 'shared' --grid-options=\"--account=nguyen_lab\" --grid-scratch /n/holyscratch01/nguyen_lab/wnickols/simulations/",
    time=1200,
	mem=20000,
	cores=26,
	partition='shared')


workflow.go()
