#!/bin/bash
#SBATCH --job-name=AMS
#SBATCH --partition=gpu
#SBATCH --gres=gpu:device:no_of_gpus
#SBATCH --array=1-19%5                  # Create 9 jobs
#SBATCH --nodes=1                       # node count
#SBATCH --ntasks=1                      # total number of tasks across all nodes
#SBATCH --cpus-per-task=8               # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=256GB                     # total memory allocated for all tasks
#SBATCH --time=48:00:00                 # total run time limit (HH:MM:SS)
#SBATCH --mail-type=end,fail            # send email when job begins,ends,fails
#SBATCH --mail-user=email@example.com   # email to send alerts to
#SBATCH --output /path/to/slurm_log.out        # Output file name and Location
module load aria2/1.36.0-GCCcore-11.3.0
module load HMMER/3.3.2-gompi-2022a
module load Python/3.10.4-GCCcore-11.3.0-bare

cd /path/to/workspace/alphanissense/
source /path/to/envs/venv/bin/activate

# Get the line of the input text file corresponding to the array task ID 
input_file="/path/to/variants.txt"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$input_file")

# Split the line by "_p."
IFS='_' read -r Gene Value2 <<< "$line" 

# Extract the necessary characters
ref_aa=${Value2:0:1}
alt_aa=${Value2: -1}
pos=${Value2:1:-1} 

echo "Running Alpha Missense test on $Gene _p.$ref_aa$pos_alt_aa"
echo "Processing Gene: $Gene"
echo "Ref AA: $ref_aa"
echo "Index: $pos"
echo "Alt AA: $alt_aa"
python ./omics_project/omics_models/01_run_AM.py \
/path/to/uniprotkb_AND_model_organism_9606.fasta \
$Gene \
$ref_aa \
$pos \
$alt_aa \
gof # Place subfolder name here e.g. either gof or lof
echo "-----------------SCRIPT COMPLETED-----------------"