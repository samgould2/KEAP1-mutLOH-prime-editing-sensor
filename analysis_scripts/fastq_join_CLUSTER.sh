#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1                      # Number of tasks (really number of CPU Cores/task). Don't specify more than 16 unless approved by the system admin
#SBATCH --array=1-16    #change this to match up with the number of parallel computing jobs; also generate a config file 
#SBATCH --mail-type=END           # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE 
#SBATCH --mail-user=samgould@mit.edu           # Email to which notifications will be sent. Equivalent to the -M option in SGE. You must replace [] with your email address.
#SBATCH --exclude=c[5-22]
#SBATCH --nice=100000
#############################################

#module load miniconda3/v4
#source /home/software/conda/miniconda3/bin/condainit

#create a new virtual environment if needed; this one is fine from previously
#conda activate /home/samgould/.conda/envs/sensor_lib_sg
cd /net/bmc-lab2/data/lab/sanchezrivera/samgould/240520San

#access the config file
#config=./config.txt
config=/net/bmc-lab2/data/lab/sanchezrivera/samgould/240520San/keap1_config.txt

# Extract R1_FILE name for the current $SLURM_ARRAY_TASK_ID
R1_FILE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
R2_FILE=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
folder_name=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)


/net/bmc-lab2/data/lab/sanchezrivera/samgould/fastq-join/fastq-join ${R1_FILE} ${R2_FILE} -o ./fastq-join/${folder_name}_%.fastq
#test 
#echo "This is array task ${SLURM_ARRAY_TASK_ID}, the sample name is ${R1_FILE} and ${R2_FILE} and the output folder is ${folder_name}" >> output.txt