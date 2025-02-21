#!/bin/bash

#SBATCH --output=filter_step1_PTV_%j.log   # Output file (%j expands to job ID)
#SBATCH --error=filter_step1_PTV_%j.err    # Error file (%j expands to job ID)
#SBATCH --time=INFINITE                      # Time limit (hh:mm:ss)
#SBATCH --nodes=1                            # Number of nodes
#SBATCH --ntasks=1                           # Number of tasks
#SBATCH --cpus-per-task=1                    # Number of CPU cores per task
#SBATCH --mem=50G                           # Memory per node
#SBATCH --partition=core                     # Partition name

# Change to the directory where the VEP output files are located
cd /media/NAS/InfoGene/Rachele_Rubiu/data/VEP_output

# Run the filtering with VEP
filter_vep \
  -i merged_VEP_output.vcf.gz\
  -filter "((gnomADe_AF <= 0.01 or not gnomADe_AF) and 
             (gnomADg_AF <= 0.01 or not gnomADg_AF)) and
             (Consequence is missense_variant)" \
  --only_matched \
  -o filtered_missense.vcf