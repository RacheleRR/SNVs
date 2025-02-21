#!/bin/bash

#SBATCH --output=filter_step1_MPC_alpha_%j.log   # Output file (%j expands to job ID)
#SBATCH --error=filter_step1_MPC_alpha_%j.err    # Error file (%j expands to job ID)
#SBATCH --time=INFINITE                          # Time limit (hh:mm:ss)
#SBATCH --nodes=1                                # Number of nodes
#SBATCH --ntasks=1                               # Number of tasks
#SBATCH --cpus-per-task=1                        # Number of CPU cores per task
#SBATCH --mem=50G                                # Memory per node
#SBATCH --partition=core                         # Partition name

# Define the input file
INPUT_FILE="/media/NAS/InfoGene/Rachele_Rubiu/data/VEP_output/merged_VEP_output.vcf"

# Create the output directories for filtered results
mkdir -p filtered_round_1 filtered_round_2

# Round 1 Filtering
echo "Starting Round 1 Filtering..."

# Run the filtering with VEP for MPC Alpha Missense
filter_vep \
  -i "$INPUT_FILE" \
  -filter "((gnomADe_AF <= 0.01 or not gnomADe_AF) and 
             (gnomADg_AF <= 0.01 or not gnomADg_AF)) and 
             (MPC_score >= 2) and 
             (am_pathogenicity >= 0.98)" \
  --only_matched \
  -o filtered_round_1/pred_patho_MPC_ALPHAMISSENSE_VEP_filtered.vcf

# Run the filtering with VEP for MPC
filter_vep \
  -i "$INPUT_FILE" \
  -filter "((gnomADe_AF <= 0.01 or not gnomADe_AF) and 
             (gnomADg_AF <= 0.01 or not gnomADg_AF)) and 
             (MPC_score < 2)" \
  --only_matched \
  -o filtered_round_1/pred_NON_patho_MPC_VEP_filtered.vcf  

# Run the filtering with VEP for PTV
filter_vep \
  -i "$INPUT_FILE" \
  -filter "((gnomADe_AF <= 0.01 or not gnomADe_AF) and 
             (gnomADg_AF <= 0.01 or not gnomADg_AF)) and
             (Consequence is frameshift_variant or 
              Consequence is stop_gained or 
              Consequence is splice_acceptor_variant or 
              Consequence is splice_donor_variant or 
              Consequence is start_lost)" \
  --only_matched \
  -o filtered_round_1/pred_PTV_VEP_filtered.vcf

echo "Round 1 Filtering completed."

# Round 2 Filtering
echo "Starting Round 2 Filtering..."

# Run the filtering with VEP for MPC Alpha Missense with pLI
filter_vep \
  -i "$INPUT_FILE" \
  -filter "((gnomADe_AF <= 0.01 or not gnomADe_AF) and 
             (gnomADg_AF <= 0.01 or not gnomADg_AF)) and 
             (pLI_gene_value > 0.9) and 
             (MPC_score >= 2) and 
             (am_pathogenicity >= 0.98)" \
  --only_matched \
  -o filtered_round_2/pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered.vcf

# Run the filtering with VEP for MPC with pLI
filter_vep \
  -i "$INPUT_FILE" \
  -filter "((gnomADe_AF <= 0.01 or not gnomADe_AF) and 
             (gnomADg_AF <= 0.01 or not gnomADg_AF)) and 
             (pLI_gene_value > 0.9) and 
             (MPC_score < 2)" \
  --only_matched \
  --force_overwrite \
  -o filtered_round_2/pred_NON_patho_MPC_pLI_VEP_filtered.vcf

# Run the filtering with VEP for PTV with pLI
filter_vep \
  -i "$INPUT_FILE" \
  -filter "((gnomADe_AF <= 0.01 or not gnomADe_AF) and 
             (gnomADg_AF <= 0.01 or not gnomADg_AF)) and 
             (pLI_gene_value > 0.9) and 
             (Consequence is frameshift_variant or 
              Consequence is stop_gained or 
              Consequence is splice_acceptor_variant or 
              Consequence is splice_donor_variant or 
              Consequence is start_lost)" \
  --only_matched \
  -o filtered_round_2/pred_PTV_pLI_VEP_filtered.vcf

echo "Round 2 Filtering completed."
