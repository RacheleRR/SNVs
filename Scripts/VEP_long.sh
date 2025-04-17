#!/bin/bash
#SBATCH --job-name=VEP_Rachele    # Nom du job
#SBATCH --output=vep_%j.out       # Nom du fichier de sortie
#SBATCH --error=vep_%j.err        # Nom du fichier d'erreur
#SBATCH --ntasks=1                # Nombre de tâches (1 tâche dans ce cas)
#SBATCH --nodelist=macgyver  
#SBATCH --cpus-per-task=4         # Nombre de CPU par tâche (4 CPU)
#SBATCH --mem=200G                # Mémoire à allouer
#SBATCH --time=INFINITE           # Temps d'exécution illimité
#SBATCH --partition=core          # Partition à utiliser (partition "core" par défaut)

# Dossier contenant vos fichiers d'entrée
input_dir="/home/benjamin/NAS/InfoGene/Rachele/SNVs/step13_rs"

# Dossier de sortie pour les fichiers VEP générés
output_dir="/home/benjamin/NAS/InfoGene/Rachele/SNVs/VEP_output"

# Créer le répertoire de sortie s'il n'existe pas
mkdir -p "$output_dir"

# Boucle sur tous les fichiers .vcf.gz dans le dossier d'entrée
for input_file in "$input_dir"/rs_normalized_clean_merged_all_samples_chr*.vcf.gz; do
    # Extraire le nom de base du fichier pour créer un nom de sortie correspondant
    base_name=$(basename "$input_file" .vcf.gz)
    output_file="$output_dir/${base_name}_VEP_output.vcf"
    
    echo "Processing input: $input_file"
    echo "Saving to output: $output_file"

    vep -i "$input_file" \
        --cache \
        --offline \
        --sift b \
        --fork 4 \
        --vcf \
        --assembly GRCh38 \
        --fasta /home/benjamin/NAS/InfoGene/Rachele/VEP_database/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
        --output_file "$output_file" \
        --dir_cache /home/benjamin/NAS/InfoGene/Rachele/VEP_database/ \
        --cache_version 112 \
        --af_gnomade \
        --af_gnomadg \
        --plugin G2P,file=/home/benjamin/NAS/InfoGene/Rachele/VEP_database/DDG2P.csv \
        --plugin TSSDistance \
        --plugin AlphaMissense,file=/home/benjamin/NAS/InfoGene/Rachele/VEP_database/AlphaMissense_hg38.tsv.gz \
        --plugin MPC \
        --plugin TissueExpression \
        --plugin Enformer,file=/home/benjamin/NAS/InfoGene/Rachele/VEP_database/enformer_grch38.vcf.gz \
        --plugin pLI,/home/benjamin/NAS/InfoGene/Rachele/VEP_database/pLI_values.txt \
        --plugin dbNSFP,/home/benjamin/NAS/InfoGene/Rachele/VEP_database/dbNSFP4.9/dbNSFP4.9a_grch38.gz,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_pred,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,MVP_score,MPC_score,PrimateAI_score,PrimateAI_pred,ClinPred_score,Aloft_pred,Aloft_Confidence,CADD_phred,GERP++_RS,phyloP100way_vertebrate,phyloP470way_mammalian,phyloP17way_primate,phastCons100way_vertebrate,phastCons470way_mammalian,phastCons17way_primate,VEST4_score,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,Reliability_index,MetaRNN_score,MetaRNN_pred,M-CAP_score,M-CAP_pred,REVEL_score,MutPred_score,gMVP_score,DEOGEN2_score,DEOGEN2_pred,BayesDel_addAF_score,BayesDel_addAF_pred,BayesDel_noAF_score,BayesDel_noAF_pred,VARITY_R_score,VARITY_ER_score,CADD_raw,DANN_score,fathmm-MKL_coding_score,fathmm-MKL_coding_rankscore,fathmm-MKL_coding_pred,fathmm-MKL_coding_group,fathmm-XF_coding_score,fathmm-XF_coding_rankscore,fathmm-XF_coding_pred,Eigen-raw_coding,Eigen-PC-raw_coding,GenoCanyon_score,integrated_fitCons_score,integrated_confidence_value,LINSIGHT,clinvar_clnsig,clinvar_id,GTEx_V8_eQTL_gene,GTEx_V8_eQTL_tissue,GTEx_V8_sQTL_gene,GTEx_V8_sQTL_tissue,LIST-S2_score,LIST-S2_pred,EVE_score,EVE_Class10_pred,ESM1b_score,ESM1b_pred,Aloft_Fraction_transcripts_affected,Aloft_pred,Aloft_Confidence,PHACTboost_score,PHACTboost_rankscore,MutFormer_score,MutFormer_rankscore,GERP_91_mammals,MutScore_score \
        --dir_plugins /opt/ensembl-vep/Plugins/ \
        --force_overwrite \
        --hgvs \
        --minimal \
        --allele_number

    echo "Done with $input_file"
done

#add transcript_match=1 to         --plugin dbNSFP,/home/benjamin/NAS/InfoGene/Rachele/VEP_database/dbNSFP4.9/dbNSFP4.9a_grch38.gz  like so --plugin dbNSFP,transcript_match=1,/home/benjamin/NAS/InfoGene/Rachele/VEP_database/dbNSFP4.9/dbNSFP4.9a_grch38.gz,SIFT_score,