# Run the filtering with bcftools to exclude Pathogenic and Likely_pathogenic
bcftools view -i 'INFO/gnomADe_AF<=0.01 & INFO/gnomADg_AF<=0.01 & INFO/acmg_classification_base!="Pathogenic" & INFO/acmg_classification_base!="Likely_pathogenic"' \
  "/media/rachele/One\ Touch/merged_genebe.vcf" -o non_pathogenic_variants.vcf


# Run the filtering with bcftools to exclude Pathogenic and Likely_pathogenic
bcftools view -i 'INFO/gnomADe_AF<=0.01 & INFO/gnomADg_AF<=0.01 & INFO/acmg_classification_base!="Pathogenic" & INFO/acmg_classification_base="Likely_pathogenic"' \
  "/media/rachele/One\ Touch/merged_genebe.vcf" -o non_pathogenic_variants.vcf  