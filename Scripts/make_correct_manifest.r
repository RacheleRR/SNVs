Liste_samples_WGS_final <- read_excel("~/Downloads/Liste_samples_WGS_final.xlsx")
View(Liste_samples_WGS_final)

# Generate the new column
Liste_samples_WGS_final$new_column <- ifelse(
Liste_samples_WGS_final$Status == "Non_Converter",
"Control",
ifelse(Liste_samples_WGS_final$Status == "UHR-NA", "UHR_NA", "Case")
)



manifest  <- Liste_samples_WGS_final[, c("Status", "Sequencing_number", "new_column")]

wgs <- read.table("~/Downloads/wgs.list", quote="\"", comment.char="")

# Filter the manifest file

filtered_manifest <- manifest[manifest$Sequencing_number %in% wgs$V1, ]
