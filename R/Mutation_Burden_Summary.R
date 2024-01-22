
library(readr)
library(stringr)
library(dplyr)


VCF_Anno_Files_Folder <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\WES\\Somatic_VCF_Filtered_Annotated_simple\\"

Aneuploidy_File <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\Clinical\\Orien_IO_Rig_somatic_CNV_aneuploidy_Data_v3.txt"

VCF_Exome_Kit_File <- "~/R/Projects/Orien_Data_Processing/Data/Orien_VCF_Kit_Type.txt"

Output_File <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\WES\\Orien_IO_Rig_VAF_Summary_20240116_v3.txt"


aneu <- as.data.frame(read_delim(Aneuploidy_File,delim = '\t', col_names = T))
VCF_Anno_files <- list.files(VCF_Anno_Files_Folder,full.names = T)
VCF_Exome_Kit_df <- as.data.frame(read_delim(VCF_Exome_Kit_File,delim = '\t', col_names = T))

WES_NIM <- VCF_Exome_Kit_df[which(VCF_Exome_Kit_df$VCF_Exome_Kit == "NIM"),"WES"]
WES_IDT <- VCF_Exome_Kit_df[which(VCF_Exome_Kit_df$VCF_Exome_Kit == "IDT"),"WES"]

mutation_cols <- c("Variant_InFrame","Variant_OutOfFrame","Variant_MISSENSE","Variant_NONSENSE","Variant_SILENT",
                   "Total_Mutation_Count","Total_Protein_Impacting_Mutation_Count")

header <- c("WES","VCF_Sample_Name","VCF_Exome_Kit",
            "Variant_InFrame","Variant_InFrame_Norm","Variant_OutOfFrame","Variant_OutOfFrame_Norm",
            "Variant_MISSENSE","Variant_MISSENSE_Norm","Variant_NONSENSE","Variant_NONSENSE_Norm","Variant_SILENT","Variant_SILENT_Norm",
            "dN.ds","dN.ds_GreaterThan_3","Aneuploidy_Score","HRD_scarHRD","Telomeric.AI_scarHRD","LST_scarHRD","HRD.sum_scarHRD",
            "Total_Mutation_Count","Total_Mutation_Count_Norm","Total_Protein_Impacting_Mutation_Count","Total_Protein_Impacting_Mutation_Count_Norm")
write(header,file = Output_File, append = T, sep = '\t', ncolumns = 24)


for (file in VCF_Anno_files) {
  
  VCF_Sample_Name <- basename(file)
  WES <- gsub("^T","",gsub("_..*","",VCF_Sample_Name))
  #WES <- gsub("Patient","",strsplit(VCF_Sample_Name,"_")[[1]][6])
  df <- as.data.frame(read_delim(file, delim = '\t', col_names = T))
  VCF_Exome_Kit <- VCF_Exome_Kit_df[which(VCF_Exome_Kit_df$WES == WES),"VCF_Exome_Kit"]
  
  inFrame_rows <- which(df$Variant_Frame_Classification=="In-Frame Event")
  outFrame_rows <- which(df$Variant_Frame_Classification=="Out-of-Frame Event")
  MISSENSE_rows <- which(df$Variant_Classification=="MISSENSE")
  NONSENSE_rows <- which(df$Variant_Classification=="NONSENSE")
  SILENT_rows <- which(df$Variant_Classification=="SILENT")
  
  Variant_InFrame <- nrow(df[which(df[inFrame_rows,"VAF"] != 0),])
  Variant_OutOfFrame <- nrow(df[which(df[outFrame_rows,"VAF"] != 0),])
  Variant_MISSENSE <- nrow(df[which(df[MISSENSE_rows,"VAF"] != 0),])
  Variant_NONSENSE <- nrow(df[which(df[NONSENSE_rows,"VAF"] != 0),])
  Variant_SILENT <- nrow(df[which(df[SILENT_rows,"VAF"] != 0),])
  
  Variant_InFrame_Norm <- ifelse(WES %in% WES_NIM,Variant_InFrame/63.377915,Variant_InFrame/39.671448)
  Variant_OutOfFrame_Norm <- ifelse(WES %in% WES_NIM,Variant_OutOfFrame/63.377915,Variant_OutOfFrame/39.671448)
  Variant_MISSENSE_Norm <- ifelse(WES %in% WES_NIM,Variant_MISSENSE/63.377915,Variant_MISSENSE/39.671448)
  Variant_NONSENSE_Norm <- ifelse(WES %in% WES_NIM,Variant_NONSENSE/63.377915,Variant_NONSENSE/39.671448)
  Variant_SILENT_Norm <- ifelse(WES %in% WES_NIM,Variant_SILENT/63.377915,Variant_SILENT/39.671448)
  
  dN.ds <- (Variant_MISSENSE+Variant_NONSENSE)/Variant_SILENT
  dN.ds_GreaterThan_3 <- ifelse(dN.ds > 3,TRUE,FALSE)
  Aneuploidy_Score <- aneu[which(aneu$WES == WES),"seg_Level_CNV_Aneuploidy"]
  HRD_scarHRD <- aneu[which(aneu$WES == WES),"HRD_scarHRD"]
  Telomeric.AI_scarHRD <- aneu[which(aneu$WES == WES),"Telomeric.AI_scarHRD"]
  LST_scarHRD <- aneu[which(aneu$WES == WES),"LST_scarHRD"]
  HRD.sum_scarHRD <- aneu[which(aneu$WES == WES),"HRD.sum_scarHRD"]
  
  Total_Mutation_Count <- nrow(df[which(df[,"VAF"] != 0),])
  Total_Protein_Impacting_Mutation_Count <- Variant_InFrame+Variant_OutOfFrame
  Total_Mutation_Count_Norm <- ifelse(WES %in% WES_NIM,Total_Mutation_Count/63.377915,Total_Mutation_Count/39.671448)
  Total_Protein_Impacting_Mutation_Count_Norm <- ifelse(WES %in% WES_NIM,Total_Protein_Impacting_Mutation_Count/63.377915,Total_Protein_Impacting_Mutation_Count/39.671448)
  
  Mutation_Burden_Summary <- c(WES,VCF_Sample_Name,VCF_Exome_Kit,
                               Variant_InFrame,Variant_InFrame_Norm,Variant_OutOfFrame,Variant_OutOfFrame_Norm,
                               Variant_MISSENSE,Variant_MISSENSE_Norm,Variant_NONSENSE,Variant_NONSENSE_Norm,Variant_SILENT,Variant_SILENT_Norm,
                               dN.ds,dN.ds_GreaterThan_3,Aneuploidy_Score,HRD_scarHRD,Telomeric.AI_scarHRD,LST_scarHRD,HRD.sum_scarHRD,
                               Total_Mutation_Count,Total_Mutation_Count_Norm,Total_Protein_Impacting_Mutation_Count,Total_Protein_Impacting_Mutation_Count_Norm)
  write(Mutation_Burden_Summary,file = Output_File, append = T, sep = '\t', ncolumns = 24)
  

  
}
