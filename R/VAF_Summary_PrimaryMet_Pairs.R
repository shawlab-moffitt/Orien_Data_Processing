# Pair Primary/Met samples and summarize mutation counts based on VAF
# Author: Alyssa Obermayer

library(readr)
library(stringr)
library(dplyr)

VAF_Summary_Files_Dir <- "M:\\dept\\Dept_BBSR\\Projects\\Shaw_Timothy\\ORIEN_Diagnosis_Metastatic_Pairs\\Orien_IO_Rig_Neoantigen_Pipeline\\VAF_Neoantigen_Patient_Summaries_v3\\"

Files_Of_Interest_Pattern <- "MHC1hits_Neoantigen_Summary.txt"

Aneuploidy_File <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\Clinical\\Orien_IO_Rig_somatic_CNV_aneuploidy_Data_v3.txt"

VCF_Exome_Kit_File <- "~/R/Projects/Orien_Data_Processing/Data/Orien_VCF_Kit_Type.txt"

Output_File_Name <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\WES\\Orien_PairedSample_VAF_Summary_MHC1hit_01192024.txt"


VAF_files <- list.files(VAF_Summary_Files_Dir,full.names = T, pattern = Files_Of_Interest_Pattern)

if (file.exists(Aneuploidy_File)) {
  aneu <- as.data.frame(read_delim(Aneuploidy_File,delim = '\t', col_names = T))
} else {
  aneu <- NA
}

VCF_Exome_Kit_df <- as.data.frame(read_delim(VCF_Exome_Kit_File,delim = '\t', col_names = T))
WES_NIM <- VCF_Exome_Kit_df[which(VCF_Exome_Kit_df$VCF_Exome_Kit == "NIM"),"WES"]
WES_IDT <- VCF_Exome_Kit_df[which(VCF_Exome_Kit_df$VCF_Exome_Kit == "IDT"),"WES"]

header <- c("ORIENAvatarKey",
            "Primary_Sample", "Metastatic_Sample",
            "Primary_Sample_Exome_Kit","Metastatic_Sample_Exome_Kit","Matching_Primary_Met_Exome_Kit",
            ## In Frame
            "Primary_Variant_INFRAME","Primary_Variant_INFRAME_Norm","Primary_Exclusive_Variant_INFRAME","Primary_Exclusive_Variant_INFRAME_Norm",
            "Metastatic_Variant_INFRAME","Metastatic_Variant_INFRAME_Norm","Metastatic_Exclusive_Variant_INFRAME","Metastatic_Exclusive_Variant_INFRAME_Norm",
            "PrimaryMet_Shared_Variant_INFRAME","PrimaryMet_Shared_Variant_INFRAME_Norm",
            "Met_minus_Prim_INFRAME","Met_minus_Prim_INFRAME_Norm","Met_minus_Prim_Exclusive_INFRAME","Met_minus_Prim_Exclusive_INFRAME_Norm",
            ## Out of Frame
            "Primary_Variant_OUTOFFRAME","Primary_Variant_OUTOFFRAME_Norm","Primary_Exclusive_Variant_OUTOFFRAME","Primary_Exclusive_Variant_OUTOFFRAME_Norm",
            "Metastatic_Variant_OUTOFFRAME","Metastatic_Variant_OUTOFFRAME_Norm","Metastatic_Exclusive_Variant_OUTOFFRAME","Metastatic_Exclusive_Variant_OUTOFFRAME_Norm",
            "PrimaryMet_Shared_Variant_OUTOFFRAME","PrimaryMet_Shared_Variant_OUTOFFRAME_Norm",
            "Met_minus_Prim_OUTOFFRAME","Met_minus_Prim_OUTOFFRAME_Norm","Met_minus_Prim_Exclusive_OUTOFFRAME","Met_minus_Prim_Exclusive_OUTOFFRAME_Norm",
            ## Missense
            "Primary_Variant_MISSENSE","Primary_Variant_MISSENSE_Norm","Primary_Exclusive_Variant_MISSENSE","Primary_Exclusive_Variant_MISSENSE_Norm",
            "Metastatic_Variant_MISSENSE","Metastatic_Variant_MISSENSE_Norm","Metastatic_Exclusive_Variant_MISSENSE","Metastatic_Exclusive_Variant_MISSENSE_Norm",
            "PrimaryMet_Shared_Variant_MISSENSE","PrimaryMet_Shared_Variant_MISSENSE_Norm",
            "Met_minus_Prim_MISSENSE","Met_minus_Prim_MISSENSE_Norm","Met_minus_Prim_Exclusive_MISSENSE","Met_minus_Prim_Exclusive_MISSENSE_Norm",
            ## Nonsense 
            "Primary_Variant_NONSENSE","Primary_Variant_NONSENSE_Norm","Primary_Exclusive_Variant_NONSENSE","Primary_Exclusive_Variant_NONSENSE_Norm",
            "Metastatic_Variant_NONSENSE","Metastatic_Variant_NONSENSE_Norm","Metastatic_Exclusive_Variant_NONSENSE","Metastatic_Exclusive_Variant_NONSENSE_Norm",
            "PrimaryMet_Shared_Variant_NONSENSE","PrimaryMet_Shared_Variant_NONSENSE_Norm",
            "Met_minus_Prim_NONSENSE","Met_minus_Prim_NONSENSE_Norm","Met_minus_Prim_Exclusive_NONSENSE","Met_minus_Prim_Exclusive_NONSENSE_Norm",
            ## Silent 
            "Primary_Variant_SILENT","Primary_Variant_SILENT_Norm","Primary_Exclusive_Variant_SILENT","Primary_Exclusive_Variant_SILENT_Norm",
            "Metastatic_Variant_SILENT","Metastatic_Variant_SILENT_Norm","Metastatic_Exclusive_Variant_SILENT","Metastatic_Exclusive_Variant_SILENT_Norm",
            "PrimaryMet_Shared_Variant_SILENT","PrimaryMet_Shared_Variant_SILENT_Norm",
            "Met_minus_Prim_SILENT","Met_minus_Prim_SILENT_Norm","Met_minus_Prim_Exclusive_SILENT","Met_minus_Prim_Exclusive_SILENT_Norm",
            ## Total Mutation Count
            "Total_Mutation_Count_Exclusive_Primary","Total_Mutation_Count_Exclusive_Primary_Norm","Total_Mutation_Count_Exclusive_Metastatic","Total_Mutation_Count_Exclusive_Metastatic_Norm",
            "Met_minus_Prim_Total_Mutations_Exclusive","Met_minus_Prim_Total_Mutations_Exclusive_Norm",
            ## Total protein impacting mutation count
            "Total_Protein_Impacting_Mutation_Count_Exclusive_Primary","Total_Protein_Impacting_Mutation_Count_Exclusive_Primary_Norm",
            "Total_Protein_Impacting_Mutation_Count_Exclusive_Metastatic","Total_Protein_Impacting_Mutation_Count_Exclusive_Metastatic_Norm",
            "Met_minus_Prim_Total_Protein_Impacting_Mutations","Met_minus_Prim_Total_Protein_Impacting_Mutations_Norm",
            ## dN/ds
            "Primary_dN.ds","Primary_dN.ds_Norm","Primary_dN.ds_Exclusive","Primary_dN.ds_Exlusive_Norm",
            "Metastatic_dN.ds","Metastatic_dN.ds_Norm","Metastatic_dN.ds_Exclusive","Metastatic_dN.ds_Exlusive_Norm",
            "Primary_dN.ds_GreaterThan_3","Primary_dN.ds_Norm_GreaterThan_3","Primary_dN.ds_Exclusive_GreaterThan_3","Primary_dN.ds_Exlusive_Norm_GreaterThan_3",
            "Metastatic_dN.ds_GreaterThan_3","Metastatic_dN.ds_Norm_GreaterThan_3","Metastatic_dN.ds_Exclusive_GreaterThan_3","Metastatic_dN.ds_Exlusive_Norm_GreaterThan_3",
            "Met_Minus_Prim_dn.ds","Met_Minus_Prim_dn.ds_Norm","Met_Minus_Prim_dn.ds_Exclusive","Met_Minus_Prim_dn.ds_Exclusive_Norm",
            ## Aneuploidy
            "Aneuploidy_Score_Primary","Aneuploidy_Score_Metastatic","Met_minus_Prim_Aneu",
            ## HRD
            "HRD_Score_Primary","HRD_Score_Metastatic","Met_minus_Prim_HRD")
write(header,file = Output_File_Name, append = T, sep = '\t', ncolumns = 114)


for (file in VAF_files) {
  
  ## Get file info
  AvatarKey <- strsplit(basename(file),"_")[[1]][4]
  df <- as.data.frame(read_delim(file, delim = '\t', col_names = T))
  df2 <- df
  col_comb <- as.data.frame(combn(grep("_VAF$",colnames(df2),value = T), 2))
  for (pair in seq_along(col_comb)) {
    df3 <- df2[,c("Gene_Mutation","Variant_Classification","Variant_Frame_Classification",col_comb[1,pair],col_comb[2,pair])]
    met_cols <- grep("Metastatic",colnames(df3),value = T)
    pri_cols <- grep("Primary",colnames(df3),value = T)
    inFram_rows <- which(df3$Variant_Frame_Classification=="In-Frame Event")
    outFram_rows <- which(df3$Variant_Frame_Classification=="Out-of-Frame Event")
    MISSENSE_rows <- which(df3$Variant_Classification=="MISSENSE")
    NONSENSE_rows <- which(df3$Variant_Classification=="NONSENSE")
    SILENT_rows <- which(df3$Variant_Classification=="SILENT")
    mut_row_list <- list(inFram_rows = inFram_rows,outFram_rows = outFram_rows,
                         MISSENSE_rows = MISSENSE_rows,NONSENSE_rows = NONSENSE_rows,SILENT_rows = SILENT_rows)
    ## normal pair
    if (length(met_cols)==1 & length(pri_cols)==1) {
      df_data <- c(AvatarKey,pri_cols,met_cols)
      pri_WES <- strsplit(pri_cols,"_",3)[[1]][3]
      met_WES <- strsplit(met_cols,"_",3)[[1]][3]
      Primary_Sample_Exome_Kit <- ifelse(pri_WES %in% WES_NIM,"NIM","IDT")
      Metastatic_Sample_Exome_Kit <- ifelse(met_WES %in% WES_NIM,"NIM","IDT")
      Matching_Primary_Met_Exome_Kit <- ifelse(Primary_Sample_Exome_Kit==Metastatic_Sample_Exome_Kit,TRUE,FALSE)
      
      df_data <- c(df_data,Primary_Sample_Exome_Kit,Metastatic_Sample_Exome_Kit,Matching_Primary_Met_Exome_Kit)
      
      for (mut in seq_along(mut_row_list)) {
        mut_rows <- mut_row_list[[mut]]
        
        pri_mut <- nrow(df3[which(df3[mut_rows,pri_cols] != 0),])
        pri_mut_norm <- ifelse(pri_WES %in% WES_NIM,pri_mut/63.377915,pri_mut/39.671448)
        pri_mut_exl <- nrow(df3[which(df3[mut_rows,pri_cols] != 0 & df3[mut_rows,met_cols] == 0),])
        pri_mut_exl_norm <- ifelse(pri_WES %in% WES_NIM,pri_mut_exl/63.377915,pri_mut_exl/39.671448)
        
        met_mut <- nrow(df3[which(df3[mut_rows,met_cols] != 0),])
        met_mut_norm <- ifelse(met_WES %in% WES_NIM,met_mut/63.377915,met_mut/39.671448)
        met_mut_exl <- nrow(df3[which(df3[mut_rows,pri_cols] == 0 & df3[mut_rows,met_cols] != 0),])
        met_mut_exl_norm <- ifelse(met_WES %in% WES_NIM,met_mut_exl/63.377915,met_mut_exl/39.671448)
        
        Shared_mut <- nrow(df3[which(df3[mut_rows,met_cols] != 0 & df3[mut_rows,pri_cols] != 0),])
        if ((pri_WES %in% WES_NIM) & (met_WES %in% WES_NIM)) {
          Shared_mut_norm <- Shared_mut/63.377915
        } else if ((pri_WES %in% WES_IDT) & (met_WES %in% WES_IDT)) {
          Shared_mut_norm <- Shared_mut/39.671448
        } else {
          Shared_mut_norm <- Shared_mut/39.671448
        }

        met_minus_pri <- met_mut-pri_mut
        met_minus_pri_norm <- met_mut_norm-pri_mut_norm
        met_minus_pri_exl <- met_mut_exl-pri_mut_exl
        met_minus_pri_exl_norm <- met_mut_exl_norm-pri_mut_exl_norm
        
        df_data <- c(df_data,pri_mut,pri_mut_norm,pri_mut_exl,pri_mut_exl_norm,met_mut,met_mut_norm,met_mut_exl,met_mut_exl_norm,
                     Shared_mut,Shared_mut_norm,met_minus_pri,met_minus_pri_norm,met_minus_pri_exl,met_minus_pri_exl_norm)
      }
      
      ## Total Mutation Count
      Total_Pri_Mut <- nrow(df3[which(df3[,met_cols] == 0 & df3[,pri_cols] != 0),])
      Total_Pri_Mut_norm <- ifelse(pri_WES %in% WES_NIM,Total_Pri_Mut/63.377915,Total_Pri_Mut/39.671448)
      Total_Met_Mut <- nrow(df3[which(df3[,met_cols] != 0 & df3[,pri_cols] == 0),])
      Total_Met_Mut_norm <- ifelse(met_WES %in% WES_NIM,Total_Met_Mut/63.377915,Total_Met_Mut/39.671448)
      Met_minus_Prim_Total_Mutations <- Total_Met_Mut-Total_Pri_Mut
      Met_minus_Prim_norm_Total_Mutations <- Total_Met_Mut_norm-Total_Pri_Mut_norm
      
      df_data <- c(df_data,Total_Pri_Mut,Total_Pri_Mut_norm,Total_Met_Mut,Total_Met_Mut_norm,Met_minus_Prim_Total_Mutations,Met_minus_Prim_norm_Total_Mutations)
      
      ## Total protein impacting mutation count
      Prim_INFRAME <- as.numeric(df_data[which(header == "Primary_Exclusive_Variant_INFRAME")])
      Prim_INFRAME_norm <- as.numeric(df_data[which(header == "Primary_Exclusive_Variant_INFRAME_Norm")])
      Prim_OUTOFFRAME <- as.numeric(df_data[which(header == "Primary_Exclusive_Variant_OUTOFFRAME")])
      Prim_OUTOFFRAME_norm <- as.numeric(df_data[which(header == "Primary_Exclusive_Variant_OUTOFFRAME_Norm")])
      Met_INFRAME <- as.numeric(df_data[which(header == "Metastatic_Exclusive_Variant_INFRAME")])
      Met_INFRAME_norm <- as.numeric(df_data[which(header == "Metastatic_Exclusive_Variant_INFRAME_Norm")])
      Met_OUTOFFRAME <- as.numeric(df_data[which(header == "Metastatic_Exclusive_Variant_OUTOFFRAME")])
      Met_OUTOFFRAME_norm <- as.numeric(df_data[which(header == "Metastatic_Exclusive_Variant_OUTOFFRAME_Norm")])
      
      Total_Pri_ProtImp_Mut <- Prim_INFRAME+Prim_OUTOFFRAME
      Total_Pri_ProtImp_Mut_norm <- Prim_INFRAME_norm+Prim_OUTOFFRAME_norm
      Total_Met_ProtImp_Mut <- Met_INFRAME+Met_OUTOFFRAME
      Total_Met_ProtImp_Mut_norm <- Met_INFRAME_norm+Met_OUTOFFRAME_norm
      Met_minus_Prim_TotalProtImp_Mutations <- Total_Met_ProtImp_Mut-Total_Pri_ProtImp_Mut
      Met_minus_Prim_TotalProtImp_Mutations_norm <- Total_Met_ProtImp_Mut_norm-Total_Pri_ProtImp_Mut_norm
      
      df_data <- c(df_data,Total_Pri_ProtImp_Mut,Total_Pri_ProtImp_Mut_norm,Total_Met_ProtImp_Mut,Total_Met_ProtImp_Mut_norm,
                   Met_minus_Prim_TotalProtImp_Mutations,Met_minus_Prim_TotalProtImp_Mutations_norm)
      
      ## dN/ds
      Prim_MISSENSE <- as.numeric(df_data[which(header == "Primary_Variant_MISSENSE")])
      Prim_MISSENSE_norm <- as.numeric(df_data[which(header == "Primary_Variant_MISSENSE_Norm")])
      Prim_MISSENSE_exl <- as.numeric(df_data[which(header == "Primary_Exclusive_Variant_MISSENSE")])
      Prim_MISSENSE_norm_exl <- as.numeric(df_data[which(header == "Primary_Exclusive_Variant_MISSENSE_Norm")])
      Met_MISSENSE <- as.numeric(df_data[which(header == "Metastatic_Variant_MISSENSE")])
      Met_MISSENSE_norm <- as.numeric(df_data[which(header == "Metastatic_Variant_MISSENSE_Norm")])
      Met_MISSENSE_exl <- as.numeric(df_data[which(header == "Metastatic_Exclusive_Variant_MISSENSE")])
      Met_MISSENSE_norm_exl <- as.numeric(df_data[which(header == "Metastatic_Exclusive_Variant_MISSENSE_Norm")])
      
      Prim_NONSENSE <- as.numeric(df_data[which(header == "Primary_Variant_NONSENSE")])
      Prim_NONSENSE_norm <- as.numeric(df_data[which(header == "Primary_Variant_NONSENSE_Norm")])
      Prim_NONSENSE_exl <- as.numeric(df_data[which(header == "Primary_Exclusive_Variant_NONSENSE")])
      Prim_NONSENSE_norm_exl <- as.numeric(df_data[which(header == "Primary_Exclusive_Variant_NONSENSE_Norm")])
      Met_NONSENSE <- as.numeric(df_data[which(header == "Metastatic_Variant_NONSENSE")])
      Met_NONSENSE_norm <- as.numeric(df_data[which(header == "Metastatic_Variant_NONSENSE_Norm")])
      Met_NONSENSE_exl <- as.numeric(df_data[which(header == "Metastatic_Exclusive_Variant_NONSENSE")])
      Met_NONSENSE_norm_exl <- as.numeric(df_data[which(header == "Metastatic_Exclusive_Variant_NONSENSE_Norm")])
      
      Prim_SILENT <- as.numeric(df_data[which(header == "Primary_Variant_SILENT")])
      Prim_SILENT_norm <- as.numeric(df_data[which(header == "Primary_Variant_SILENT_Norm")])
      Prim_SILENT_exl <- as.numeric(df_data[which(header == "Primary_Exclusive_Variant_SILENT")])
      Prim_SILENT_norm_exl <- as.numeric(df_data[which(header == "Primary_Exclusive_Variant_SILENT_Norm")])
      Met_SILENT <- as.numeric(df_data[which(header == "Metastatic_Variant_SILENT")])
      Met_SILENT_norm <- as.numeric(df_data[which(header == "Metastatic_Variant_SILENT_Norm")])
      Met_SILENT_exl <- as.numeric(df_data[which(header == "Metastatic_Exclusive_Variant_SILENT")])
      Met_SILENT_norm_exl <- as.numeric(df_data[which(header == "Metastatic_Exclusive_Variant_SILENT_Norm")])
      
      dnds_pri <- ifelse(Prim_SILENT != 0,(Prim_MISSENSE+Prim_NONSENSE)/Prim_SILENT,NA)
      dnds_pri_GreaterThan_3 <- ifelse(dnds_pri > 3,TRUE,FALSE)
      dnds_pri_norm <- ifelse(Prim_SILENT_norm != 0,(Prim_MISSENSE_norm+Prim_NONSENSE_norm)/Prim_SILENT_norm,NA)
      dnds_pri_norm_GreaterThan_3 <- ifelse(dnds_pri_norm > 3,TRUE,FALSE)
      dnds_pri_exl <- ifelse(Prim_SILENT_exl != 0,(Prim_MISSENSE_exl+Prim_NONSENSE_exl)/Prim_SILENT_exl,NA)
      dnds_pri_exl_GreaterThan_3 <- ifelse(dnds_pri_exl > 3,TRUE,FALSE)
      dnds_pri_norm_exl <- ifelse(Prim_SILENT_norm_exl != 0,(Prim_MISSENSE_norm_exl+Prim_NONSENSE_norm_exl)/Prim_SILENT_norm_exl,NA)
      dnds_pri_norm_exl_GreaterThan_3 <- ifelse(dnds_pri_norm_exl > 3,TRUE,FALSE)
      
      dnds_met <- ifelse(Met_SILENT != 0,(Met_MISSENSE+Met_NONSENSE)/Met_SILENT,NA)
      dnds_met_GreaterThan_3 <- ifelse(dnds_met > 3,TRUE,FALSE)
      dnds_met_norm <- ifelse(Met_SILENT_norm != 0,(Met_MISSENSE_norm+Met_NONSENSE_norm)/Met_SILENT_norm,NA)
      dnds_met_norm_GreaterThan_3 <- ifelse(dnds_met_norm > 3,TRUE,FALSE)
      dnds_met_exl <- ifelse(Met_SILENT_exl != 0,(Met_MISSENSE_exl+Met_NONSENSE_exl)/Met_SILENT_exl,NA)
      dnds_met_exl_GreaterThan_3 <- ifelse(dnds_met_exl > 3,TRUE,FALSE)
      dnds_met_norm_exl <- ifelse(Met_SILENT_norm_exl != 0,(Met_MISSENSE_norm_exl+Met_NONSENSE_norm_exl)/Met_SILENT_norm_exl,NA)
      dnds_met_norm_exl_GreaterThan_3 <- ifelse(dnds_met_norm_exl > 3,TRUE,FALSE)
      
      Met_minus_Prim_dnds <- dnds_met-dnds_pri
      Met_minus_Prim_dnds_norm <- dnds_met_norm-dnds_pri_norm
      Met_minus_Prim_dnds_exl <- dnds_met_exl-dnds_pri_exl
      Met_minus_Prim_dnds_norm_exl <- dnds_met_norm_exl-dnds_pri_norm_exl
      
      df_data <- c(df_data,dnds_pri,dnds_pri_norm,dnds_pri_exl,dnds_pri_norm_exl,dnds_met,dnds_met_norm,dnds_met_exl,dnds_met_norm_exl,
                   dnds_pri_GreaterThan_3,dnds_pri_norm_GreaterThan_3,dnds_pri_exl_GreaterThan_3,dnds_pri_norm_exl_GreaterThan_3,
                   dnds_met_GreaterThan_3,dnds_met_norm_GreaterThan_3,dnds_met_exl_GreaterThan_3,dnds_met_norm_exl_GreaterThan_3,
                   Met_minus_Prim_dnds,Met_minus_Prim_dnds_norm,Met_minus_Prim_dnds_exl,Met_minus_Prim_dnds_norm_exl)
      
      if (exists("aneu")) {
        ## Aneuploidy
        pri_aneu <- aneu[which(aneu$WES == strsplit(pri_cols,"_")[[1]][3]),"seg_Level_CNV_Aneuploidy"]
        met_aneu <- aneu[which(aneu$WES == strsplit(met_cols,"_")[[1]][3]),"seg_Level_CNV_Aneuploidy"]
        Met_minus_Prim_Aneu <- met_aneu-pri_aneu
        ## HRD
        pri_HRD <- aneu[which(aneu$WES == strsplit(pri_cols,"_")[[1]][3]),"HRD_scarHRD"]
        met_HRD <- aneu[which(aneu$WES == strsplit(met_cols,"_")[[1]][3]),"HRD_scarHRD"]
        Met_minus_Prim_HRD <- met_HRD-pri_HRD
        
        df_data <- c(df_data,pri_aneu,met_aneu,Met_minus_Prim_Aneu,pri_HRD,met_HRD,Met_minus_Prim_HRD)
      }
      
    } 
    
    write(df_data,file = Output_File_Name, append = T, sep = '\t', ncolumns = 114)
    
  }
  
}
