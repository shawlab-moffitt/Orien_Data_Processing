

library(readr)
library(dplyr)
library(parallel)
library(purrr)
library(stringr)


Neoantigen_Folders_List <- "~/R/Projects/Orien_Data_Processing/Data/Completed_NeoAntigen_Folders.txt"

Path_To_Parsed_VCF_Files <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\WES\\Somatic_VCF_Filtered_Annotated\\Patient_VCF_Annotated\\"

Output_Folder <- "M:\\dept\\Dept_BBSR\\Projects\\Shaw_Timothy\\ORIEN_Diagnosis_Metastatic_Pairs\\Orien_IO_Rig_Neoantigen_Pipeline\\VAF_Neoantigen_Patient_Summaries_v3\\"

neoant_folders <- as.data.frame(read_delim(Neoantigen_Folders_List,delim = '\t', col_names = F))
nrow(neoant_folders)
vcf_files <- list.files(Path_To_Parsed_VCF_Files,full.names = T)
length(vcf_files)
suppressWarnings(dir.create(Output_Folder))

empty <- mclapply(vcf_files, function(file) {
  
  AvatarKey <- strsplit(basename(file),"_")[[1]][4]
  
  # Check f AvatarKey was ran throuhg neoantigen pipeline
  if (any(grepl(AvatarKey,neoant_folders[,1]))) {
    print(paste0("Merging: ",AvatarKey))
    
    vcf <- as.data.frame(suppressWarnings(read_delim(file, delim = '\t', col_names = T, show_col_types = FALSE)))
    
    # Create dna_key column from mutation information for merge
    vcf$chr <- apply(vcf[,"Mutation",drop = F],1,function(x){strsplit(x,"[[:punct:]]")[[1]][2]})
    vcf$pos <- apply(vcf[,"Mutation",drop = F],1,function(x){strsplit(x,"[[:punct:]]")[[1]][3]})
    vcf$pos <- sub("[A-Za-z]*$","",vcf$pos)
    vcf$dna_key <- paste(vcf$chr,vcf$pos,vcf$REF,vcf$ALT,sep = "_")
    vcf <- vcf[,which(!colnames(vcf) %in% c("chr","pos"))]
    
    # Get neoantigen folders
    folders <- grep(AvatarKey,neoant_folders[,1],value = T)
    if (length(folders) > 0) {
      full_vcf <- list()
      mhc1_vcf <- list()
      mhc2_vcf <- list()
      j <- 1
      for (folder in folders) {
        # Get folder information and read through required files
        RNASeq <- strsplit(basename(folder),"_")[[1]][2]
        WES <- strsplit(basename(folder),"_")[[1]][3]
        pm <- strsplit(basename(folder),"_")[[1]][4]
        folder_files <- list.files(folder,full.names = T)
        mhc1_file <- grep("mhc1.netmhc.txt",folder_files, value = T)
        mhc2_file <- grep("mhc2.netmhc.txt",folder_files, value = T)
        hla1_file <- grep("HLAI.txt",folder_files, value = T)
        hla2_file <- grep("HLAII.txt",folder_files, value = T)
        mhc1 <- as.data.frame(suppressWarnings(read_delim(mhc1_file, delim = '\t', col_names = T, col_select = c("#HLA_allele","dna_key","prot_change","var_ic50"), show_col_types = FALSE)))
        colnames(mhc1)[1] <- "HLA_allele"
        mhc2 <- as.data.frame(suppressWarnings(read_delim(mhc2_file, delim = '\t', col_names = T, col_select = c("#HLA_allele","dna_key","prot_change","var_ic50"), show_col_types = FALSE)))
        colnames(mhc2)[1] <- "HLA_allele"
        hla1 <- as.data.frame(suppressWarnings(read_delim(hla1_file, delim = '\t', col_names = F, show_col_types = FALSE)))
        hla1 <- paste(hla1[,1],collapse = ",")
        hla2 <- as.data.frame(suppressWarnings(read_delim(hla2_file, delim = '\t', col_names = F, show_col_types = FALSE)))
        hla2 <- paste(hla2[,1],collapse = ",")
        
        
        mhc1_pc <- merge(mhc1,vcf[,c("dna_key","Protein_Change")], all = T)
        mhc2_pc <- merge(mhc2,vcf[,c("dna_key","Protein_Change")], all = T)
        
        mhc1_pc_match <- mhc1_pc %>% 
          filter(map2_lgl(prot_change, gsub("fs","",Protein_Change),  str_detect)) %>%
          group_by(dna_key) %>%
          mutate(HLA_allele = paste(unique(HLA_allele),collapse = ","),
                 prot_change = paste(unique(prot_change),collapse = ","),
                 var_ic50 = paste(unique(var_ic50),collapse = ",")) %>%
          unique() %>%
          select(!Protein_Change) %>%
          as.data.frame()
        
        mhc2_pc_match <- mhc2_pc %>% 
          filter(map2_lgl(prot_change, gsub("fs","",Protein_Change),  str_detect)) %>%
          group_by(dna_key) %>%
          mutate(HLA_allele = paste(unique(HLA_allele),collapse = ","),
                 prot_change = paste(unique(prot_change),collapse = ","),
                 var_ic50 = paste(unique(var_ic50),collapse = ",")) %>%
          unique() %>%
          select(!Protein_Change) %>%
          as.data.frame()
        
        
        
        # Combine HLA allele and var IC50 information based on same mutation
        mhc1 <- mhc1 %>%
          group_by(dna_key) %>%
          mutate(HLA_allele = paste(unique(HLA_allele),collapse = ","),
                 prot_change = paste(unique(prot_change),collapse = ","),
                 var_ic50 = paste(unique(var_ic50),collapse = ",")) %>%
          unique() %>%
          as.data.frame()
        mhc1$MHC1_Hit <- TRUE
        colnames(mhc1)[c(1,3,4,5)] <- c(paste(AvatarKey,pm,WES,"MHC1_HLA_allele",sep = "_"),
                                        paste(AvatarKey,pm,WES,"MHC1_Protein_Change",sep = "_"),
                                        paste(AvatarKey,pm,WES,"MHC1_Var_IC50",sep = "_"),
                                        paste(AvatarKey,pm,WES,"MHC1_Hit",sep = "_"))
        mhc2 <- mhc2 %>%
          group_by(dna_key) %>%
          mutate(HLA_allele = paste(unique(HLA_allele),collapse = ","),
                 prot_change = paste(unique(prot_change),collapse = ","),
                 var_ic50 = paste(unique(var_ic50),collapse = ",")) %>%
          unique() %>%
          as.data.frame()
        mhc2$MHC2_Hit <- TRUE
        colnames(mhc2)[c(1,3,4,5)] <- c(paste(AvatarKey,pm,WES,"MHC2_HLA_allele",sep = "_"),
                                        paste(AvatarKey,pm,WES,"MHC2_Protein_Change",sep = "_"),
                                        paste(AvatarKey,pm,WES,"MHC2_Var_IC50",sep = "_"),
                                        paste(AvatarKey,pm,WES,"MHC2_Hit",sep = "_"))
        
        # Merge data
        vcf_full <- merge(vcf,mhc1,all = T)
        vcf_full <- merge(vcf_full,mhc2,all = T)
        vcf_full[,paste(AvatarKey,pm,WES,"HLAI",sep = "_")] <- hla1
        vcf_full[,paste(AvatarKey,pm,WES,"HLAII",sep = "_")] <- hla2
        full_vcf[[paste("vcf_full_",j)]] <- vcf_full
        
        if (nrow(mhc1_pc_match) > 0) {
          mhc1_pc_match$MHC1_Hit <- TRUE
          colnames(mhc1_pc_match)[c(2:5)] <- c(paste(AvatarKey,pm,WES,"MHC1_HLA_allele",sep = "_"),
                                               paste(AvatarKey,pm,WES,"MHC1_Protein_Change",sep = "_"),
                                               paste(AvatarKey,pm,WES,"MHC1_Var_IC50",sep = "_"),
                                               paste(AvatarKey,pm,WES,"MHC1_Hit",sep = "_"))
          vcf_mhc1 <- merge(vcf,mhc1_pc_match,all.y = T)
          vcf_mhc1 <- merge(vcf_mhc1,mhc2_pc_match,all.x = T)
          vcf_mhc1[,paste(AvatarKey,pm,WES,"HLAI",sep = "_")] <- hla1
          vcf_mhc1[,paste(AvatarKey,pm,WES,"HLAII",sep = "_")] <- hla2
          mhc1_vcf[[paste("mhc1_vcf_",j)]] <- vcf_mhc1
        }
        
        
        if (nrow(mhc2_pc_match) > 0) {
          mhc2_pc_match$MHC2_Hit <- TRUE
          colnames(mhc2_pc_match)[c(2:5)] <- c(paste(AvatarKey,pm,WES,"MHC2_HLA_allele",sep = "_"),
                                               paste(AvatarKey,pm,WES,"MHC2_Protein_Change",sep = "_"),
                                               paste(AvatarKey,pm,WES,"MHC2_Var_IC50",sep = "_"),
                                               paste(AvatarKey,pm,WES,"MHC2_Hit",sep = "_"))
          vcf_mhc2 <- merge(vcf,mhc2_pc_match,all.y = T)
          vcf_mhc2 <- merge(vcf_mhc2,mhc1_pc_match,all.x = T)
          vcf_mhc2[,paste(AvatarKey,pm,WES,"HLAI",sep = "_")] <- hla1
          vcf_mhc2[,paste(AvatarKey,pm,WES,"HLAII",sep = "_")] <- hla2
          mhc2_vcf[[paste("mhc2_vcf_",j)]] <- vcf_mhc2
        }
        
        j <- j+1
      }
      
      full_vcf_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all = TRUE),
                                full_vcf)
      mhc1_vcf_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all = TRUE),
                                mhc1_vcf)
      mhc2_vcf_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, all = TRUE),
                                mhc2_vcf)
      
      #mhc1_hit_cols <- grep("MHC1_Hit",colnames(vcf_full),value = T)
      #mhc2_hit_cols <- grep("MHC2_Hit",colnames(vcf_full),value = T)
      #
      #vcf_mhc1_hit <- vcf_full[which(rowSums(vcf_full[,mhc1_hit_cols, drop = F], na.rm = T) > 0),]
      #vcf_mhc2_hit <- vcf_full[which(rowSums(vcf_full[,mhc2_hit_cols, drop = F], na.rm = T) > 0),]
      
      # Write data
      write.table(full_vcf_merged,paste0(Output_Folder,"Orien_IO_Rig_",AvatarKey,"_VAF_Neoantigen_Summary.txt"), sep = '\t', row.names = F)
      write.table(mhc1_vcf_merged,paste0(Output_Folder,"Orien_IO_Rig_",AvatarKey,"_VAF_MHC1hits_Neoantigen_Summary.txt"), sep = '\t', row.names = F)
      write.table(mhc2_vcf_merged,paste0(Output_Folder,"Orien_IO_Rig_",AvatarKey,"_VAF_MHC2hits_Neoantigen_Summary.txt"), sep = '\t', row.names = F)
      
    }
    
  }
  
  return(NULL)
  
}#, mc.cores = 8
)
