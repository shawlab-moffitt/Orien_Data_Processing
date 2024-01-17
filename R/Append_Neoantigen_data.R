

library(readr)
library(dplyr)
library(parallel)


Neoantigen_Folders_List <- "~/R/Projects/Orien_Data_Processing/Data/Completed_NeoAntigen_Folders.txt"

Path_To_Parsed_VCF_Files <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\WES\\Somatic_VCF_Filtered_Annotated\\Patient_VCF_Annotated\\"

Output_Folder <- "M:\\dept\\Dept_BBSR\\Projects\\Shaw_Timothy\\ORIEN_Diagnosis_Metastatic_Pairs\\Orien_IO_Rig_Neoantigen_Pipeline\\VAF_Neoantigen_Patient_Summaries\\"

neoant_folders <- as.data.frame(read_delim(Neoantigen_Folders_List,delim = '\t', col_names = F))
nrow(neoant_folders)
vcf_files <- list.files(Path_To_Parsed_VCF_Files,full.names = T)
length(vcf_files)
dir.create(Output_Folder)

j <- 1
empty <- mcmapply(vcf_files, function(file) {
  
  AvatarKey <- strsplit(basename(file),"_")[[1]][4]
  
  if (any(grepl(AvatarKey,neoant_folders[,1]))) {
    print(paste0("Merging ",j,": ",AvatarKey))
    
    vcf <- as.data.frame(suppressWarnings(read_delim(file, delim = '\t', col_names = T, show_col_types = FALSE)))
    
    vcf$chr <- apply(vcf[,"Mutation",drop = F],1,function(x){strsplit(x,"[[:punct:]]")[[1]][2]})
    vcf$pos <- apply(vcf[,"Mutation",drop = F],1,function(x){strsplit(x,"[[:punct:]]")[[1]][3]})
    vcf$pos <- sub("[A-Za-z]*$","",vcf$pos)
    vcf$dna_key <- paste(vcf$chr,vcf$pos,vcf$REF,vcf$ALT,sep = "_")
    vcf <- vcf[,which(!colnames(vcf) %in% c("chr","pos"))]
    
    folders <- grep(AvatarKey,neoant_folders[,1],value = T)
    
    for (folder in folders) {
      RNASeq <- strsplit(basename(folder),"_")[[1]][2]
      WES <- strsplit(basename(folder),"_")[[1]][3]
      pm <- strsplit(basename(folder),"_")[[1]][4]
      folder_files <- list.files(folder,full.names = T)
      mhc1_file <- grep("mhc1.netmhc.txt",folder_files, value = T)
      mhc2_file <- grep("mhc2.netmhc.txt",folder_files, value = T)
      hla1_file <- grep("HLAI.txt",folder_files, value = T)
      hla2_file <- grep("HLAII.txt",folder_files, value = T)
      mhc1 <- as.data.frame(suppressWarnings(read_delim(mhc1_file, delim = '\t', col_names = T, col_select = c("#HLA_allele","dna_key","var_ic50"), show_col_types = FALSE)))
      colnames(mhc1)[1] <- "HLA_allele"
      mhc2 <- as.data.frame(suppressWarnings(read_delim(mhc2_file, delim = '\t', col_names = T, col_select = c("#HLA_allele","dna_key","var_ic50"), show_col_types = FALSE)))
      colnames(mhc2)[1] <- "HLA_allele"
      hla1 <- as.data.frame(suppressWarnings(read_delim(hla1_file, delim = '\t', col_names = F, show_col_types = FALSE)))
      hla1 <- paste(hla1[,1],collapse = ",")
      hla2 <- as.data.frame(suppressWarnings(read_delim(hla2_file, delim = '\t', col_names = F, show_col_types = FALSE)))
      hla2 <- paste(hla2[,1],collapse = ",")
      
      mhc1 <- mhc1 %>%
        group_by(dna_key) %>%
        mutate(HLA_allele = paste(HLA_allele,collapse = ","),
               var_ic50 = paste(var_ic50,collapse = ",")) %>%
        unique() %>%
        as.data.frame()
      mhc1$MHC1_Hit <- TRUE
      colnames(mhc1)[c(1,3,4)] <- c(paste(AvatarKey,pm,WES,"MHC1_HLA_allele",sep = "_"),
                                    paste(AvatarKey,pm,WES,"MHC1_Var_IC50",sep = "_"),
                                    paste(AvatarKey,pm,WES,"MHC1_Hit",sep = "_"))
      mhc2 <- mhc2 %>%
        group_by(dna_key) %>%
        mutate(HLA_allele = paste(HLA_allele,collapse = ","),
               var_ic50 = paste(var_ic50,collapse = ",")) %>%
        unique() %>%
        as.data.frame()
      mhc2$MHC2_Hit <- TRUE
      colnames(mhc2)[c(1,3,4)] <- c(paste(AvatarKey,pm,WES,"MHC2_HLA_allele",sep = "_"),
                                    paste(AvatarKey,pm,WES,"MHC2_Var_IC50",sep = "_"),
                                    paste(AvatarKey,pm,WES,"MHC2_Hit",sep = "_"))
      
      vcf <- merge(vcf,mhc1,all = T)
      vcf <- merge(vcf,mhc2,all = T)
      vcf[,paste(AvatarKey,pm,WES,"HLAI",sep = "_")] <- hla1
      vcf[,paste(AvatarKey,pm,WES,"HLAII",sep = "_")] <- hla2
      
    }
    
    write.table(vcf,paste0(Output_Folder,"Orien_IO_Rig_",AvatarKey,"_VAF_Neoantigen_Summary.txt"), sep = '\t', row.names = F)
  }
  
  j <- j+1
  return(NULL)
  
}#, mc.cores = 8
)
