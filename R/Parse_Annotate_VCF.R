# Parse and Annotate Orien AVATAR Somatic filtered VCF Files
# Author: Alyssa Obermayer

library(vcfR)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(data.table)
library(tidyr)


VCF_Files_Folder <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\WES\\Somatic_VCF_Filtered\\"

VCF_File_Pattern <- ".AF04_F1R21_F2R11_5reads.vcf.gz$"

OncoDriver_File <- "~/R/Projects/Orien_Data_Processing/Data/OncoDrivers.txt"

Output_Path <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\WES\\Somatic_VCF_Filtered_Annotated_simple\\"




dir.create(Output_Path)

VCF_Files <- list.files(VCF_Files_Folder,pattern = VCF_File_Pattern, full.names = T)
Onco_Drivers <- as.data.frame(read_delim(OncoDriver_File,delim = '\t', col_names = T))



InFrameKeys <- c("MISSENSE","IN_FRAME_DEL","IN_FRAME_INS","DE_NOVO_START_IN_FRAME")
OutFrameKeys <- c("FRAME_SHIFT_DEL","FRAME_SHIFT_INS","NONSTOP","DE_NOVO_START_OUT_FRAME ","DE_NOVO_START_OUT_FRAME","NONSENSE")
Other <- c("SPLICE_SITE","START_CODON_SNP")
check_all <- c(InFrameKeys,OutFrameKeys,Other)

split_and_name <- function(vector) {
  split_elements <- strsplit(vector, "=")
  
  # Create a named vector
  named_vector <- setNames(sapply(split_elements, function(x) {
    if (length(x) == 2) {
      x[2]  # Use the second part as the value
    } else {
      x[1]  # Use the entire element as both name and value
    }
  }), sapply(split_elements, "[", 1))
  
  return(as.data.frame(t(named_vector)))
}



k <- 1
for (file in VCF_Files) {
  
  Sample_Name <- gsub(".vcf.gz$","",basename(file))

  ## extract VCF
  vcf_in <- read.vcfR(file, verbose = F)
  vcf <- as.data.frame(cbind(vcf_in@fix,vcf_in@gt))
  vcfMeta <- as.data.frame(vcf_in@meta)
  
  ## Split INFO column to data table
  vcf_info <- strsplit(vcf$INFO,";")
  vcf_info <- lapply(vcf_info, grep, pattern = "FUNCOTATION", invert = TRUE, value = TRUE)
  vcf_info <- lapply(vcf_info, split_and_name)
  vcf_info <- rbindlist(vcf_info,fill = T)
  
  ## extract the row with the annotation labels
  vcf_func_row <- vcfMeta[which(grepl("FUNCOTATION",vcfMeta[,1])),]
  ## clean the annotation labels into column names
  vcf_func_IDs <- gsub(" ","",strsplit(vcf_func_row,":")[[1]][2])
  vcf_func_IDs <- sub("\">","",vcf_func_IDs)
  vcf_func_ColN <- strsplit(vcf_func_IDs,"[|]")[[1]]
  
  ## separate the info column to the individual columns labled with the column names derived
  vcf[vcf_func_ColN] <- str_split_fixed(vcf$INFO,"[|]",length(vcf_func_ColN))
  ## Remove initial bracket?
  vcf[,vcf_func_ColN[1]] <- sub("^..*\\[","",vcf[,vcf_func_ColN[1]])
  vcf[,vcf_func_ColN[length(vcf_func_ColN)]] <- sub("]..*..$","",vcf[,vcf_func_ColN[length(vcf_func_ColN)]])
  
  ## spread format column to individual columns
  vcf_format <- str_split(vcf[,10],":")
  vcf_format_names <- str_split(vcf[,9],":")
  for(i in seq_along(vcf_format)) names(vcf_format[[i]]) <- vcf_format_names[[i]]
  vcf_format_df <- rbindlist(lapply(vcf_format, as.data.frame.list),fill = T)
  vcf <- cbind(vcf,vcf_info,vcf_format_df)
  
  ## subset columns of interest from VCF
  vcf_sub <- vcf %>%
    select(Gencode_32_hugoSymbol,Gencode_32_genomeChange,REF,ALT,Gencode_32_variantClassification,Gencode_32_proteinChange,Gencode_32_referenceContext,
           DP,AD,SA_MAP_AF,SA_POST_PROB)
  colnames(vcf_sub)[c(1,2,5,6,7)] <- c("Gene","Mutation","Variant_Classification","Protein_Change","Reference_Context")
  vcf_sub <- vcf_sub %>%
    mutate(Variant_Frame_Classification = case_when(
      Variant_Classification %in% InFrameKeys ~ "In-Frame Event",
      Variant_Classification %in% OutFrameKeys ~ "Out-of-Frame Event",
      Variant_Classification %in% Other ~ "Other",
      !Variant_Classification %in% check_all ~ "UnalteredProtein"
    ))
  
  ## Add Onco drivers
  vcf_sub <- merge(vcf_sub,Onco_Drivers, by = "Gene", all.x = T)
  
  ## Arrange Columns of interest
  vcf_sub <- vcf_sub %>%
    relocate(Variant_Frame_Classification, .after = Variant_Classification)
  vcf_sub <- vcf_sub %>%
    relocate(any_of(c("DP","AD","SA_MAP_AF","SA_POST_PROB")), .after = last_col())
  vcf_sub$Gene_Mutation <- paste0(vcf_sub$Gene,"_",vcf_sub$Mutation)
  vcf_sub <- vcf_sub %>%
    relocate(Gene_Mutation,Gene)
  
  ## Split allele count column (AD)
  vcf_sub <- vcf_sub %>% 
    mutate(AD = strsplit(AD, ",")) %>%
    unnest(AD) %>%
    group_by(Gene_Mutation) %>%
    mutate(row = row_number()) %>%
    spread(row, AD)
  
  ## Rename and derive read count and AF
  ref_col_ind <- which(colnames(vcf_sub) == "1")
  colnames(vcf_sub)[ref_col_ind] <- "REF_READ_COUNT"
  if (length((ref_col_ind+1):(ncol(vcf_sub))) > 1) {
    alt_col_names <- paste0("ALT_READ_COUNT_",seq(length((ref_col_ind+1):(ncol(vcf_sub)))))
  } else {
    alt_col_names <- "ALT_READ_COUNT"
  }
  colnames(vcf_sub)[c((ref_col_ind+1):(ncol(vcf_sub)))] <- alt_col_names
  vcf_sub[c("REF_READ_COUNT",alt_col_names)] <- lapply(vcf_sub[c("REF_READ_COUNT",alt_col_names)], as.numeric)
  vcf_sub$ALT_READ_COUNT_SUM <- rowSums(vcf_sub[,alt_col_names],na.rm = T)
  vcf_sub$VAF <- as.numeric(sub('.*\\,', '', vcf_sub$SA_MAP_AF))
  write_delim(vcf_sub,paste0(Output_Path,Sample_Name,"_VCF_Annotated.txt"), delim = '\t')
  
  print(paste("VCFs Processed:",k,Sample_Name))
  k <- k+1
  
}



