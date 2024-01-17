
library(readr)

VAF_Summary_File <- "~/R/Projects/Orien_Data_Processing/Data/Orien_IO_Rig_VAF_Summary_20240116.txt"

Orien_TMB_File <- "~/R/Projects/Orien_Data_Processing/Data/21PRJNOVA009MCC_20230731_TMB_MSI_table.csv"


vaf_summ <- as.data.frame(read_delim(VAF_Summary_File,delim = '\t', col_names = T))
Orien_TMB <- as.data.frame(read_delim(Orien_TMB_File,delim = ',', col_names = T))
colnames(Orien_TMB) <- paste0("Orien_Derived_",colnames(Orien_TMB))

Orien_TMB$WES <- gsub("^T","",str_split_fixed(Orien_TMB$Orien_Derived_Ttumor_Nnormal,"_",2)[,1])

vaf_summ2 <- merge(vaf_summ,Orien_TMB, all.x = T)



write.table(vaf_summ2,"~/R/Projects/Orien_Data_Processing/Data/Orien_IO_Rig_VAF_Summary_20240116_wMSI.txt", sep = '\t', row.names = F)
