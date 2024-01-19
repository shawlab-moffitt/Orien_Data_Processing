# Libraries --------------------------------------------------------------------
library(readr)
library(dplyr)
library(stringr)


# Input Files ------------------------------------------------------------------

VAF_Summary <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\WES\\Orien_PairedSample_VAF_Summary_01192024.txt"

Clinical_Data <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\Clinical\\Orien_IO_Rig_NovaGrant_Merged_Clinical_All_v4.txt"

Specimen_Site_Data <- "~/R/ShinyAppDev/Orien_IO_Rig_Pateint_Event_Timeline_App/Orien_specimensite_simplified.txt"

ER_Stress_ssGSEA_Scores <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\RNAseq\\ERStress_Scores\\Orien_IO_Rig_Paired_ERstress_Sig_noNorm.txt"

Output_File <- "J:\\Proj_IO_RIG_NOVA_Grant\\Harmonized_Softlinks\\work\\Alyssa_Work\\WES\\Orien_PairedSample_VAF_Summary_Extended_01182024.txt"



# Read In Files ----------------------------------------------------------------

vaf <- as.data.frame(read_delim(VAF_Summary, delim = '\t', col_names = T))
clin <- as.data.frame(read_delim(Clinical_Data, delim = '\t', col_names = T))
site <- as.data.frame(read_delim(Specimen_Site_Data, delim = '\t', col_names = F))
ER_Stress <- as.data.frame(read_delim(ER_Stress_ssGSEA_Scores,delim = '\t', col_names = T))



# Append Disease Type ----------------------------------------------------------

clin_dis <- clin %>%
  filter(!is.na(WES)) %>%
  select(WES,`Disease Type`) %>%
  group_by(WES) %>%
  mutate(`Disease Type` = paste(unique(`Disease Type`), collapse = ',')) %>%
  unique() %>%
  as.data.frame()

clin_dis$`Disease Type` <- gsub("&","",gsub(" ","",gsub(" - ","_",clin_dis$`Disease Type`)))

vaf_dis <- merge(vaf,clin_dis,by.x = "Primary_WES",by.y = "WES", all.x = T)
vaf_dis <- vaf_dis %>%
  relocate(`Disease Type`, .after = Metastatic_Sample) %>%
  rename("DiseaseType_Primary" = `Disease Type`)
vaf_dis <- merge(vaf_dis,clin_dis,by.x = "Metastatic_WES",by.y = "WES", all.x = T)
vaf_dis <- vaf_dis %>%
  relocate(`Disease Type`, .after = DiseaseType_Primary) %>%
  rename("DiseaseType_Metastatic" = `Disease Type`)
vaf_dis$Matching_Primary_Met_DiseaseType <- ifelse(vaf_dis$DiseaseType_Primary == vaf_dis$DiseaseType_Metastatic,TRUE,FALSE)
vaf_dis <- vaf_dis %>%
  relocate(Matching_Primary_Met_DiseaseType, .after = DiseaseType_Metastatic)

clin_dis_pat <- clin %>%
  filter(!is.na(AvatarKey)) %>%
  select(AvatarKey,`Disease Type`) %>%
  group_by(AvatarKey) %>%
  mutate(`Disease Type` = paste(unique(`Disease Type`), collapse = ',')) %>%
  unique() %>%
  as.data.frame()
clin_dis_pat$`Disease Type` <- gsub("&","",gsub(" ","",gsub(" - ","_",clin_dis_pat$`Disease Type`)))

vaf_dis <- merge(vaf_dis,clin_dis_pat,by.x = "ORIENAvatarKey",by.y = "AvatarKey", all.x = T)
vaf_dis <- vaf_dis %>%
  relocate(`Disease Type`, .after = Metastatic_Sample) %>%
  rename("DiseaseType" = `Disease Type`) %>%
  as.data.frame()

# Append Patient Age -----------------------------------------------------------

clin$`Age At Specimen Collection` <- ifelse(clin$`Age At Specimen Collection` == "Age 90 or older", 91,clin$`Age At Specimen Collection`)

clin_age <- clin %>%
  filter(!is.na(WES)) %>%
  select(WES,`Age At Specimen Collection`) %>%
  unique() %>%
  as.data.frame()

vaf_dis_age <- merge(vaf_dis,clin_age,by.x = "Primary_WES",by.y = "WES", all.x = T)
vaf_dis_age$Primary_Sample <- paste0(vaf_dis_age$Primary_Sample,"_Age",vaf_dis_age$`Age At Specimen Collection`,"Years")
vaf_dis_age <- vaf_dis_age %>%
  select(!`Age At Specimen Collection`)
vaf_dis_age <- merge(vaf_dis_age,clin_age,by.x = "Metastatic_WES",by.y = "WES", all.x = T)
vaf_dis_age$Metastatic_Sample <- paste0(vaf_dis_age$Metastatic_Sample,"_Age",vaf_dis_age$`Age At Specimen Collection`,"Years")
vaf_dis_age <- vaf_dis_age %>%
  select(!`Age At Specimen Collection`)



# Append ICI data --------------------------------------------------------------

## subset clinical data to patients of interest
clin_vaf <- clin[which(clin$AvatarKey %in% vaf_dis_age$ORIENAvatarKey),]

## subset ICI med info
ici_drugs <- c("Ipilimumab", "Atezolizumab", "Avelumab", "Durvalumab", "Nivolumab", "Pembrolizumab")
ici_drug_cols <- grep(paste(ici_drugs,collapse = "|"),colnames(clin_vaf))[c(15:20)]
clin_vaf2 <- clin_vaf[,c(1,2,3,10,11,17,ici_drug_cols)]

## Generate T/F column of ICI treatment and age at first ICI treatment
clin_vaf2$ICI_Treated <- as.logical(rowSums(clin_vaf2[,c(7:12)],na.rm = T))
my_min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
clin_vaf2$Age_At_First_ICI_Treatment <- apply(clin_vaf2[,c(7:12)],1,my_min)

## separate primary and met samples to look at specifc tissue source and ici data
clin_vaf2_prim <- clin_vaf2[which(clin_vaf2$`Primary/Met` == "Primary"),]
colnames(clin_vaf2_prim)[c(3,4,6)] <- c("Primary_WES","SpecimenSitOFCollection_Primary","AgeAtSpecimenCollection_Primary")
clin_vaf2_met <- clin_vaf2[which(clin_vaf2$`Primary/Met` == "Metastatic"),]
colnames(clin_vaf2_met)[c(3,4,6)] <- c("Metastatic_WES","SpecimenSitOFCollection_Metastatic","AgeAtSpecimenCollection_Metastatic")

## Merge primary sample data to make column for sample collection source
clin_vaf2_prim$ICI_BeforeOrAfter_Primary_Sample_Collection <- ifelse(clin_vaf2_prim$AgeAtSpecimenCollection_Primary<clin_vaf2_prim$Age_At_First_ICI_Treatment,"After","Before")
vaf_dis_age_ici <- merge(vaf_dis_age,clin_vaf2_prim[c(3,4,15)],all.x = T, sort = F)

## Annotate if ICI treatment was before or after metastatic sample collection
clin_vaf2_met$ICI_BeforeOrAfter_Met_Sample_Collection <- ifelse(clin_vaf2_met$AgeAtSpecimenCollection_Metastatic<clin_vaf2_met$Age_At_First_ICI_Treatment,"After","Before")
vaf_dis_age_ici <- merge(vaf_dis_age_ici,clin_vaf2_met[c(3,4,7:15)],all.x = T, sort = F)
ICI_toLogical <- function(x) ifelse(is.na(x),FALSE,TRUE)

ICI_Drug_Cols_inVAF <- grep(paste(ici_drugs, collapse = "|"),colnames(vaf_dis_age_ici))
vaf_dis_age_ici[,ICI_Drug_Cols_inVAF] <- apply(vaf_dis_age_ici[,ICI_Drug_Cols_inVAF], 2, ICI_toLogical)
colnames(vaf_dis_age_ici)[ICI_Drug_Cols_inVAF] <- gsub("StartAge","Treated",colnames(vaf_dis_age_ici)[ICI_Drug_Cols_inVAF])


vaf_dis_age_ici <- vaf_dis_age_ici %>%
  relocate(ORIENAvatarKey,Primary_WES,Metastatic_WES,
           Primary_Sample,Metastatic_Sample,
           DiseaseType,DiseaseType_Primary,DiseaseType_Metastatic,Matching_Primary_Met_DiseaseType,
           SpecimenSitOFCollection_Primary,SpecimenSitOFCollection_Metastatic,
           Primary_Sample_Exome_Kit, Metastatic_Sample_Exome_Kit, Matching_Primary_Met_Exome_Kit,
           ICI_Treated,Age_At_First_ICI_Treatment,ICI_BeforeOrAfter_Primary_Sample_Collection,ICI_BeforeOrAfter_Met_Sample_Collection,
           Atezolizumab_Treated,Avelumab_Treated,Durvalumab_Treated,Ipilimumab_Treated,Nivolumab_Treated,Pembrolizumab_Treated)


# Calculate Percent InFrame and Out of Frame variants --------------------------

vaf_dis_age_ici[,paste0("Percent_InFrame_Primary_Variant")] <- vaf_dis_age_ici$PrimaryMet_Shared_Variant_INFRAME/(vaf_dis_age_ici$Primary_Exclusive_Variant_INFRAME+vaf_dis_age_ici$PrimaryMet_Shared_Variant_INFRAME)
vaf_dis_age_ici[,paste0("Percent_OutOfFrame_Primary_Variant")] <- vaf_dis_age_ici$PrimaryMet_Shared_Variant_OUTOFFRAME/(vaf_dis_age_ici$Primary_Exclusive_Variant_OUTOFFRAME+vaf_dis_age_ici$PrimaryMet_Shared_Variant_OUTOFFRAME)

vaf_dis_age_ici[,paste0("Percent_InFrame_Metastatic_Variant")] <- vaf_dis_age_ici$PrimaryMet_Shared_Variant_INFRAME/(vaf_dis_age_ici$Metastatic_Variant_INFRAME+vaf_dis_age_ici$PrimaryMet_Shared_Variant_INFRAME)
vaf_dis_age_ici[,paste0("Percent_OutOfFrame_Metastatic_Variant")] <- vaf_dis_age_ici$PrimaryMet_Shared_Variant_OUTOFFRAME/(vaf_dis_age_ici$Metastatic_Variant_OUTOFFRAME+vaf_dis_age_ici$PrimaryMet_Shared_Variant_OUTOFFRAME)



# Add Site of collection major category ----------------------------------------

site_prim <- site[which(site$X1 %in% vaf_dis_age_ici$SpecimenSitOFCollection_Primary),]
colnames(site_prim) <- c("SpecimenSitOFCollection_Primary","SpecimenSitOFCollection_Primary_Major")
site_met <- site[which(site$X1 %in% vaf_dis_age_ici$SpecimenSitOFCollection_Metastatic),]
colnames(site_met) <- c("SpecimenSitOFCollection_Metastatic","SpecimenSitOFCollection_Metastatic_Major")

vaf_dis_age_ici_site <- merge(vaf_dis_age_ici,site_prim,all.x = T)
vaf_dis_age_ici_site <- merge(vaf_dis_age_ici_site,site_met,all.x = T)
df_check <- vaf_dis_age_ici_site %>% select(SpecimenSitOFCollection_Metastatic,SpecimenSitOFCollection_Primary,SpecimenSitOFCollection_Primary_Major,SpecimenSitOFCollection_Metastatic_Major)

## Manual Check of MHC1hit and MHC2hit 01192024 Specimen sites -----------------
vaf_dis_age_ici_site[274,"SpecimenSitOFCollection_Primary_Major"] <- "Other"
vaf_dis_age_ici_site[69,"SpecimenSitOFCollection_Metastatic_Major"] <- "Liver"
vaf_dis_age_ici_site[267,"SpecimenSitOFCollection_Metastatic_Major"] <- "Lymph node"
vaf_dis_age_ici_site[442,"SpecimenSitOFCollection_Metastatic_Major"] <- "Other"
vaf_dis_age_ici_site[455,"SpecimenSitOFCollection_Metastatic_Major"] <- "Other"


vaf_dis_age_ici_site$Matching_Primary_Met_SpecimenSitOFCollection <- ifelse(vaf_dis_age_ici_site$SpecimenSitOFCollection_Primary == vaf_dis_age_ici_site$SpecimenSitOFCollection_Metastatic,
                                                                            TRUE,FALSE)
vaf_dis_age_ici_site$Matching_Primary_Met_SpecimenSitOFCollection_Major <- ifelse(vaf_dis_age_ici_site$SpecimenSitOFCollection_Primary_Major == vaf_dis_age_ici_site$SpecimenSitOFCollection_Metastatic_Major,
                                                                            TRUE,FALSE)

vaf_dis_age_ici_site <- vaf_dis_age_ici_site %>%
  relocate(ORIENAvatarKey,Primary_WES,Metastatic_WES,
           Primary_Sample,Metastatic_Sample,
           DiseaseType,DiseaseType_Primary,DiseaseType_Metastatic,Matching_Primary_Met_DiseaseType,
           SpecimenSitOFCollection_Primary_Major,SpecimenSitOFCollection_Primary,
           SpecimenSitOFCollection_Metastatic_Major,SpecimenSitOFCollection_Metastatic,
           Matching_Primary_Met_SpecimenSitOFCollection_Major,Matching_Primary_Met_SpecimenSitOFCollection)


# Add tumor content ------------------------------------------------------------

tc <- clin_vaf[,c(3,16)]
tc[,2] <- gsub("> 30","31", tc[,2])
tc[,2] <- as.numeric(tc[,2])

vaf_dis_age_ici_site_tc <- merge(vaf_dis_age_ici_site,tc, by.x = "Primary_WES", by.y = "WES",all.x = T)

vaf_dis_age_ici_site_tc <- vaf_dis_age_ici_site_tc %>%
  relocate(`% Tumor Content`, .after = Matching_Primary_Met_DiseaseType) %>%
  rename("Percent_Tumor_Content_Primary" = `% Tumor Content`)
vaf_dis_age_ici_site_tc <- merge(vaf_dis_age_ici_site_tc,tc,by.x = "Metastatic_WES",by.y = "WES", all.x = T)
vaf_dis_age_ici_site_tc <- vaf_dis_age_ici_site_tc %>%
  relocate(`% Tumor Content`, .after = Percent_Tumor_Content_Primary) %>%
  rename("Percent_Tumor_Content_Metastatic" = `% Tumor Content`)

vaf_dis_age_ici_site_tc <- vaf_dis_age_ici_site_tc %>%
  mutate(Met_Minus_Prim_Tumor_Content = Percent_Tumor_Content_Metastatic-Percent_Tumor_Content_Primary) %>%
  relocate(Met_Minus_Prim_Tumor_Content, .after = Percent_Tumor_Content_Metastatic)


# Add RNAseq ID to table -------------------------------------------------------

nm <- clin_vaf[,c(2,3)]

vaf_dis_age_ici_site_tc_rna <- merge(nm,vaf_dis_age_ici_site_tc,by.x = "WES",by.y = "Primary_WES", all.y = T)
colnames(vaf_dis_age_ici_site_tc_rna)[c(1,2)] <- c("Primary_WES","Primary_RNASeq")

vaf_dis_age_ici_site_tc_rna <- merge(nm,vaf_dis_age_ici_site_tc_rna,by.x = "WES",by.y = "Metastatic_WES", all.y = T)
colnames(vaf_dis_age_ici_site_tc_rna)[c(1,2)] <- c("Metastatic_WES","Metastatic_RNASeq")

vaf_dis_age_ici_site_tc_rna <- vaf_dis_age_ici_site_tc_rna %>%
  relocate(ORIENAvatarKey,Primary_RNASeq,Primary_WES,Metastatic_RNASeq,Metastatic_WES)

write.table(vaf_dis_age_ici_site_tc_rna,Output_File,sep = '\t', col.names = T)



# Add ER Stress ssGSEA ---------------------------------------------------------

ER_Stress_prim <- ER_Stress[which(ER_Stress$RNASeq %in% vaf_dis_age_ici_site_tc_rna$Primary_RNASeq),]
colnames(ER_Stress_prim)[-1] <- paste0(colnames(ER_Stress_prim)[-1],"_Prim")
colnames(ER_Stress_prim)[1] <- "Primary_RNASeq"
ER_Stress_met <- ER_Stress[which(ER_Stress$RNASeq %in% vaf_dis_age_ici_site_tc_rna$Metastatic_RNASeq),]
colnames(ER_Stress_met)[-1] <- paste0(colnames(ER_Stress_met)[-1],"_Met")
colnames(ER_Stress_met)[1] <- "Metastatic_RNASeq"

vaf_dis_age_ici_site_tc_rna_erP <- merge(vaf_dis_age_ici_site_tc_rna,ER_Stress_prim, by = "Primary_RNASeq", all.x = T)
vaf_dis_age_ici_site_tc_rna_erP_erM <- merge(vaf_dis_age_ici_site_tc_rna_erP,ER_Stress_met, by = "Metastatic_RNASeq", all.x = T)


ER_Stress_prim <- as.matrix(vaf_dis_age_ici_site_tc_rna_erP_erM[,grep("_Prim$",colnames(vaf_dis_age_ici_site_tc_rna_erP_erM),value = T)])
ER_Stress_met <- as.matrix(vaf_dis_age_ici_site_tc_rna_erP_erM[,grep("_Met$",colnames(vaf_dis_age_ici_site_tc_rna_erP_erM),value = T)])

ER_Stress_Sub <- ER_Stress_met - ER_Stress_prim
colnames(ER_Stress_Sub) <- gsub("_Met$","_MetMinusPrim",colnames(ER_Stress_Sub))
ER_Stress_Sub <- as.data.frame(ER_Stress_Sub)

vaf_dis_age_ici_site_tc_rna_erP_erM_erSub <- cbind(vaf_dis_age_ici_site_tc_rna_erP_erM,ER_Stress_Sub)

vaf_dis_age_ici_site_tc_rna_erP_erM_erSub <- vaf_dis_age_ici_site_tc_rna_erP_erM_erSub %>%
  relocate(ORIENAvatarKey,Primary_RNASeq,Primary_WES,Metastatic_RNASeq,Metastatic_WES)

write.table(vaf_dis_age_ici_site_tc_rna_erP_erM_erSub,gsub(paste0("\\.",tools::file_ext(Output_File)),"_ERStress.txt",Output_File),
            sep = '\t', row.names = F)



