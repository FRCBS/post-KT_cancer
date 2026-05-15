
library(tidyverse)
library(data.table)
library(corrplot)
library(knitr)
library(kableExtra)
library(ggsci)
library(ggpubr)
library(patchwork)
library(readxl)


## ========================================================
##
## FinnGen tx cancer cohort summary stats
##
## ========================================================


## --------------------------------------------------------
## data read-in
## --------------------------------------------------------

# longitudinal event data
cancer.kidney.first    <- fread('data/phenotypes/cancer.kidney.first')
cancer.kidney.first.tx <- fread('data/phenotypes/cancer.kidney.first.tx')

# Pancancer SNP risk scores
risk.scores <- rbind(
  fread('data/FG_pancancer_scores_01.tsv'),
  fread('data/FG_pancancer_scores_02.tsv')
)

# primary diag
tx.longi.primdiag <- fread('data/phenotypes/TX_LONG_phe_primary_2.tsv')



## --------------------------------------------------------
## cohort info
## --------------------------------------------------------

# KT age distribution
cancer.kidney.first.tx$EVENT_AGE.kidney %>% hist(20)


## --------------------------------------------------------
## cohort data
## --------------------------------------------------------

# data sets
dall <- fread('data/phenotypes/TX_LONG_phe_primary_counts.tsv')
dall <- cancer.kidney.first.tx[, c(1:6, 149:151, 165:168)] # all

# get primary diagnoses for each patient
dall$Diabetes <- map(dall$FINNGENID, function(x) {
  filter(tx.longi.primdiag, FINNGENID==x, Primary=='Diabetes') %>% nrow
}) %>% unlist
dall$RTID <- map(dall$FINNGENID, function(x) {
  filter(tx.longi.primdiag, FINNGENID==x, Primary=='Renal tubulo-interstitial disease') %>% nrow
}) %>% unlist
dall$Hypertension <- map(dall$FINNGENID, function(x) {
  filter(tx.longi.primdiag, FINNGENID==x, Primary=='Hypertension') %>% nrow
}) %>% unlist
dall$CKD <- map(dall$FINNGENID, function(x) {
  filter(tx.longi.primdiag, FINNGENID==x, Primary=='Chronic kidney disease') %>% nrow
}) %>% unlist
dall$Glomerulonephritis <- map(dall$FINNGENID, function(x) {
  filter(tx.longi.primdiag, FINNGENID==x, Primary=='Glomerulonephritis') %>% nrow
}) %>% unlist
dall$Congenital <- map(dall$FINNGENID, function(x) {
  filter(tx.longi.primdiag, FINNGENID==x, Primary=='Congenital malformation of the urinary system') %>% nrow
}) %>% unlist
dall$Cystic <- map(dall$FINNGENID, function(x) {
  filter(tx.longi.primdiag, FINNGENID==x, Primary=='Cystic kidney disease') %>% nrow
}) %>% unlist


# get the fist primary diagnosis for each patient
dall$FirstPrim <- map(dall$FINNGENID, function(x) {
  res <- filter(tx.longi.primdiag, FINNGENID==x)
  dis <- c('Diabetes', 'Renal tubulo-interstitial disease', 
           'Hypertension', 'Chronic kidney disease', 'Glomerulonephritis',
           'Congenital malformation of the urinary system', 'Cystic kidney disease')
  res <- c(which(res$Primary==dis[1])[1],
           which(res$Primary==dis[2])[1],
           which(res$Primary==dis[3])[1],
           which(res$Primary==dis[4])[1],
           which(res$Primary==dis[5])[1],
           which(res$Primary==dis[6])[1],
           which(res$Primary==dis[7])[1])
  ind <- which.min(res)
  if(length(ind)>0) dis[ind] %>% return() else return('Other')
  
}) %>% unlist


# write out primary counts per patient
fwrite(dall, 'data/phenotypes/TX_LONG_phe_primary_counts.tsv', sep='\t') 
dall <- fread('data/phenotypes/TX_LONG_phe_primary_counts.tsv', data.table = F)

# subset by cancer status
dcan <- dall %>% filter(EVENT==1) # with cancer
dnoc <- dall %>% filter(EVENT==0) # without cancer


## --------------------------------------------------------
## summary stat table info
## --------------------------------------------------------

# patient numbers
c(
  dall %>% nrow, dcan %>% nrow, dnoc %>% nrow
)

# male patients 
c(
  dall %>% filter(SEX=='male') %>% nrow,
  (dall %>% filter(SEX=='male') %>% nrow) / (dall %>% nrow),
  dcan %>% filter(SEX=='male') %>% nrow,
  (dcan %>% filter(SEX=='male') %>% nrow) / (dcan %>% nrow),
  dnoc %>% filter(SEX=='male') %>% nrow,
  (dnoc %>% filter(SEX=='male') %>% nrow) / (dnoc %>% nrow),
  prop.test(c(dcan %>% filter(SEX=='male') %>% nrow, dnoc %>% filter(SEX=='male') %>% nrow), 
            c(dcan %>% nrow,                         dnoc %>% nrow))$p.value
)

# tx ages
c(
  dall %>% .$EVENT_AGE.kidney %>% median,
  dall %>% .$EVENT_AGE.kidney %>% IQR,
  dcan %>% .$EVENT_AGE.kidney %>% median,
  dcan %>% .$EVENT_AGE.kidney %>% IQR,
  dnoc %>% .$EVENT_AGE.kidney %>% median,
  dnoc %>% .$EVENT_AGE.kidney %>% IQR,
  t.test(dcan %>% .$EVENT_AGE.kidney, dnoc %>% .$EVENT_AGE.kidney)$p.value
)

# time from tx
c(
  dall %>% .$TIME_TX_FU %>% median,
  dall %>% .$TIME_TX_FU %>% IQR,
  dcan %>% .$TIME_TX_FU %>% median,
  dcan %>% .$TIME_TX_FU %>% IQR,
  dnoc %>% .$TIME_TX_FU %>% median,
  dnoc %>% .$TIME_TX_FU %>% IQR,
  t.test(dcan %>% .$TIME_TX_FU, dnoc %>% .$TIME_TX_FU)$p.value
)

# # diabetes primary diag  
# c(
#   dall %>% filter(Diabetes>0) %>% nrow,
#   (dall %>% filter(Diabetes>0) %>% nrow) / (dall %>% nrow),
#   dcan %>% filter(Diabetes>0) %>% nrow,
#   (dcan %>% filter(Diabetes>0) %>% nrow) / (dcan %>% nrow),
#   dnoc %>% filter(Diabetes>0) %>% nrow,
#   (dnoc %>% filter(Diabetes>0) %>% nrow) / (dnoc %>% nrow),
#   prop.test(c(dcan %>% filter(Diabetes>0) %>% nrow, dnoc %>% filter(Diabetes>0) %>% nrow), 
#             c(dcan %>% nrow,                         dnoc %>% nrow))$p.value
# ) %>% round(4)
# 
# # glomerulonephritis primary diag  
# c(
#   dall %>% filter(Glomerulonephritis>0) %>% nrow,
#   (dall %>% filter(Glomerulonephritis>0) %>% nrow) / (dall %>% nrow),
#   dcan %>% filter(Glomerulonephritis>0) %>% nrow,
#   (dcan %>% filter(Glomerulonephritis>0) %>% nrow) / (dcan %>% nrow),
#   dnoc %>% filter(Glomerulonephritis>0) %>% nrow,
#   (dnoc %>% filter(Glomerulonephritis>0) %>% nrow) / (dnoc %>% nrow),
#   prop.test(c(dcan %>% filter(Glomerulonephritis>0) %>% nrow, dnoc %>% filter(Glomerulonephritis>0) %>% nrow), 
#             c(dcan %>% nrow,                         dnoc %>% nrow))$p.value
# ) %>% round(4)
# 
# # chronic KD primary diag  
# c(
#   dall %>% filter(CKD>0) %>% nrow,
#   (dall %>% filter(CKD>0) %>% nrow) / (dall %>% nrow),
#   dcan %>% filter(CKD>0) %>% nrow,
#   (dcan %>% filter(CKD>0) %>% nrow) / (dcan %>% nrow),
#   dnoc %>% filter(CKD>0) %>% nrow,
#   (dnoc %>% filter(CKD>0) %>% nrow) / (dnoc %>% nrow),
#   prop.test(c(dcan %>% filter(CKD>0) %>% nrow, dnoc %>% filter(CKD>0) %>% nrow), 
#             c(dcan %>% nrow,                         dnoc %>% nrow))$p.value
# ) %>% round(4)
# 
# # Congenital primary diag  
# c(
#   dall %>% filter(Congenital>0) %>% nrow,
#   (dall %>% filter(Congenital>0) %>% nrow) / (dall %>% nrow),
#   dcan %>% filter(Congenital>0) %>% nrow,
#   (dcan %>% filter(Congenital>0) %>% nrow) / (dcan %>% nrow),
#   dnoc %>% filter(Congenital>0) %>% nrow,
#   (dnoc %>% filter(Congenital>0) %>% nrow) / (dnoc %>% nrow),
#   prop.test(c(dcan %>% filter(Congenital>0) %>% nrow, dnoc %>% filter(Congenital>0) %>% nrow), 
#             c(dcan %>% nrow,                         dnoc %>% nrow))$p.value
# ) %>% round(4)
# 
# # Hypertension  
# c(
#   dall %>% filter(Hypertension>0) %>% nrow,
#   (dall %>% filter(Hypertension>0) %>% nrow) / (dall %>% nrow),
#   dcan %>% filter(Hypertension>0) %>% nrow,
#   (dcan %>% filter(Hypertension>0) %>% nrow) / (dcan %>% nrow),
#   dnoc %>% filter(Hypertension>0) %>% nrow,
#   (dnoc %>% filter(Hypertension>0) %>% nrow) / (dnoc %>% nrow),
#   prop.test(c(dcan %>% filter(Hypertension>0) %>% nrow, dnoc %>% filter(Hypertension>0) %>% nrow), 
#             c(dcan %>% nrow,                         dnoc %>% nrow))$p.value
# ) %>% round(4)
# 
# # Renal tubulo-interstitial disease  
# c(
#   dall %>% filter(RTID>0) %>% nrow,
#   (dall %>% filter(RTID>0) %>% nrow) / (dall %>% nrow),
#   dcan %>% filter(RTID>0) %>% nrow,
#   (dcan %>% filter(RTID>0) %>% nrow) / (dcan %>% nrow),
#   dnoc %>% filter(RTID>0) %>% nrow,
#   (dnoc %>% filter(RTID>0) %>% nrow) / (dnoc %>% nrow),
#   prop.test(c(dcan %>% filter(RTID>0) %>% nrow, dnoc %>% filter(RTID>0) %>% nrow), 
#             c(dcan %>% nrow,                         dnoc %>% nrow))$p.value
# ) %>% round(4)
# 
# # Cystic kidney disease  
# c(
#   dall %>% filter(Cystic>0) %>% nrow,
#   (dall %>% filter(Cystic>0) %>% nrow) / (dall %>% nrow),
#   dcan %>% filter(Cystic>0) %>% nrow,
#   (dcan %>% filter(Cystic>0) %>% nrow) / (dcan %>% nrow),
#   dnoc %>% filter(Cystic>0) %>% nrow,
#   (dnoc %>% filter(Cystic>0) %>% nrow) / (dnoc %>% nrow),
#   prop.test(c(dcan %>% filter(Cystic>0) %>% nrow, dnoc %>% filter(Cystic>0) %>% nrow), 
#             c(dcan %>% nrow,                         dnoc %>% nrow))$p.value
# ) %>% round(4)
# 
# 
# ## first occurrencies
# 
# statFirstDiag <- function(diag, rounding=4) {
#   c(
#     dall %>% filter(FirstPrim==diag) %>% nrow,
#     (dall %>% filter(FirstPrim==diag) %>% nrow) / (dall %>% nrow),
#     dcan %>% filter(FirstPrim==diag) %>% nrow,
#     (dcan %>% filter(FirstPrim==diag) %>% nrow) / (dcan %>% nrow),
#     dnoc %>% filter(FirstPrim==diag) %>% nrow,
#     (dnoc %>% filter(FirstPrim==diag) %>% nrow) / (dnoc %>% nrow),
#     prop.test(c(dcan %>% filter(FirstPrim==diag) %>% nrow, 
#                 dnoc %>% filter(FirstPrim==diag) %>% nrow), 
#               c(dcan %>% nrow,
#                 dnoc %>% nrow))$p.value
#   ) %>% round(rounding) %>% return()
# }
# 
# c('Diabetes', 'Glomerulonephritis', 'Chronic kidney disease', 
#   'Congenital malformation of the urinary system', 'Hypertension',  
#   'Renal tubulo-interstitial disease', 'Cystic kidney disease', 'Other') %>% 
#   map(function(x) statFirstDiag(x))
# 
# statFirstDiag('Diabetes', 10)
# 



## --------------------------------------------------------
## Kidney disease register
## --------------------------------------------------------

cancer.kidney.first.tx$FINNGENID %>% unique %>% length

# "KIDNEYREG_T1DM"	Type 1 diabetes with a causal role to some kidney disease. E10
# "KIDNEYREG_T2DM"	Type 2 diabetes with a causal role to some kidney disease. E11
# "KIDNEYREG_POLYCYS"	Polycystic kidney disease. Q61[1-3] 
# "KIDNEYREG_GLOMNEPH"	Glomerulonephritis. N0[0|1|3]
# "KIDNEYREG_NEPHROSCLER"	Nephrosclerosis. I12|I13|I701|N280
# "KIDNEYREG_TUBULOINT"	Tubulointerstitial nephritis. N1[0-2]
# "KIDNEYREG_AMYLOIDOSIS"	Amyloidosis. E85
# kd <- fread2('/finngen/library-red/finngen_R12/kidney_disease_register_1.0/data/finngen_R12_kidney_disease_register_1.0.txt') %>% 
#   filter(FINNGENID %in% cancer.kidney.first.tx$FINNGENID) %>% 
#   dplyr::select(-KIDNEYREG_ACTIVE)
# kd %>% dim
# kd <- left_join(cancer.kidney.first.tx[, 1:2], kd)
# kd %>% dim
# kd[, -c(1:3)] %>% colSums(na.rm = T)

##
kd <- fread2('/finngen/library-red/finngen_R12/kidney_disease_register_1.0/data/finngen_R12_kidney_combined_1.0.txt') %>% 
  filter(FINNGENID %in% cancer.kidney.first.tx$FINNGENID)
kd %>% dim
kd$FINNGENID %>% unique %>% length
kd.diag <- map(kd$FINNGENID %>% unique, function(x) {
  tmp <- filter(kd, FINNGENID==x)
  data.frame(FINNGENID=x, 
             KIDNEY_DISEASE_DIAGNOSIS=tmp$KIDNEY_DISEASE_DIAGNOSIS_1[1])
}) %>% do.call(rbind, .)
kd.diag %>% dim
kd.diag <- data.frame(kd.diag,
                      kd.diag$KIDNEY_DISEASE_DIAGNOSIS %>% str_split_fixed('\\*', 2))
kd.diag$X1 %>% table
kd.diag <- left_join(cancer.kidney.first.tx %>% dplyr::select(c(FINNGENID, EVENT)), 
                     kd.diag)
icd <- read_xlsx('data/ICD-10_MIT_2021_Excel_16-March_2021.xlsx')
icd$ICD_newcode <- icd$ICD10_Code %>% gsub('\\.', '', .)
kd.diag <- left_join(kd.diag, 
                     icd %>% dplyr::select(Group_Desc, ICD10_3_Code_Desc, ICD10_Code, WHO_Full_Desc, ICD_newcode),
                     by = c('X1'='ICD_newcode'))
kd.diag %>% head

kd.diag.can <- kd.diag %>% filter(EVENT==1) %>% .$ICD10_3_Code_Desc %>% 
  table(useNA = 'always') %>% sort %>% data.frame
kd.diag.noc <- kd.diag %>% filter(EVENT==0) %>% .$ICD10_3_Code_Desc %>% 
  table(useNA = 'always') %>% sort %>% data.frame
colnames(kd.diag.can)[1] <- colnames(kd.diag.noc)[1] <- 'Diag'
kd.diag.can[, 1] <- as.character(kd.diag.can[, 1])
kd.diag.noc[, 1] <- as.character(kd.diag.noc[, 1])

kd.diag.can$Diag[
  is.na(kd.diag.can$Diag)] <- 'Unknown'
kd.diag.can$Diag[
  kd.diag.can$Diag == 'Chronic kidney disease'] <- 'Unknown'
kd.diag.can$Diag[
  kd.diag.can$Diag=='Glomerular disorders in diseases classified elsewhere'] <-
  'Glomerular disorders in diabetes mellitus'

kd.diag.noc$Diag[
  is.na(kd.diag.noc$Diag)] <- 'Unknown'
kd.diag.noc$Diag[
  kd.diag.noc$Diag == 'Chronic kidney disease'] <- 'Unknown'
kd.diag.noc$Diag[
  kd.diag.noc$Diag=='Glomerular disorders in diseases classified elsewhere'] <-
  'Glomerular disorders in diabetes mellitus'

cat.other <- c(filter(kd.diag.can, Freq<5)$Diag, 
               filter(kd.diag.noc, Freq<5)$Diag,
               'Rapidly progressive nephritic syndrome',
               'Hypertensive renal disease') %>%  unique
glonef <- c(cat.other[c(1,13,22,23,26)], 
            'Chronic nephritic syndrome', 'Nephrotic syndrome')
caku <- cat.other[c(4,7,8,12,14,19)]
cat.other <- cat.other[!(cat.other %in% c(glonef, caku))]

kd.diag.can$Diag[kd.diag.can$Diag %in% cat.other] <- 'Other'
kd.diag.noc$Diag[kd.diag.noc$Diag %in% cat.other] <- 'Other'

kd.diag.can$Diag[kd.diag.can$Diag %in% glonef] <- 'Glomerulonephritis'
kd.diag.noc$Diag[kd.diag.noc$Diag %in% glonef] <- 'Glomerulonephritis'
kd.diag.noc$Diag[kd.diag.noc$Diag %in% caku] <- 'Congenital malformations'
kd.diag.can$Diag[kd.diag.can$Diag %in% caku] <- 'Congenital malformations'

kd.diag.can <- kd.diag.can %>% group_by(Diag) %>% summarise(SUM=sum(Freq)) 
kd.diag.noc <- kd.diag.noc %>% group_by(Diag) %>% summarise(SUM=sum(Freq)) 
kd.diag.can

kd.diags <- full_join(kd.diag.can, kd.diag.noc, by = 'Diag')
kd.diags[is.na(kd.diags)] <- 0
colnames(kd.diags)[2:3] <- c('Cancer', 'NoCancer')
kd.diags[,2:3] %>% colSums
kd.diags$CancerPer <- ((kd.diags[, 2] / (kd.diags[, 2] %>% colSums))*100) %>% signif(2)
kd.diags$NoCancerPer <- ((kd.diags[, 3] / (kd.diags[, 3] %>% colSums))*100) %>% signif(2)
kd.diags$All <- kd.diags$Cancer + kd.diags$NoCancer
kd.diags$AllPer <- ((kd.diags$All / (kd.diags$All %>% sum))*100) %>% signif(2)
kd.diags$pvalue <- 
  map(1:nrow(kd.diags), function(i) {
    prop.test(x = kd.diags[i, 2:3] %>% unlist, 
              n = kd.diags[,2:3] %>% colSums)$p.value
  }) %>% unlist %>% signif(3)


data.frame(
  Diagnosis = kd.diags$Diag,
  All = paste0(kd.diags$All %>% unlist, ' (', kd.diags$AllPer %>% unlist, ')'),
  Cancer = paste0(kd.diags$Cancer %>% unlist, ' (', kd.diags$CancerPer %>% unlist, ')'),
  NoCancer = paste0(kd.diags$NoCancer %>% unlist, ' (', kd.diags$NoCancerPer %>% unlist, ')'),
  Pvalue = kd.diags$pvalue) %>% 
  fwrite('results/kidney_cancer/KT_groups_diagnoses.tsv', sep = '\t')






# ## --------------------------------------------------------
# ## Summary of TX patients' cancer diagnoses
# ## --------------------------------------------------------
# 
# # cancers in tx patients
# phe <- fread('data/phenotypes/TX_EVENT_cancers.tsv', data.table=F)
# 
# # counts of different cancer types
# phe <- data.frame(Endpoint = phe[, grepl('CD2_|C3_', colnames(phe))] %>% colnames,
#                   Count = phe[, grepl('CD2_|C3_', colnames(phe))] %>% colSums(na.rm=T)) %>% 
#   filter(Count>4, !(Endpoint %in% c('CD2_NEOPLASM', 'C3_CANCER_WIDE', 'C3_CANCER', 'C3_SKIN', 
#                                     'C3_OTHER_SKIN', 'CD2_INSITU'))) %>% 
#   arrange(Count)
# phe$Endpoint <- factor(phe$Endpoint, levels = phe$Endpoint)
# 
# p.tx.diag <- ggplot(phe, aes(Count, Endpoint)) +
#   geom_bar(stat='identity') +
#   theme_minimal() +
#   theme(axis.text.y=element_text(size=7))
# p.tx.diag
# 
# jpeg('results/kidney_cancer/Tx_diagnoses.jpg', width=5, height=7, res=600, units='in')
# p.tx.diag
# dev.off()
# 
