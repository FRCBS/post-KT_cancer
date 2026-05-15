
library(tidyverse)
library(dplyr)
library(data.table)
library(R.utils)
library(corrplot)


## ========================================================
##
## Extract tx cancer longitudinal pheno data
##
## ========================================================


## --------------------------------------------------------
## longitudinal phenotype data
## --------------------------------------------------------

# get kidney tx longitudinal endpoints
system(
  paste0("zcat /finngen/library-red/finngen_R12/phenotype_1.0/data/finngen_R12_endpoint_longitudinal_1.0.txt.gz ",
         "| grep KIDNEY_TRANSPLANT > data/phenotypes/R12_longitudinal_1.0_KidneyTx.txt") 
)
kidney <- fread("data/phenotypes/R12_longitudinal_1.0_KidneyTx.txt", data.table = F)
colnames(kidney) <- c('FINNGENID',	'EVENT_TYPE',	'EVENT_AGE',	'EVENT_YEAR',	'ICDVER',	'ENDPOINT')
kidney$EVENT_YEAR %>% min() # "1995-12-20"
kidney$EVENT_YEAR %>% max() # "2023-04-15"


# get liver tx IDs
system(
  paste0('zcat /finngen/library-red/finngen_R13/phenotype_1.0/data/',
         'finngen_R13_detailed_longitudinal_1.0.txt.gz | grep JJC[01]0',
         ' > data/phenotypes/R13_longitudinal_LiverTxIDs.txt'))

# get heart tx IDs
system(
  paste0('zcat /finngen/library-red/finngen_R13/phenotype_1.0/data/',
         'finngen_R13_detailed_longitudinal_1.0.txt.gz | grep -e FQA10 -e FQ00',
         ' > data/phenotypes/R13_longitudinal_HeartTxIDs.txt'))

# get all transplants
system(
  paste0("zcat /finngen/library-red/finngen_R12/phenotype_1.0/data/finngen_R12_endpoint_longitudinal_1.0.txt.gz ",
         "| grep TRANSPLANT > data/phenotypes/R12_longitudinal_1.0_AllTx.txt") 
)
alltx <- fread('data/phenotypes/R12_longitudinal_1.0_AllTx.txt', data.table = F)
colnames(alltx) <- c('FINNGENID',	'EVENT_TYPE',	'EVENT_AGE',	'EVENT_YEAR',	'ICDVER',	'ENDPOINT')
# endpoints in kidney-only tx
alltx %>% filter(ENDPOINT != 'Z21_LUNG_TRANSPLANT_STATUS',
                 ENDPOINT != 'LUNG_TRANSPLANT',
                 ENDPOINT != 'LUNG_TRANSPLANTATION',
                 ENDPOINT != 'Z21_TRANSPLANTED_ORGAN_TISSUE_STATUS') %>% 
  .$ENDPOINT %>% table
# only non-kidneys
alltx %>% filter(!(FINNGENID %in% kidney$FINNGENID)) %>% 
  .$FINNGENID %>% unique %>% length
alltx %>% filter(!(FINNGENID %in% kidney$FINNGENID)) %>% 
  .$ENDPOINT %>% table
alltx %>% filter(ENDPOINT != 'Z21_KIDNEY_TRANSPLANT_STATUS',
                 ENDPOINT != 'Z21_TRANSPLANTED_ORGAN_TISSUE_STATUS') %>% 
  .$FINNGENID %>% unique %>% length
alltx %>% filter(ENDPOINT != 'Z21_KIDNEY_TRANSPLANT_STATUS') %>% 
  .$FINNGENID %>% unique %>% length
# lung tx
lung <- alltx %>% filter(ENDPOINT == 'LUNG_TRANSPLANTATION',
                 ENDPOINT != 'Z21_KIDNEY_TRANSPLANT_STATUS')
lung %>% .$FINNGENID %>% unique %>% length
# other tx
othertx <- alltx %>% filter(ENDPOINT != 'LUNG_TRANSPLANTATION',
                 ENDPOINT != 'LUNG_TRANSPLANT',
                 ENDPOINT != 'Z21_LUNG_TRANSPLANT_STATUS',
                 ENDPOINT != 'Z21_KIDNEY_TRANSPLANT_STATUS')
othertx %>% .$FINNGENID %>% unique %>% length
othertx %>% .$ENDPOINT %>% table

# get cancer incidence years from cancer registry
cancer <- fread("/finngen/library-red/finngen_R12/cancer_detailed_1.0/data/finngen_R12_cancer_detailed_1.0.txt", 
                data.table = F) 

# get immunosuppression longitudinal endpoints
system(
  paste0("zcat /finngen/library-red/finngen_R12/phenotype_1.0/data/finngen_R12_endpoint_longitudinal_1.0.txt.gz ",
         "| grep IMMUNOS > data/phenotypes/R12_longitudinal_1.0_Immunosupp.txt") 
)
immuno <- fread("data/phenotypes/R12_longitudinal_1.0_Immunosupp.txt", data.table = F)
colnames(immuno) <- c('FINNGENID',	'EVENT_TYPE',	'EVENT_AGE',	'EVENT_YEAR',	'ICDVER',	'ENDPOINT')

# first event age of kidney tx
# select for having immunosuppression
kidney.first <- map(kidney$FINNGENID %>% unique, function(x) {
  filter(kidney, FINNGENID==x, FINNGENID %in% immuno$FINNGENID) %>% 
    arrange(EVENT_AGE) %>% .[1, ]
}) %>% do.call(rbind, .)
fwrite(kidney.first, 'results/kidney_cancer/kidney.first.tsv', sep='\t')

# first event age of lung tx
# select for having immunosuppression
lung.first <- map(lung$FINNGENID %>% unique, function(x) {
  filter(lung, FINNGENID==x, FINNGENID %in% immuno$FINNGENID) %>% 
    arrange(EVENT_AGE) %>% .[1, ]
}) %>% do.call(rbind, .)
fwrite(lung.first, 'results/kidney_cancer/lung.first.tsv', sep='\t')

# first event age of other tx
# select for having immunosuppression
other.first <- map(othertx$FINNGENID %>% unique, function(x) {
  filter(othertx, FINNGENID==x, FINNGENID %in% immuno$FINNGENID) %>% 
    arrange(EVENT_AGE) %>% .[1, ]
}) %>% do.call(rbind, .)
fwrite(other.first, 'results/kidney_cancer/othertx.first.tsv', sep='\t')

# first cancer incidence
cancer.first <- map(cancer$FINNGENID %>% unique, function(x) {
  filter(cancer, FINNGENID==x) %>% 
    arrange(EVENT_AGE) %>% .[1, ]
}) %>% do.call(rbind, .)
fwrite(cancer.first, 'results/kidney_cancer/cancer.first.tsv', sep='\t')


## --------------------------------------------------------
## 1st cancer diagnosis
## --------------------------------------------------------

fglongi <- fread('/finngen/library-red/finngen_R12/phenotype_1.0/data/finngen_R12_endpoint_longitudinal_1.0.txt.gz')
fglongi <- filter(fglongi, FINNGENID %in% cancer$FINNGENID) 
fglongi %>% dim
fglongi <- left_join(fglongi, cancer[, c(1,9)], by = 'FINNGENID')
fglongi <- filter(fglongi, EVENT_TYPE == 'CANCER')
fglongi <- fglongi[round(fglongi$EVENT_AGE.x) == round(fglongi$EVENT_AGE.y), ]
fglongi <- fglongi[grepl('CD2_|C3_', fglongi$ENDPOINT), ]
fglongi <- fglongi[!grepl('BENIGN|CD2_NEOPLASM|C3_CANCER|C3_CANCER_WIDE', fglongi$ENDPOINT), ]

fglongi1 <- map(fglongi$FINNGENID %>% unique, function(x) {
  filter(fglongi, FINNGENID==x)[1, ] %>% return()
}) %>% do.call(rbind, .)
fglongi1$ENDPOINT %>% table %>% sort
fglongi1$FINNGENID %>% unique %>% length
cancer$FINNGENID %>% unique %>% length
fglongi1$ENDPOINT_class01 <- gsub('_NAS|HLP_', '', fglongi1$ENDPOINT)
fwrite(fglongi1, 'data/tx_cancer_1stCancer.tsv', sep = '\t')

#
fglongi1 <- fread('data/tx_cancer_1stCancer.tsv', data.table = F)


## --------------------------------------------------------
## survival analysis data preparation
## --------------------------------------------------------

# FG covariate data with FU info
cvr2 <- fread('/finngen/library-red/finngen_R12/analysis_covariates/R12_COV_V2.FID.txt.gz')
colnames(cvr2)
cvr2 <- cvr2[, c(2:3, 9, 14:154, 156:157, 158, 165, 167, 177, 180:189)] 

# read first events
kidney.first <- fread('results/kidney_cancer/kidney.first.tsv', data.table=F)
kidney.first %>% na.omit %>% dim # 1801

lung.first <- fread('results/kidney_cancer/lung.first.tsv', data.table=F)
lung.first %>% na.omit %>% dim # 233

other.first <- fread('results/kidney_cancer/othertx.first.tsv', data.table=F)
other.first %>% na.omit %>% dim # 2938

cancer.first <- fread('results/kidney_cancer/cancer.first.tsv', data.table=F)


# combine cancer and tx data 
# full join = take all kidney and cancer events
cancer.kidney.first <- full_join(kidney.first, cancer.first, by='FINNGENID', suffix=c('.kidney', '.cancer'))
cancer.kidney.first %>% filter(!is.na(EVENT_YEAR.kidney)) %>% dim # 1801
cancer.kidney.first %>% dim
cancer.kidney.first %>% head()
cancer.kidney.first <- dplyr::select(cancer.kidney.first, c(FINNGENID, EVENT_AGE.kidney, EVENT_YEAR.kidney,
                                                            EVENT_AGE.cancer, EVENT_YEAR.cancer))

# combine cancer and tx data 
# full join = take all lung and cancer events
cancer.lung.first <- full_join(lung.first, cancer.first, by='FINNGENID', suffix=c('.lung', '.cancer'))
cancer.lung.first %>% filter(!is.na(EVENT_YEAR.lung)) %>% dim # 233
cancer.lung.first %>% dim
cancer.lung.first %>% head()
cancer.lung.first <- dplyr::select(cancer.lung.first, c(FINNGENID, EVENT_AGE.lung, EVENT_YEAR.lung,
                                                            EVENT_AGE.cancer, EVENT_YEAR.cancer))

# combine cancer and tx data 
# full join = take all othertx and cancer events
cancer.other.first <- full_join(other.first, cancer.first, by='FINNGENID', suffix=c('.other', '.cancer'))
cancer.other.first %>% filter(!is.na(EVENT_YEAR.other)) %>% dim # 2938
cancer.other.first %>% dim
cancer.other.first %>% head()
cancer.other.first <- dplyr::select(cancer.other.first, c(FINNGENID, EVENT_AGE.other, EVENT_YEAR.other,
                                                        EVENT_AGE.cancer, EVENT_YEAR.cancer))

# combine cancer, tx, and FG data
# right join = take all covariates, incl. missing BMI, SMOKE
# kidney
cancer.kidney.first <- right_join(cancer.kidney.first, cvr2[, ], by=c('FINNGENID'='IID'))
cancer.kidney.first <- filter(cancer.kidney.first, 
                              !(FINNGENID %in% filter(cancer.kidney.first, EVENT_AGE.cancer<EVENT_AGE.kidney)$FINNGENID))
cancer.kidney.first %>% head
cancer.kidney.first %>% dim()
cancer.kidney.first %>% filter(!is.na(EVENT_YEAR.kidney)) %>% dim # 1546

# lung
cancer.lung.first <- right_join(cancer.lung.first, cvr2[, ], by=c('FINNGENID'='IID'))
cancer.lung.first <- filter(cancer.lung.first, 
                              !(FINNGENID %in% filter(cancer.lung.first, EVENT_AGE.cancer<EVENT_AGE.lung)$FINNGENID))
cancer.lung.first %>% head
cancer.lung.first %>% dim()
cancer.lung.first %>% filter(!is.na(EVENT_YEAR.lung)) %>% dim # 206

# other
cancer.other.first <- right_join(cancer.other.first, cvr2[, ], by=c('FINNGENID'='IID'))
cancer.other.first <- filter(cancer.other.first, 
                            !(FINNGENID %in% filter(cancer.other.first, EVENT_AGE.cancer<EVENT_AGE.other)$FINNGENID))
cancer.other.first %>% head
cancer.other.first %>% dim()
cancer.other.first %>% filter(!is.na(EVENT_YEAR.other)) %>% dim # 2260

# age at cancer, death or end of FU: time to event
# kidney
cancer.kidney.first$EVENT_OR_FU <- cancer.kidney.first %>% apply(., 1, function(x) {
  res <- x[c('AGE_AT_DEATH_OR_END_OF_FOLLOWUP', 'EVENT_AGE.cancer')]
  res[which.min(res)] %>% as.numeric
})
cancer.kidney.first %>% head
cancer.kidney.first %>% dim
# add cancer event 
cancer.kidney.first$EVENT <- ifelse(!is.na(cancer.kidney.first$EVENT_YEAR.cancer), 1, 0)
cancer.kidney.first %>% dim

# lung
cancer.lung.first$EVENT_OR_FU <- cancer.lung.first %>% apply(., 1, function(x) {
  res <- x[c('AGE_AT_DEATH_OR_END_OF_FOLLOWUP', 'EVENT_AGE.cancer')]
  res[which.min(res)] %>% as.numeric
})
cancer.lung.first %>% head
cancer.lung.first %>% dim # 206
# add cancer event 
cancer.lung.first$EVENT <- ifelse(!is.na(cancer.lung.first$EVENT_YEAR.cancer), 1, 0)
cancer.lung.first <- cancer.lung.first %>% filter(!is.na(EVENT_YEAR.lung))
# time from tx to fu or can
cancer.lung.first$TIME_TX_FU <- 
  cancer.lung.first$EVENT_OR_FU - cancer.lung.first$EVENT_AGE.lung 
# write
fwrite(cancer.lung.first, 'data/phenotypes/cancer.lung.first', sep='\t', na='NA')

# other
cancer.other.first$EVENT_OR_FU <- cancer.other.first %>% apply(., 1, function(x) {
  res <- x[c('AGE_AT_DEATH_OR_END_OF_FOLLOWUP', 'EVENT_AGE.cancer')]
  res[which.min(res)] %>% as.numeric
})
cancer.other.first %>% head
cancer.other.first %>% dim # 2260
# add cancer event 
cancer.other.first$EVENT <- ifelse(!is.na(cancer.other.first$EVENT_YEAR.cancer), 1, 0)
cancer.other.first <- cancer.other.first %>% filter(!is.na(EVENT_YEAR.other))
# time from tx to fu or can
cancer.other.first$TIME_TX_FU <- 
  cancer.other.first$EVENT_OR_FU - cancer.other.first$EVENT_AGE.other 
# write
fwrite(cancer.other.first, 'data/phenotypes/cancer.other.first', sep='\t', na='NA')


# add variable for having had tx
cancer.kidney.first$HASTX <- ifelse(!is.na(cancer.kidney.first$EVENT_YEAR.kidney), 1, 0)

# data subset for TX samples only
cancer.kidney.first.tx <- filter(cancer.kidney.first, HASTX==1)
cancer.kidney.first.tx %>% dim
# time from tx to fu or can
cancer.kidney.first.tx$TIME_TX_FU <- 
  cancer.kidney.first.tx$EVENT_OR_FU - cancer.kidney.first.tx$EVENT_AGE.kidney 

# data subset for non-TX samples only
cancer.kidney.first <- filter(cancer.kidney.first, HASTX==0)
cancer.kidney.first %>% dim
# endpoints for any transplantation
tmp <- fread('data/phenotypes/pheno.header') %>% colnames
tmp[which(grepl('transplant', tmp, ignore.case=T))]
which(grepl('transplant', tmp, ignore.case=T))[c(1, 5, 9, 13, 17, 21)]
system(
  paste0("zcat /finngen/library-red/finngen_R12/phenotype_1.0/data/finngen_R12_endpoint_1.0.txt.gz ",
         "| cut -f1,11149,12155,12187,12597,13559,17322 > data/phenotypes/R12_transpl.txt") 
)
transp <- fread("data/phenotypes/R12_transpl.txt", data.table = F)
transp$anytransp <- rowSums(transp[, -1], na.rm=T)
# keep samples that have no transplants
cancer.kidney.first <- filter(cancer.kidney.first, 
                              FINNGENID %in% (transp %>% filter(anytransp==0) %>% .$FINNGENID))


# write
fwrite(cancer.kidney.first,    'data/phenotypes/cancer.kidney.first', sep='\t', na='NA')
fwrite(cancer.kidney.first.tx, 'data/phenotypes/cancer.kidney.first.tx', sep='\t', na='NA')



## -----------------------------------------------
## FinnGen phenos to use in adjusting cox models
## -----------------------------------------------

# Adjust for cancer relevant endpoints:
# Z21_FAMILY_HISTORY_MALIG_NEOPLASM   12447
# Q17_CYSTIC_KIDNEY_DISEA   9896
# E4_DIABETES   17893
# Q17_OTHER_CONGEN_MALFO_KIDNEY   9904
# Q17_OTHER_CONGEN_MALFO_URINARY_SYSTEM   9908

tmp <- fread('data/phenotypes/pheno.header') %>% colnames
which(tmp %in% c('Z21_FAMILY_HISTORY_MALIG_NEOPLASM', 'Q17_CYSTIC_KIDNEY_DISEA', 'E4_DIABETES',
                 'Q17_OTHER_CONGEN_MALFO_KIDNEY', 'Q17_OTHER_CONGEN_MALFO_URINARY_SYSTEM'))
tmp[which(tmp %in% c('Z21_FAMILY_HISTORY_MALIG_NEOPLASM', 'Q17_CYSTIC_KIDNEY_DISEA', 'E4_DIABETES',
                     'Q17_OTHER_CONGEN_MALFO_KIDNEY', 'Q17_OTHER_CONGEN_MALFO_URINARY_SYSTEM'))
]
system(
  paste0("zcat /finngen/library-red/finngen_R12/phenotype_1.0/data/finngen_R12_endpoint_1.0.txt.gz ",
         "| cut -f1,9896,9904,9908,12447,17893 > data/phenotypes/TxCan_cox_adj.tsv") 
)


## -----------------------------------------------
## Longitudinal data
## -----------------------------------------------

# all tx patients ID's
write(cancer.kidney.first.tx %>% .$FINNGENID, 'data/tx_all_ids', ncolumns=1)

# extract tx longitudinal data
system("zgrep -F -f data/tx_all_ids /finngen/library-red/finngen_R12/phenotype_1.0/data/finngen_R12_endpoint_longitudinal_1.0.txt.gz > data/phenotypes/TX_LONG_phe.tsv")
tx.longi <- fread('data/phenotypes/TX_LONG_phe.tsv')
colnames(tx.longi) <- c('FINNGENID',	'EVENT_TYPE',	'EVENT_AGE',	'EVENT_YEAR',	'ICDVER',	'ENDPOINT')
fwrite(tx.longi, 'data/phenotypes/TX_LONG_phe.tsv')

# extract non-kidney tx longitudinal data
system("zgrep -F -f data/tx_ids2 /finngen/library-red/finngen_R12/phenotype_1.0/data/finngen_R12_endpoint_longitudinal_1.0.txt.gz > data/phenotypes/TX_LONG_phe2.tsv")
tx.longi <- fread('data/phenotypes/TX_LONG_phe2.tsv')
colnames(tx.longi) <- c('FINNGENID',	'EVENT_TYPE',	'EVENT_AGE',	'EVENT_YEAR',	'ICDVER',	'ENDPOINT')
fwrite(tx.longi, 'data/phenotypes/TX_LONG_phe2.tsv')


# extract immunosuppression med 
# L04AA Selective immunisuppressants
# L04AB TNFa inhibitors
# L04AB01 Etanercept
# L04AB02 Infliximab
# L04AB04 Adalimumab
# L04AC Interleukin inhibitors
# L04AD Calcineurin inhibitors
# L04AX01 Azathioprine
tx.longi.immsupp <- map_dfr(tx.longi$FINNGENID %>% unique, function(x) {
  tx.longi %>% filter(FINNGENID==x, 
                      grepl('IMMUNOS|L04A|AZATHIOPRINE|ILD_CORTISONE', ENDPOINT)) %>% 
    .$ENDPOINT %>% table
}) %>% data.frame
tx.longi.immsupp <- data.frame(FINNGENID=tx.longi$FINNGENID %>% unique,
                               tx.longi.immsupp)
tx.longi.immsupp[is.na(tx.longi.immsupp)] <- 0
# how many entries per ID?
tx.longi %>% group_by(FINNGENID) %>% summarise(N=n())
# correlation
tx.longi.immsupp[, -1] %>% cor() %>% corrplot::corrplot()
# remove redundants
tx.longi.immsupp <- tx.longi.immsupp %>% dplyr::select(-c(ILD_AZATHIOPRINE, RX_L04AX))
# write
fwrite(tx.longi.immsupp, 'data/phenotypes/TX_LONG_phe_immsupp.tsv')

# extract endpoints before the first tx
tx.longi.primdiag <- tx.longi$FINNGENID %>% unique %>% map_dfr(function(x) {
  res <- filter(tx.longi, FINNGENID==x)
  ind <- which(res$ENDPOINT=="Z21_KIDNEY_TRANSPLANT_STATUS")[1]
  res <- res[1:(ind), ]
  res$EVENT_YEAR <- as.Date(res$EVENT_YEAR)
  res$TimeDiffD_fromTX <- (res$EVENT_YEAR[ind] - res$EVENT_YEAR) %>% as.numeric
  res$TimeDiffY_fromTX <- res$TimeDiffD_fromTX/365 %>% as.numeric
  res %>% return()
})
fwrite(tx.longi.primdiag, 'data/phenotypes/TX_LONG_phe_primary.tsv', sep='\t')


## --------------------------------------------------------
## TX patients cancer diagnoses
## --------------------------------------------------------

# tx patients with cancer event
write(cancer.kidney.first.tx %>% filter(EVENT==1) %>% .$FINNGENID, 'data/tx_event_ids', ncolumns=1)

# extract phenotype data from tx patients with cancer event
system(paste0("zgrep -F -f data/tx_event_ids ",
              "/finngen/library-red/finngen_R12/phenotype_1.0/data/finngen_R12_endpoint_1.0.txt.gz > ",
              "data/phenotypes/TX_EVENT_phe.tsv"))
phe <- fread('data/phenotypes/TX_EVENT_phe.tsv', data.table=F)
colnames(phe) <- fread('data/phenotypes/pheno.header') %>% colnames
phe <- phe[, !grepl('_FU_AGE$|_APPROX_EVENT_DAY$|_NEVT$|_EXALLC$|BENIGN|INCLAVO', colnames(phe))]
phe <- phe[, grepl('FINNGENID|CD2_|C3_', colnames(phe))] 
fwrite(phe, 'data/phenotypes/TX_EVENT_cancers.tsv', sep='\t')


# cancers in tx patients
phe <- fread('data/phenotypes/TX_EVENT_cancers.tsv', data.table=F)

# melanoma vs other
# phe[, grepl('melano', phe %>% colnames(), ignore.case=T)] %>% rowSums()
# phe[, !grepl('melano', phe %>% colnames(), ignore.case=T)][, -1] %>% rowSums(na.rm=T)
hasMelanoma <- ifelse(phe$C3_MELANOMA==1, 'Melanoma', 'OtherCancer')
hasMelanoma <- data.frame(FINNGENID=phe$FINNGENID, Melanoma=hasMelanoma)
fwrite(hasMelanoma, 'data/phenotypes/hasMelanoma.tsv', sep='\t')



## --------------------------------------------------------
## TX patients dialysis prior to KT
## --------------------------------------------------------

# computes both the time in years from dialysis start age to end age
# and the number of unique dialysis treatments per patient

# read longitudinal pheno data for KT patients
tx.longi <- fread('data/phenotypes/TX_LONG_phe.tsv', data.table = F)

# extract outpatient register start and end ages of dialysis
dialysis.times <- map_dfr(tx.longi$FINNGENID %>% unique, function(x) {
  print(x)
  
  res <- tx.longi %>% filter(FINNGENID == x, 
                             grepl('DIALY', ENDPOINT), 
                             EVENT_TYPE %in% c('ERIK_AVO', 'ERIK_OPER')) %>% 
    arrange(EVENT_AGE)
  if(nrow(res)==0) diatime <- 0 else diatime <- res$EVENT_AGE[nrow(res)] - res$EVENT_AGE[1]
  
  res2 <- tx.longi %>% filter(FINNGENID == x, 
                              grepl('DIALY', ENDPOINT)) %>% 
    arrange(EVENT_AGE)
  
  data.frame(FINNGENID = x,
             Dialysis_time = diatime,
             Unique_dialysis_entries = res2$EVENT_AGE %>% unique %>% length) %>% return()
})
fwrite(dialysis.times, 'results/kidney_cancer/KT_dialysis_times.tsv', sep = '\t', quote = F)

dialysis.times <- fread('results/kidney_cancer/KT_dialysis_times.tsv', data.table = F)

dialysis.times$Dialysis_time %>% hist
dialysis.times$Unique_dialysis_entries %>% hist

(dialysis.times$Dialysis_time+1) %>% log %>% hist
(dialysis.times$Unique_dialysis_entries+1) %>% log %>% hist

cor(dialysis.times$Dialysis_time, dialysis.times$Unique_dialysis_entries)
plot(dialysis.times$Dialysis_time, dialysis.times$Unique_dialysis_entries)

cor((dialysis.times$Dialysis_time+1) %>% log, (dialysis.times$Unique_dialysis_entries+1) %>% log)
plot((dialysis.times$Dialysis_time+1) %>% log, (dialysis.times$Unique_dialysis_entries+1) %>% log)



## --------------------------------------------------------
## cystic kidney diagnoses
## --------------------------------------------------------

dialysis.times <- fread('results/kidney_cancer/KT_dialysis_times.tsv', data.table = F)
fglongi <- fread('/finngen/library-red/finngen_R12/phenotype_1.0/data/finngen_R12_endpoint_longitudinal_1.0.txt.gz')
fglongi <- filter(fglongi, FINNGENID %in% dialysis.times$FINNGENID) 
fglongi %>% dim

tx.ckd.data <- map_dfr(fglongi$FINNGENID %>% unique, function(x) {
  print(x)
  tmp <- filter(fglongi, FINNGENID == x)
  tx.age <- filter(tmp, ENDPOINT == 'Z21_KIDNEY_TRANSPLANT_STATUS') %>% .$EVENT_AGE %>% .[1]
  tmp <- filter(tmp, ENDPOINT == 'Q17_CYSTIC_KIDNEY_DISEA', EVENT_AGE < tx.age)
  if(nrow(tmp) > 0) cystic.period <- tmp$EVENT_AGE[nrow(tmp)] - tmp$EVENT_AGE[1] else cystic.period <- 0
  data.frame(FINNGENID = x,
             Years_ckd = cystic.period,
             Number_uniq_ckd_entries = tmp$EVENT_AGE %>% unique %>% length) %>% 
    return()
})
fwrite(tx.ckd.data, 'results/kidney_cancer/KT_CKD_times.tsv', sep = '\t', quote = F)

tx.ckd.data$Years_ckd %>% hist
tx.ckd.data$Number_uniq_ckd_entries %>% hist

(tx.ckd.data$Years_ckd+1) %>% log %>% hist
(tx.ckd.data$Number_uniq_ckd_entries+1) %>% log %>% hist

cor((tx.ckd.data$Years_ckd), (tx.ckd.data$Number_uniq_ckd_entries))
cor((tx.ckd.data$Years_ckd+1) %>% log, (tx.ckd.data$Number_uniq_ckd_entries+1) %>% log)

plot((tx.ckd.data$Years_ckd), (tx.ckd.data$Number_uniq_ckd_entries))
plot((tx.ckd.data$Years_ckd+1) %>% log, (tx.ckd.data$Number_uniq_ckd_entries+1) %>% log)

plot((dialysis.times$Dialysis_time), (tx.ckd.data$Number_uniq_ckd_entries))
plot((dialysis.times$Dialysis_time), (tx.ckd.data$Years_ckd))


inner_join(dialysis.times, tx.ckd.data) %>% 
  mutate(Dialysis_time = Dialysis_time + 1,
         Unique_dialysis_entries = Unique_dialysis_entries +1,
         Years_ckd = Years_ckd + 1,
         Number_uniq_ckd_entries = Number_uniq_ckd_entries + 1) %>% 
  mutate(Dialysis_time_log = Dialysis_time %>% log,
         Unique_dialysis_entries_log = Unique_dialysis_entries %>% log,
         Years_ckd_log = Years_ckd %>% log,
         Number_uniq_ckd_entries_log = Number_uniq_ckd_entries %>% log) %>% 
  dplyr::select(-FINNGENID) %>% cor %>% corrplot

