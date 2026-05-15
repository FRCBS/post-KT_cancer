
library(tidyverse)
library(data.table)
library(corrplot)
library(survival)
library(survminer)
library(knitr)
library(kableExtra)
library(ggsci)
library(ggpubr)
library(patchwork)
library(gtsummary)
library(readxl)
library(broom)


## ========================================================
##
## FinnGen survival analysis for cancer SNP score
##
## ========================================================

fread2 <- function(x) fread(x, data.table = F)

signif2 <- function(x) {
  data.frame(variable = x[, 1], 
             x[, -1] %>% signif(3))
}


## --------------------------------------------------------
## data read-in
## --------------------------------------------------------

# longitudinal event data
cancer.kidney.first    <- fread('data/phenotypes/cancer.kidney.first')
cancer.kidney.first.tx <- fread('data/phenotypes/cancer.kidney.first.tx')

# live tx IDs
liverid <- fread2('data/phenotypes/R13_longitudinal_LiverTxIDs.txt')$V1

# heart and liver IDs
# lhid <- c(fread2('data/phenotypes/R13_longitudinal_LiverTxIDs.txt')$V1,
#           fread2('data/phenotypes/R13_longitudinal_HeartTxIDs.txt')$V1)
lhid <- c(fread2('data/phenotypes/R13_longitudinal_HeartTxIDs.txt')$V1)

# longitudinal event data
cancer.lung.first <- fread2('data/phenotypes/cancer.lung.first')
cancer.other.first <- fread2('data/phenotypes/cancer.other.first') %>% 
  filter(FINNGENID %in% lhid)
cancer.liver.first <- fread2('data/phenotypes/cancer.other.first') %>% 
  filter(FINNGENID %in% liverid)

cancer.comb.first <- rbind( 
  cancer.lung.first %>% rename(EVENT_AGE.comb = EVENT_AGE.lung) %>%
    select(-EVENT_YEAR.lung) %>% mutate(TX = 'lung'),
  cancer.liver.first %>% rename(EVENT_AGE.comb = EVENT_AGE.other) %>%
    select(-EVENT_YEAR.other) %>% mutate(TX = 'liver')
)


# pancancer SNP risk scores
risk.scores <- rbind(
  fread('data/FG_pancancer_scores_01.tsv'),
  fread('data/FG_pancancer_scores_02.tsv')
)

# adjustment pehotypes
adj.phenos <- fread('data/phenotypes/TxCan_cox_adj.tsv')

# tx longitudinal endpoints
tx.longi <- fread('data/phenotypes/TX_LONG_phe.tsv')
tx.longi.immsupp <- fread('data/phenotypes/TX_LONG_phe_immsupp.tsv')

# tx melanoma vs other cancers
hasMelanoma <- right_join(fread('data/phenotypes/hasMelanoma.tsv'), tx.longi.immsupp[, 1:2])[, -3]
hasMelanoma$Melanoma[is.na(hasMelanoma$Melanoma)] <- 'NoCancer' 

# FinnGen pheno names
fgp <- read_xlsx('data/phenotypes/FINNGEN_ENDPOINTS_DF12_Final_2023-05-17_public.xlsx')[, 1:4]



## --------------------------------------------------------
## FinnGen 
## validation of cancer scores
## --------------------------------------------------------

## cancer score cox survival models
## All adjustments, including 
## Z21_FAMILY_HISTORY_MALIG_NEOPLASM & DIAB

dat <- rbind(
  cancer.kidney.first %>% 
    dplyr::select(c(FINNGENID, BL_YEAR, BL_AGE, cohort,
                    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                    EVENT_OR_FU, EVENT,
                    SEX)) %>% na.omit %>% mutate(TX = 0),
  cancer.kidney.first.tx %>% 
    dplyr::select(c(FINNGENID, BL_YEAR, BL_AGE, cohort,
                    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                    EVENT_OR_FU, EVENT,
                    SEX)) %>% na.omit %>% mutate(TX = 1),
  cancer.comb.first %>% 
    dplyr::select(c(FINNGENID, BL_YEAR, BL_AGE, cohort,
                    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                    EVENT_OR_FU, EVENT,
                    SEX)) %>% na.omit %>% mutate(TX = 1)
)


# prepare data for cox
cancer.kidney.first.model <- left_join(dat, adj.phenos)
cancer.kidney.first.model <- left_join(cancer.kidney.first.model, risk.scores,
                                       by=c('FINNGENID'='IID')) %>% na.omit

# # adjustment phenotypes
# adj.phenos <- fread2('data/phenotypes/TxCan_cox_adj.tsv')
# cancer.kidney.first.model <- inner_join(cancer.kidney.first.model, adj.phenos)

# numeric matrix
cancer.kidney.first.model <- model.matrix( ~. ,
                                           data = cancer.kidney.first.model %>% 
                                             dplyr::select(-FINNGENID)) %>% .[, -1] %>% data.frame()


# survival
res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ BL_AGE + BL_YEAR + sum + SEXmale +  
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + E4_DIABETES + TX*sum +
                    Z21_FAMILY_HISTORY_MALIG_NEOPLASM,
                  data = cancer.kidney.first.model, iter.max = 1000)
res.surv %>% summary %>% coef %>% na.omit
fwrite(tbl_regression(res.surv, exp=T)$table_body %>% data.frame() %>%  
         dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                         statistic, conf.low, conf.high, p.value)) %>% signif2,
       'results/kidney_cancer/fulldata_cox_noBMI_sumXtx.tsv', sep = '\t'
)


res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ BL_AGE + BL_YEAR + sum + SEXmale +  
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + E4_DIABETES + TX*sum +
                    Z21_FAMILY_HISTORY_MALIG_NEOPLASM,
                  data = cancer.kidney.first.model %>% filter(BL_AGE < 40), iter.max = 1000)
res.surv %>% summary %>% coef %>% na.omit
fwrite(tbl_regression(res.surv, exp=T)$table_body %>% data.frame() %>%  
         dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                         statistic, conf.low, conf.high, p.value)) %>% signif2,
       'results/kidney_cancer/fulldata_cox_under40blage_noBMI_sumXtx.tsv', sep = '\t'
)





