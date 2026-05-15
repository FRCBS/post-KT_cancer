
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


## ========================================================
##
## TX cancer survival analysis 
## PRS and malignant cancer 
##
## ========================================================


## functions
fread2 <- function(x) fread(x, data.table = F)

signif2 <- function(x, s=3) {
  data.frame(variable = x[, 1], 
             x[, -1] %>% signif(s))
}


## --------------------------------------------------------
## data read-in
## --------------------------------------------------------

# longitudinal event data
cancer.kidney.first    <- fread2('data/phenotypes/cancer.kidney.first')
cancer.kidney.first.tx <- fread2('data/phenotypes/cancer.kidney.first.tx')

# pancancer SNP risk scores
risk.scores <- rbind(
  fread2('data/FG_pancancer_scores_01.tsv'),
  fread2('data/FG_pancancer_scores_02.tsv')
)

# adjustment phenotypes
adj.phenos <- fread2('data/phenotypes/TxCan_cox_adj.tsv')

# tx longitudinal endpoints
tx.longi <- fread2('data/phenotypes/TX_LONG_phe.tsv')
tx.longi.immsupp <- fread2('data/phenotypes/TX_LONG_phe_immsupp.tsv')

# get cancer incidence years from cancer registry
cancer <- fread("/finngen/library-red/finngen_R12/cancer_detailed_1.0/data/finngen_R12_cancer_detailed_1.0.txt", 
                data.table = F) 
cancer %>% head


## --------------------------------------------------------
## PRS and risk of malignant cancer 
## within non-tx FG cancer patients
## --------------------------------------------------------

cancer.kidney.first.model <- cancer.kidney.first %>% filter(BL_AGE<100) %>% 
  dplyr::select(c(FINNGENID, BL_AGE, BL_YEAR, 
                  PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                  EVENT_OR_FU, EVENT, SEX, BMI, SMOKE2))
cancer.kidney.first.model %>% dim
cancer.kidney.first.model$EVENT %>% table
cancer.kidney.first.model %>% str

# join with other covars
cancer.kidney.first.model <- inner_join(cancer.kidney.first.model, risk.scores[, c('IID', 'sum')], 
                                        by=c('FINNGENID'='IID'))
# cancer data
cancer.kidney.first.model <- left_join(cancer.kidney.first.model,
                                       dplyr::select(cancer, c(FINNGENID, beh)))

# family history data
cancer.kidney.first.model <- left_join(cancer.kidney.first.model, adj.phenos[, c(1, 5:6)])
cancer.kidney.first.model <- cancer.kidney.first.model %>% dplyr::select(-FINNGENID)
# NA = 0 in fam history data
cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM[
  is.na(cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM)] <- 0
cancer.kidney.first.model$EVENT %>% table

# numeric
cancer.kidney.first.model <- model.matrix( ~. , data = cancer.kidney.first.model) %>% 
  .[, -1] %>% data.frame()
cancer.kidney.first.model %>% dim
cancer.kidney.first.model %>% head
cancer.kidney.first.model$EVENT %>% table # only events

# Cancer classification
# beh					ICD-O-3 behaviour
# 0 benign
# 1 semimalignant
# 2 carcinoma in situ
# 3 malignant

cancer.kidney.first.model$behBvsM <- NA
cancer.kidney.first.model$behBvsM[cancer.kidney.first.model$beh<2] <- 0 
cancer.kidney.first.model$behBvsM[cancer.kidney.first.model$beh==3] <- 1 
cancer.kidney.first.model$behBvsM %>% table
cancer.kidney.first.model$beh %>% table

# cox model
# benign/insitu vs. malignant
res.surv <- coxph(Surv(EVENT_OR_FU, behBvsM)  ~ ., 
                  data = cancer.kidney.first.model %>% 
                    dplyr::select(-c(EVENT, beh)))
res.surv %>% summary %>% coef

tbl_regression(res.surv, exp=T)$table_body %>% data.frame() %>% 
  dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                  statistic, conf.low, conf.high, p.value)) %>% signif2(4)


