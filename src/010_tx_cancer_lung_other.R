
library(tidyverse)
library(dplyr)
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
## TX cancer survival analysis F_2023_043
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
# cancer.kidney.first    <- fread2('data/phenotypes/cancer.kidney.first')
cancer.kidney.first.tx <- fread2('data/phenotypes/cancer.kidney.first.tx')

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
  # cancer.kidney.first.tx %>% rename(EVENT_AGE.comb = EVENT_AGE.kidney) %>%
  #   select(-EVENT_YEAR.kidney, -HASTX) %>% mutate(TX = 'kidney'),
  cancer.lung.first %>% rename(EVENT_AGE.comb = EVENT_AGE.lung) %>%
    select(-EVENT_YEAR.lung) %>% mutate(TX = 'lung'),
  # cancer.other.first %>% rename(EVENT_AGE.comb = EVENT_AGE.other) %>%
  #   select(-EVENT_YEAR.other) %>% mutate(TX = 'heart'),
  cancer.liver.first %>% rename(EVENT_AGE.comb = EVENT_AGE.other) %>%
    select(-EVENT_YEAR.other) %>% mutate(TX = 'liver')
  )
cancer.comb.first$FINNGENID %>% write('data/tx_ids2', ncolumns = 1)
cancer.comb.first %>% dim # 2022

# pancancer SNP risk scores
risk.scores <- rbind(
  fread2('data/FG_pancancer_scores_01.tsv'),
  fread2('data/FG_pancancer_scores_02.tsv')
)

# adjustment phenotypes
adj.phenos <- fread2('data/phenotypes/TxCan_cox_adj.tsv')

# tx longitudinal endpoints
# tx.longi <- fread2('data/phenotypes/TX_LONG_phe.tsv')
tx.longi.immsupp <- fread2('data/phenotypes/TX_LONG_phe_immsupp.tsv')

# tx melanoma vs other cancers
hasMelanoma <- right_join(fread2('data/phenotypes/hasMelanoma.tsv'), tx.longi.immsupp[, 1:2])[, -3]
hasMelanoma$Melanoma[is.na(hasMelanoma$Melanoma)] <- 'NoCancer' 

# # cancers in tx patients
# # all cancers per individual
# phe <- fread2('data/phenotypes/TX_EVENT_cancers.tsv')
# 
# # finngen pheno descriptions
# fgp <- read_xlsx('data/phenotypes/FINNGEN_ENDPOINTS_DF12_Final_2023-05-17_public.xlsx')[, 1:4]
# 
# # cancer diag grouping
# cangroup <- fread2('data/phenotypes/cancer_grouping.tsv')
# 
# # dialysis times
# dialysis.times <- fread('results/kidney_cancer/KT_dialysis_times.tsv', data.table = F)
# dialysis.times$Unique_dialysis_entries_log <- (dialysis.times$Unique_dialysis_entries+1) %>% log


## --------------------------------------------------------
## BMI and smoke info available
## --------------------------------------------------------

# cancer.kidney.first.tx %>%  dplyr::select(c(SMOKE2, BMI)) %>% na.omit %>% dim




## --------------------------------------------------------
## survival; genetic risks in lung and liver tx
## --------------------------------------------------------

## risk scores in tx age <40
## no adj for BMI, SMO
# family history data

# data
cancer.comb.first.model <- cancer.comb.first %>% filter(EVENT_AGE.comb>1) %>% 
  dplyr::select(c(FINNGENID, BL_AGE, BL_YEAR, EVENT_AGE.comb,
                  PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                  EVENT_OR_FU, EVENT, # TIME_TX_FU
                  SEX, TX))#, SMOKE2, BMI))


# # number of ILD_IMMUNOSUPPRESSIVE events per individual
# immuno.comb <- filter(immuno, FINNGENID %in% cancer.comb.first.model$FINNGENID)
# immuno.comb$FINNGENID %>% unique %>% length
# immuno.comb %>% head
# cancer.comb.first.model <- immuno.comb$FINNGENID %>% unique %>% map_dfr(function(x) {
#   tmp01 <- immuno.comb %>% filter(FINNGENID == x)
#   tmp02 <- cancer.comb.first.model %>% filter(FINNGENID == x)
#   tmp01 <- tmp01 %>% filter(EVENT_AGE > tmp02$EVENT_AGE.comb[1])
#   tmp02$NUM_IMMUNOSUPP <- nrow(tmp01)
#   tmp02
# })
# cancer.comb.first.model$NUM_IMMUNOSUPP %>% hist


# join genetic risk score data
cancer.comb.first.model <- inner_join(cancer.comb.first.model, risk.scores[, c('IID', 'sum')], 
                                      by=c('FINNGENID'='IID'))
# family history data
cancer.comb.first.model <- left_join(cancer.comb.first.model, adj.phenos[, c(1, 5:6)])
# NA = 0 in fam history data
cancer.comb.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM[
  is.na(cancer.comb.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM)] <- 0
cancer.comb.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM %>% sum

# check
cancer.comb.first.model %>% dim
cancer.comb.first.model %>% colnames

# cancer types in DIAB/non-DIAB
tx.longi <- fread('data/phenotypes/TX_LONG_phe2.tsv')

can.diab <- tx.longi %>% filter(FINNGENID %in% 
                      (cancer.comb.first.model %>% filter(E4_DIABETES == 1) %>% .$FINNGENID)) %>% 
  filter(grepl('C3|CD2', ENDPOINT)) %>% .$ENDPOINT %>% table %>% sort

can.nondiab <- tx.longi %>% filter(FINNGENID %in% 
                      (cancer.comb.first.model %>% filter(E4_DIABETES == 0) %>% .$FINNGENID)) %>% 
  filter(grepl('C3|CD2', ENDPOINT)) %>% .$ENDPOINT %>% table %>% sort

inner_join(
  (sort(can.diab)/sum(can.diab)) %>% tail(100) %>% data.frame %>% rename('Can' = '.'),
  (sort(can.nondiab)/sum(can.nondiab)) %>% tail(100) %>% data.frame %>% rename('Can' = '.'),
  by = 'Can', suffix = c('_Diab', '_NonDiab')
) %>% ggplot(aes(Freq_Diab %>% log10, Freq_NonDiab %>% log10,  label = Can)) +
  #geom_text(size = 1.8) +
  geom_point() +
  geom_abline(slope = 1) +
  geom_label(size = 2)


# remove ID from orig
cancer.comb.first.model <- cancer.comb.first.model %>% dplyr::select(-FINNGENID)

# format variable matrix to numeric
cancer.comb.first.model <- model.matrix( ~. , data = cancer.comb.first.model) %>% 
  .[, -1] %>% data.frame()
cancer.comb.first.model %>% dim
cancer.comb.first.model %>% head
cancer.comb.first.model %>% colnames

# cox model
res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ EVENT_AGE.comb + sum + SEXmale +  TXlung + #sum*TXliver + sum*TXlung + 
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + E4_DIABETES + E4_DIABETES*sum +
                    Z21_FAMILY_HISTORY_MALIG_NEOPLASM,
                  data = cancer.comb.first.model %>%
                    filter(EVENT_AGE.comb>1, EVENT_AGE.comb<40), iter.max = 1000)
res.surv %>% summary %>% coef %>% na.omit
fwrite(tbl_regression(res.surv, exp=T)$table_body %>% data.frame() %>%  
         dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                         statistic, conf.low, conf.high, p.value)) %>% signif2,
       'results/kidney_cancer/lung_liver_cox_under40tx_noBMI_sumXdiab.tsv', sep = '\t'
)


cancer.comb.first.model %>% filter(EVENT_AGE.comb < 400) %>% .$E4_DIABETES %>% table # 0:313 1:152
cancer.comb.first.model %>% filter(EVENT_AGE.comb < 400) %>% .$EVENT %>% table # 0:366 1:99
cancer.comb.first.model %>% filter(EVENT_AGE.comb < 400) %>% .$TXlung %>% table # 0:262 1:203

cancer.comb.first.model %>% filter(EVENT_AGE.comb < 40) %>% .$E4_DIABETES %>% table # 0:85 1:23
cancer.comb.first.model %>% filter(EVENT_AGE.comb < 40) %>% .$EVENT %>% table # 0:95 1:13
cancer.comb.first.model %>% filter(EVENT_AGE.comb < 40) %>% .$TXlung %>% table # 0:70 1:38

jpeg('results/kidney_cancer/Liver_lung_cox.jpg', width = 6, height = 4.5, res = 1000, units = 'in')
res.surv %>% tidy(exponentiate = T) %>% na.omit %>% filter(!grepl('PC', term)) %>%
  arrange(estimate) %>% 
  mutate(term = ifelse(term == 'EVENT_AGE.comb', 'Event_age', term), 
         term = factor(term, levels = term),
         Significance = ifelse(p.value < 0.05, 'p < 0.05', 'n.s.')) %>% 
  ggplot(aes(estimate,  term, xmin = estimate-(std.error*1.96), 
             xmax = estimate+(std.error*1.96), shape = Significance)) +
  geom_pointrange() +
  xlab('HR') +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  #scale_x_continuous(transform = 'log10') +
  scale_shape_manual(values = c(1, 19)) +
  ggpubr::theme_classic2() +
  theme(axis.text = element_text(color = 'black'))
dev.off()


# res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ EVENT_AGE.comb + sum + SEXmale +  TXlung + #sum*TXliver + sum*TXlung + 
#                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
#                   data = cancer.comb.first.model %>%
#                     filter(EVENT_AGE.comb>1, EVENT_AGE.comb<400, E4_DIABETES==1), iter.max = 10000)
# res.surv %>% summary %>% coef %>% na.omit


# primary diagnoses
tx.longi <- fread('data/phenotypes/TX_LONG_phe2.tsv')
cancer.comb.first %>% head

# diagnoses
diag.llt <- tx.longi$FINNGENID %>% unique %>% map_dfr(function(x) {
  tmp <- tx.longi %>% filter(FINNGENID == x, EVENT_TYPE=='HILMO') %>% .$ENDPOINT
  data.frame(FINNGENID = x,
             Hepatitis = grepl('K11_CHRONHEP', tmp, ignore.case = T) %>% any,
             Cirrhosis = grepl('cirrhosis', tmp, ignore.case = T) %>% any,
             Alcohol = grepl('alcohol', tmp, ignore.case = T) %>% any,
             COPD = grepl('copd', tmp, ignore.case = T) %>% any,
             ILD = grepl('ILD', tmp, ignore.case = F) %>% any,
             Bronchiectasis = grepl('J10_BRONCHIECTASIS', tmp, ignore.case = F) %>% any,
             Diabetes = grepl('diab|t1d|t2d', tmp, ignore.case = T) %>% any)
})
cancer.comb.first <- left_join(cancer.comb.first, diag.llt)

# lung tx patients
all.llt <- data.frame(male_sex_n = sum(filter(cancer.comb.first, TX =='lung')$SEX=='male'),
                      male_sex_perc = 100*(sum(filter(cancer.comb.first, TX =='lung')$SEX=='male') / 
                                             length(filter(cancer.comb.first, TX =='lung')$SEX)),
                      tx_age_median = filter(cancer.comb.first, TX =='lung')$EVENT_AGE.comb %>% median,
                      tx_age_iqr = filter(cancer.comb.first, TX =='lung')$EVENT_AGE.comb %>% IQR,
                      years_from_tx_median = filter(cancer.comb.first, TX =='lung')$TIME_TX_FU %>% median,
                      years_from_tx_iqr = filter(cancer.comb.first, TX =='lung')$TIME_TX_FU %>% IQR,
                      Alcohol_n = sum(filter(cancer.comb.first, TX =='lung')$Alcohol==T),
                      Alcohol_perc = 100*(sum(filter(cancer.comb.first, TX =='lung')$Alcohol==T) / 
                                            length(filter(cancer.comb.first, TX =='lung')$Alcohol)),
                      COPD_n = sum(filter(cancer.comb.first, TX =='lung')$COPD==T),
                      COPD_perc = 100*(sum(filter(cancer.comb.first, TX =='lung')$COPD==T) / 
                                         length(filter(cancer.comb.first, TX =='lung')$COPD)),
                      ILD_n = sum(filter(cancer.comb.first, TX =='lung')$ILD==T),
                      ILD_perc = 100*(sum(filter(cancer.comb.first, TX =='lung')$ILD==T) / 
                                        length(filter(cancer.comb.first, TX =='lung')$ILD)),
                      Cirrhosis_n = sum(filter(cancer.comb.first, TX =='lung')$Cirrhosis==T),
                      Cirrhosis_perc = 100*(sum(filter(cancer.comb.first, TX =='lung')$Cirrhosis==T) / 
                                        length(filter(cancer.comb.first, TX =='lung')$Cirrhosis)),
                      Hepatitis_n = sum(filter(cancer.comb.first, TX =='lung')$Hepatitis==T),
                      Hepatitis_perc = 100*(sum(filter(cancer.comb.first, TX =='lung')$Hepatitis==T) / 
                                              length(filter(cancer.comb.first, TX =='lung')$Hepatitis)),
                      Bronchiectasis_n = sum(filter(cancer.comb.first, TX =='lung')$Bronchiectasis==T),
                      Bronchiectasis_perc = 100*(sum(filter(cancer.comb.first, TX =='lung')$Bronchiectasis==T) / 
                                              length(filter(cancer.comb.first, TX =='lung')$Bronchiectasis)),
                      Diabetes_n = sum(filter(cancer.comb.first, TX =='lung')$Diabetes==T),
                      Diabetes_perc = 100*(sum(filter(cancer.comb.first, TX =='lung')$Diabetes==T) / 
                                                   length(filter(cancer.comb.first, TX =='lung')$Diabetes))
)
canc.llt <- data.frame(male_sex_n = sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$SEX=='male'),
                       male_sex_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$SEX=='male') / 
                                              length(filter(cancer.comb.first, EVENT==1, TX =='lung')$SEX)),
                       tx_age_median = filter(cancer.comb.first, EVENT==1, TX =='lung')$EVENT_AGE.comb %>% median,
                       tx_age_iqr = filter(cancer.comb.first, EVENT==1, TX =='lung')$EVENT_AGE.comb %>% IQR,
                       years_from_tx_median = filter(cancer.comb.first, EVENT==1, TX =='lung')$TIME_TX_FU %>% median,
                       years_from_tx_iqr = filter(cancer.comb.first, EVENT==1, TX =='lung')$TIME_TX_FU %>% IQR,
                       Alcohol_n = sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$Alcohol==T),
                       Alcohol_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$Alcohol==T) / 
                                             length(filter(cancer.comb.first, EVENT==1, TX =='lung')$Alcohol)),
                       COPD_n = sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$COPD==T),
                       COPD_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$COPD==T) / 
                                          length(filter(cancer.comb.first, EVENT==1, TX =='lung')$COPD)),
                       ILD_n = sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$ILD==T),
                       ILD_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$ILD==T) / 
                                         length(filter(cancer.comb.first, EVENT==1, TX =='lung')$ILD)),
                       Cirrhosis_n = sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$Cirrhosis==T),
                       Cirrhosis_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$Cirrhosis==T) / 
                                               length(filter(cancer.comb.first, EVENT==1, TX =='lung')$Cirrhosis)),
                       Hepatitis_n = sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$Hepatitis==T),
                       Hepatitis_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$Hepatitis==T) / 
                                               length(filter(cancer.comb.first, EVENT==1, TX =='lung')$Hepatitis)),
                       Bronchiectasis_n = sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$Bronchiectasis==T),
                       Bronchiectasis_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$Bronchiectasis==T) / 
                                                    length(filter(cancer.comb.first, EVENT==1, TX =='lung')$Bronchiectasis)),
                       Diabetes_n = sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$Diabetes==T),
                       Diabetes_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='lung')$Diabetes==T) / 
                                                    length(filter(cancer.comb.first, EVENT==1, TX =='lung')$Diabetes))
)
no.llt <- data.frame(male_sex_n = sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$SEX=='male'),
                     male_sex_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$SEX=='male') / 
                                            length(filter(cancer.comb.first, EVENT==0, TX =='lung')$SEX)),
                     tx_age_median = filter(cancer.comb.first, EVENT==0, TX =='lung')$EVENT_AGE.comb %>% median,
                     tx_age_iqr = filter(cancer.comb.first, EVENT==0, TX =='lung')$EVENT_AGE.comb %>% IQR,
                     years_from_tx_median = filter(cancer.comb.first, EVENT==0, TX =='lung')$TIME_TX_FU %>% median,
                     years_from_tx_iqr = filter(cancer.comb.first, EVENT==0, TX =='lung')$TIME_TX_FU %>% IQR,
                     Alcohol_n = sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$Alcohol==T),
                     Alcohol_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$Alcohol==T) / 
                                           length(filter(cancer.comb.first, EVENT==0, TX =='lung')$Alcohol)),
                     COPD_n = sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$COPD==T),
                     COPD_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$COPD==T) / 
                                        length(filter(cancer.comb.first, EVENT==0, TX =='lung')$COPD)),
                     ILD_n = sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$ILD==T),
                     ILD_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$ILD==T) / 
                                       length(filter(cancer.comb.first, EVENT==0, TX =='lung')$ILD)),
                     Cirrhosis_n = sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$Cirrhosis==T),
                     Cirrhosis_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$Cirrhosis==T) / 
                                             length(filter(cancer.comb.first, EVENT==0, TX =='lung')$Cirrhosis)),
                     Hepatitis_n = sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$Hepatitis==T),
                     Hepatitis_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$Hepatitis==T) / 
                                             length(filter(cancer.comb.first, EVENT==0, TX =='lung')$Hepatitis)),
                     Bronchiectasis_n = sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$Bronchiectasis==T),
                     Bronchiectasis_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$Bronchiectasis==T) / 
                                                  length(filter(cancer.comb.first, EVENT==0, TX =='lung')$Bronchiectasis)),
                     Diabetes_n = sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$Diabetes==T),
                     Diabetes_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='lung')$Diabetes==T) / 
                                            length(filter(cancer.comb.first, EVENT==0, TX =='lung')$Diabetes))
)
lung.data <- data.frame(TX = 'lung', Group = c('all', 'cancer', 'no cancer'), rbind(all.llt, canc.llt, no.llt))

# liver tx
all.llt <- data.frame(male_sex_n = sum(filter(cancer.comb.first, TX =='liver')$SEX=='male'),
                      male_sex_perc = 100*(sum(filter(cancer.comb.first, TX =='liver')$SEX=='male') / 
                                             length(filter(cancer.comb.first, TX =='liver')$SEX)),
                      tx_age_median = filter(cancer.comb.first, TX =='liver')$EVENT_AGE.comb %>% median,
                      tx_age_iqr = filter(cancer.comb.first, TX =='liver')$EVENT_AGE.comb %>% IQR,
                      years_from_tx_median = filter(cancer.comb.first, TX =='liver')$TIME_TX_FU %>% median,
                      years_from_tx_iqr = filter(cancer.comb.first, TX =='liver')$TIME_TX_FU %>% IQR,
                      Alcohol_n = sum(filter(cancer.comb.first, TX =='liver')$Alcohol==T),
                      Alcohol_perc = 100*(sum(filter(cancer.comb.first, TX =='liver')$Alcohol==T) / 
                                            length(filter(cancer.comb.first, TX =='liver')$Alcohol)),
                      COPD_n = sum(filter(cancer.comb.first, TX =='liver')$COPD==T),
                      COPD_perc = 100*(sum(filter(cancer.comb.first, TX =='liver')$COPD==T) / 
                                         length(filter(cancer.comb.first, TX =='liver')$COPD)),
                      ILD_n = sum(filter(cancer.comb.first, TX =='liver')$ILD==T),
                      ILD_perc = 100*(sum(filter(cancer.comb.first, TX =='liver')$ILD==T) / 
                                        length(filter(cancer.comb.first, TX =='liver')$ILD)),
                      Cirrhosis_n = sum(filter(cancer.comb.first, TX =='liver')$Cirrhosis==T),
                      Cirrhosis_perc = 100*(sum(filter(cancer.comb.first, TX =='liver')$Cirrhosis==T) / 
                                              length(filter(cancer.comb.first, TX =='liver')$Cirrhosis)),
                      Hepatitis_n = sum(filter(cancer.comb.first, TX =='liver')$Hepatitis==T),
                      Hepatitis_perc = 100*(sum(filter(cancer.comb.first, TX =='liver')$Hepatitis==T) / 
                                              length(filter(cancer.comb.first, TX =='liver')$Hepatitis)),
                      Bronchiectasis_n = sum(filter(cancer.comb.first, TX =='liver')$Bronchiectasis==T),
                      Bronchiectasis_perc = 100*(sum(filter(cancer.comb.first, TX =='liver')$Bronchiectasis==T) / 
                                                   length(filter(cancer.comb.first, TX =='liver')$Bronchiectasis)),
                      Diabetes_n = sum(filter(cancer.comb.first, TX =='liver')$Diabetes==T),
                      Diabetes_perc = 100*(sum(filter(cancer.comb.first, TX =='liver')$Diabetes==T) / 
                                             length(filter(cancer.comb.first, TX =='liver')$Diabetes))
                      
)
canc.llt <- data.frame(male_sex_n = sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$SEX=='male'),
                       male_sex_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$SEX=='male') / 
                                              length(filter(cancer.comb.first, EVENT==1, TX =='liver')$SEX)),
                       tx_age_median = filter(cancer.comb.first, EVENT==1, TX =='liver')$EVENT_AGE.comb %>% median,
                       tx_age_iqr = filter(cancer.comb.first, EVENT==1, TX =='liver')$EVENT_AGE.comb %>% IQR,
                       years_from_tx_median = filter(cancer.comb.first, EVENT==1, TX =='liver')$TIME_TX_FU %>% median,
                       years_from_tx_iqr = filter(cancer.comb.first, EVENT==1, TX =='liver')$TIME_TX_FU %>% IQR,
                       Alcohol_n = sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$Alcohol==T),
                       Alcohol_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$Alcohol==T) / 
                                             length(filter(cancer.comb.first, EVENT==1, TX =='liver')$Alcohol)),
                       COPD_n = sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$COPD==T),
                       COPD_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$COPD==T) / 
                                          length(filter(cancer.comb.first, EVENT==1, TX =='liver')$COPD)),
                       ILD_n = sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$ILD==T),
                       ILD_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$ILD==T) / 
                                         length(filter(cancer.comb.first, EVENT==1, TX =='liver')$ILD)),
                       Cirrhosis_n = sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$Cirrhosis==T),
                       Cirrhosis_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$Cirrhosis==T) / 
                                               length(filter(cancer.comb.first, EVENT==1, TX =='liver')$Cirrhosis)),
                       Hepatitis_n = sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$Hepatitis==T),
                       Hepatitis_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$Hepatitis==T) / 
                                               length(filter(cancer.comb.first, EVENT==1, TX =='liver')$Hepatitis)),
                       Bronchiectasis_n = sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$Bronchiectasis==T),
                       Bronchiectasis_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$Bronchiectasis==T) / 
                                                    length(filter(cancer.comb.first, EVENT==1, TX =='liver')$Bronchiectasis)),
                       Diabetes_n = sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$Diabetes==T),
                       Diabetes_perc = 100*(sum(filter(cancer.comb.first, EVENT==1, TX =='liver')$Diabetes==T) / 
                                              length(filter(cancer.comb.first, EVENT==1, TX =='liver')$Diabetes))
)
no.llt <- data.frame(male_sex_n = sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$SEX=='male'),
                     male_sex_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$SEX=='male') / 
                                            length(filter(cancer.comb.first, EVENT==0, TX =='liver')$SEX)),
                     tx_age_median = filter(cancer.comb.first, EVENT==0, TX =='liver')$EVENT_AGE.comb %>% median,
                     tx_age_iqr = filter(cancer.comb.first, EVENT==0, TX =='liver')$EVENT_AGE.comb %>% IQR,
                     years_from_tx_median = filter(cancer.comb.first, EVENT==0, TX =='liver')$TIME_TX_FU %>% median,
                     years_from_tx_iqr = filter(cancer.comb.first, EVENT==0, TX =='liver')$TIME_TX_FU %>% IQR,
                     Alcohol_n = sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$Alcohol==T),
                     Alcohol_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$Alcohol==T) / 
                                           length(filter(cancer.comb.first, EVENT==0, TX =='liver')$Alcohol)),
                     COPD_n = sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$COPD==T),
                     COPD_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$COPD==T) / 
                                        length(filter(cancer.comb.first, EVENT==0, TX =='liver')$COPD)),
                     ILD_n = sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$ILD==T),
                     ILD_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$ILD==T) / 
                                       length(filter(cancer.comb.first, EVENT==0, TX =='liver')$ILD)),
                     Cirrhosis_n = sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$Cirrhosis==T),
                     Cirrhosis_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$Cirrhosis==T) / 
                                             length(filter(cancer.comb.first, EVENT==0, TX =='liver')$Cirrhosis)),
                     Hepatitis_n = sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$Hepatitis==T),
                     Hepatitis_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$Hepatitis==T) / 
                                             length(filter(cancer.comb.first, EVENT==0, TX =='liver')$Hepatitis)),
                     Bronchiectasis_n = sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$Bronchiectasis==T),
                     Bronchiectasis_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$Bronchiectasis==T) / 
                                                  length(filter(cancer.comb.first, EVENT==0, TX =='liver')$Bronchiectasis)),
                     Diabetes_n = sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$Diabetes==T),
                     Diabetes_perc = 100*(sum(filter(cancer.comb.first, EVENT==0, TX =='liver')$Diabetes==T) / 
                                            length(filter(cancer.comb.first, EVENT==0, TX =='liver')$Diabetes))
)
liver.data <- data.frame(TX = 'liver', Group = c('all', 'cancer', 'no cancer'), rbind(all.llt, canc.llt, no.llt))

rbind(lung.data, liver.data) %>% fwrite(., 'results/kidney_cancer/LLT_table.tsv', sep = '\t')

##
tx.longi$FINNGENID %>% unique %>% map(function(x) {
  tmp <- tx.longi %>% filter(FINNGENID == x) %>% .$ENDPOINT
  tmp[grepl('hepatitis|cirrhosis|cholangitis|atresia|alcohol', tmp, ignore.case = T)] %>% unique
}) %>% unlist %>% table %>% sort %>% tail(25)

tx.longi$FINNGENID %>% unique %>% map(function(x) {
  tmp <- tx.longi %>% filter(FINNGENID == x) %>% .$ENDPOINT
  tmp[grepl('copd|ild|interstitial|obstructive|cystic|pulmonary|fibrosis|emphysema', tmp, ignore.case = T)] %>% unique
}) %>% unlist %>% table %>% sort %>% tail(25)



