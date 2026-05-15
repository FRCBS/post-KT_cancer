
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
cancer.kidney.first    <- fread2('data/phenotypes/cancer.kidney.first')
cancer.kidney.first.tx <- fread2('data/phenotypes/cancer.kidney.first.tx')

# pancancer SNP risk scores
risk.scores <- rbind(
  fread2('data/FG_pancancer_scores_01.tsv'),
  fread2('data/FG_pancancer_scores_02.tsv')
)
# add sum risk that omits
risk.scores$sum2 <- risk.scores %>% dplyr::select(-c(IID, cervix, leukemia, gastricesoph, bladder)) %>% 
  scale %>% apply(., 1, sum)

# adjustment phenotypes
adj.phenos <- fread2('data/phenotypes/TxCan_cox_adj.tsv')

# tx longitudinal endpoints
tx.longi <- fread2('data/phenotypes/TX_LONG_phe.tsv')
tx.longi.immsupp <- fread2('data/phenotypes/TX_LONG_phe_immsupp.tsv')

# tx melanoma vs other cancers
hasMelanoma <- right_join(fread2('data/phenotypes/hasMelanoma.tsv'), tx.longi.immsupp[, 1:2])[, -3]
hasMelanoma$Melanoma[is.na(hasMelanoma$Melanoma)] <- 'NoCancer' 

# cancers in tx patients
# all cancers per individual
phe <- fread2('data/phenotypes/TX_EVENT_cancers.tsv')

# finngen pheno descriptions
fgp <- read_xlsx('data/phenotypes/FINNGEN_ENDPOINTS_DF12_Final_2023-05-17_public.xlsx')[, 1:4]

# cancer diag grouping
cangroup <- fread2('data/phenotypes/cancer_grouping.tsv')

# dialysis times
dialysis.times <- fread('results/kidney_cancer/KT_dialysis_times.tsv', data.table = F)
dialysis.times$Unique_dialysis_entries_log <- (dialysis.times$Unique_dialysis_entries+1) %>% log


## --------------------------------------------------------
## BMI and smoke info available
## --------------------------------------------------------

cancer.kidney.first.tx %>%  dplyr::select(c(SMOKE2, BMI)) %>% na.omit %>% dim


## --------------------------------------------------------
## survival; genetic risks in tx
## --------------------------------------------------------

## risk scores in tx age <40
## no adj for BMI, SMO
# family history data

# data
cancer.kidney.first.model <- cancer.kidney.first.tx %>% filter(EVENT_AGE.kidney<120) %>% 
  dplyr::select(c(FINNGENID, BL_AGE, BL_YEAR, EVENT_AGE.kidney,
                  PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                  EVENT_OR_FU, EVENT, # TIME_TX_FU
                  SEX))#, SMOKE2, BMI))

# join genetic risk score data
cancer.kidney.first.model <- inner_join(cancer.kidney.first.model, risk.scores[, c('IID', 'sum')], 
                                        by=c('FINNGENID'='IID'))
# family history data
cancer.kidney.first.model <- left_join(cancer.kidney.first.model, adj.phenos[, c(1, 5:6)])
# diag data
cancer.kidney.first.model.diag <- left_join(cancer.kidney.first.model, phe)
# dialysis time data
cancer.kidney.first.model <- left_join(cancer.kidney.first.model, dialysis.times[, c(1, 4)])
# remove ID from orig
cancer.kidney.first.model <- cancer.kidney.first.model %>% dplyr::select(-FINNGENID)
# NA = 0 in fam history data
cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM[
  is.na(cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM)] <- 0
cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM %>% sum

# check
cancer.kidney.first.model %>% dim
cancer.kidney.first.model %>% colnames

# format variable matrix to numeric
cancer.kidney.first.model <- model.matrix( ~. , data = cancer.kidney.first.model) %>% 
  .[, -1] %>% data.frame()
cancer.kidney.first.model %>% dim


# cox for tx age limits
cox.agelim <- map(seq(20, 80, by = 10), function(a) {
  res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., 
                    data = cancer.kidney.first.model %>% 
                      filter(EVENT_AGE.kidney < a) %>% 
                      dplyr::select(-Unique_dialysis_entries_log))
  # res.surv %>% tidy
             tbl_regression(res.surv, exp=T)$table_body %>% data.frame %>% 
               dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                               statistic, conf.low, conf.high, p.value)) %>% signif2
})
names(cox.agelim) <- seq(20, 80, by = 10)
map(cox.agelim, function(x) x$p.value[15]) %>% unlist %>% p.adjust(., method = 'BH')
map2_dfr(cox.agelim, cox.agelim %>% names, function(x, y) {
  data.frame(agelim = y, x)
}) %>% fwrite('results/kidney_cancer/Agelimits_KT_cox.tsv', sep = '\t')


# collect p-values from KT cox
pvalues <- map2_dfr(cox.agelim, cox.agelim %>% names, function(x, y) {
  data.frame(Variable = c(paste0('PRS KT agelimit = ', y), 
                          paste0('Family history of malign neoplasm KT agelimit = ', y)), 
             p.value = c(x[15, 'p.value'], x[16, 'p.value']))  
})
# add other tests
pvalues %>% rbind(., 
                  data.frame(Variable = c('FinnGen non-TX PRS sum', 'FinnGen non-TX PRS pos',
                                          'FinnGen non-TX PRS max', 'FinnGen non-TX PRS min',
                                          'Agelimit slope in KT', 'SOT non-KT PRSxDiab', 'SOT non-KT PRS',
                                          'KT vs FinnGen <40y', 'Age at KT'),
                             p.value = c(.650937e-50, 1.240482e-43, 3.538661e-12, 7.345429e-11,
                                         2.66e-18, 0.01, 0.207,
                                         0.00175, 4.4e-07))) %>% 
  mutate(FDR_BY = p.adjust(p.value, method = 'BY')) %>% 
  fwrite('results/kidney_cancer/p_values_adjusted.tsv', sep = '\t')
  
  


# cox model for <40y
res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., 
                  data = cancer.kidney.first.model %>% 
                    filter(EVENT_AGE.kidney<40) %>% 
                    dplyr::select(-Unique_dialysis_entries_log))
res.surv %>% summary %>% coef
res.surv %>% tidy
fwrite(tbl_regression(res.surv, exp=T)$table_body %>% data.frame() %>% 
         dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                         statistic, conf.low, conf.high, p.value)) %>% signif2,
       'results/kidney_cancer/cox_under40tx_noBMI.tsv', sep = '\t'
)

# cox model for <40y
# PRS interaction with diabetes
res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ BL_AGE + BL_YEAR + EVENT_AGE.kidney + 
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEXmale + 
                    sum + sum*E4_DIABETES + Z21_FAMILY_HISTORY_MALIG_NEOPLASM + E4_DIABETES, 
                  data = cancer.kidney.first.model %>% 
                    filter(EVENT_AGE.kidney<40) %>% 
                    dplyr::select(-Unique_dialysis_entries_log))
res.surv %>% summary %>% coef
res.surv %>% tidy
fwrite(tbl_regression(res.surv, exp=T)$table_body %>% data.frame() %>% 
         dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                         statistic, conf.low, conf.high, p.value)) %>% signif2,
       'results/kidney_cancer/cox_under40tx_noBMI_diabXsum.tsv', sep = '\t'
)


# compute genetic score risk groups
cancer.kidney.first.model$HIGHRISK <- 0
cancer.kidney.first.model$HIGHRISK[(cancer.kidney.first.model$sum %>% scale) >  1] <-  1 
cancer.kidney.first.model$HIGHRISK[(cancer.kidney.first.model$sum %>% scale) < -1] <- -1 

# cox model for <40y for groups
res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., 
                  data = cancer.kidney.first.model %>% 
                    dplyr::select(-c(sum)) %>% filter(EVENT_AGE.kidney<40) %>% 
                    dplyr::select(-Unique_dialysis_entries_log))
res.surv %>% tidy

# risk group variables including family history
cancer.kidney.first.model$RiskGroup <- NA
cancer.kidney.first.model$RiskGroup2 <- NA
cancer.kidney.first.model$RiskGroup[cancer.kidney.first.model$HIGHRISK==  1] <- 'High PRS'
cancer.kidney.first.model$RiskGroup[cancer.kidney.first.model$HIGHRISK== -1] <- 'Low PRS'
cancer.kidney.first.model$RiskGroup2[
  cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM== 1] <- 'Family history of malignancy'
# set age KT age limit 40 for PRS group
cancer.kidney.first.model$RiskGroup[cancer.kidney.first.model$EVENT_AGE.kidney>39] <- NA
cancer.kidney.first.model$RiskGroup %>% table




## --------------------------------------------------------
## survival; risk group comparison
## adj for DIAB, sum, FAM
## --------------------------------------------------------

# group n
cancer.kidney.first.model$RiskGroup %>% table
cancer.kidney.first.model$RiskGroup2 %>% table

# cox
res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., 
                  data = cancer.kidney.first.model %>% filter(EVENT_AGE.kidney<40) %>% 
                    dplyr::select(-c(RiskGroup, RiskGroup2, HIGHRISK)))
res.surv %>% summary %>% coef

# plot KM
survfit(Surv(EVENT_OR_FU, EVENT)  ~ RiskGroup, 
        data = cancer.kidney.first.model) %>% surv_pvalue()
p.tx.scores <- list(
  PRS = survfit(Surv(EVENT_OR_FU, EVENT)  ~ RiskGroup, 
                data = cancer.kidney.first.model)
  # Z21 = survfit(Surv(EVENT_OR_FU, EVENT)  ~ RiskGroup2, 
  #               data = cancer.kidney.first.model)
) %>% ggsurvplot_combine(data = cancer.kidney.first.model,
                         legend.labs = c('High PRS (n = 60)', 'Low PRS (n = 64)'),
                         #'Family history of malignant neoplasm (n = 5)'),
                         legend.title = 'Genetic risk groups in KT', # in KT patients age <40
                         palette = c('#00BA38', '#F8766D') %>% rev,
                         risk.table=F, size=0.8, censor = F) %>% .[[1]] # censor.size=3.0

p.tx.scores <- p.tx.scores + 
  xlab('Age (years)') +
  annotate(geom = 'text', 
           # label = substitute(paste('coxph ', italic('p'), ' = 0.0362')),
           label = 'coxph\np-value = 0.036\nHR = 1.083\n95% CI = 1.0051-1.1670',
           size = 3.0, x = 5, y = 0.3, hjust = 0) +
  theme_minimal() +
  theme(legend.direction = 'vertical',
        legend.position = 'top',
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8),#element_blank(),
        legend.key.size = unit(0.72, 'lines'),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.title.x = element_text(size = 10, color = 'black'),
        axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.3),
        axis.line = element_line(linewidth = 0.3),
        panel.spacing = unit(0, 'lines'),
        panel.grid = element_blank())

p.tx.scores



## --------------------------------------------------------
## survival; age group comparisons
## adj for BMI, SMO, DIAB, sum, FAM
## --------------------------------------------------------

## risk scores in tx age groups

cancer.kidney.first.model <- rbind(
  cancer.kidney.first.tx %>% filter(EVENT_AGE.kidney<40),
  cancer.kidney.first.tx %>% filter(EVENT_AGE.kidney>40)) %>% 
  dplyr::select(c(FINNGENID, BL_AGE, BL_YEAR, EVENT_AGE.kidney,
                  PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                  EVENT_OR_FU, EVENT, # TIME_TX_FU
                  SEX, BMI, SMOKE2))
cancer.kidney.first.model$AGEYOUNG <- ifelse(cancer.kidney.first.model$EVENT_AGE.kidney<40, 1, 0)
cancer.kidney.first.model <- dplyr::select(cancer.kidney.first.model, -EVENT_AGE.kidney)

# join with other covars
cancer.kidney.first.model <- inner_join(cancer.kidney.first.model, risk.scores[, c('IID', 'sum')], 
                                        by=c('FINNGENID'='IID'))
# family history data
cancer.kidney.first.model <- left_join(cancer.kidney.first.model, adj.phenos[, c(1, 5:6)])
cancer.kidney.first.model <- cancer.kidney.first.model %>% dplyr::select(-FINNGENID)
# NA = 0 in fam history data
cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM[
  is.na(cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM)] <- 0

# numeric
cancer.kidney.first.model <- model.matrix( ~. , data = cancer.kidney.first.model) %>% 
  .[, -1] %>% data.frame()
cancer.kidney.first.model %>% dim
filter(cancer.kidney.first.model, AGEYOUNG==1) %>% dim
filter(cancer.kidney.first.model, AGEYOUNG==0) %>% dim

# cox model
res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., data = cancer.kidney.first.model)
res.surv %>% summary %>% coef

# KM plot
p.tx.ages <- survfit(Surv(EVENT_OR_FU, EVENT)  ~ AGEYOUNG, 
                     data = cancer.kidney.first.model %>% na.omit) %>% 
  ggsurvplot(legend.labs = c('< 40 (n = 175)', 
                             '> 40 (n = 513)') %>% rev,
             palette = c('#233047', 'cadetblue4'),
             legend.title = 'Patient age at KT',
             risk.table=F, size=0.8, censor = F) %>% .[[1]]

p.tx.ages <- p.tx.ages + xlab('Age (years)') +
  annotate(geom = 'text', 
           # label = substitute(paste('coxph ', italic('p'), ' = 4.43e-07')), 
           label = 'coxph\np-value = 4.4e-07\nHR = 4.34\n95% CI = 2.45-7.66',
           size = 3.0, x = 10, y = 0.3, hjust = 0) +
  theme_minimal() +
  theme(legend.direction = 'vertical',
        legend.position = 'top',
        # legend.position.inside = c(0.3, 0.2),
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8),#element_blank(),
        legend.key.size = unit(0.72, 'lines'),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.title.x = element_text(size = 10, color = 'black'),
        axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.3),
        axis.line = element_line(linewidth = 0.3),
        panel.spacing = unit(0, 'lines'),
        panel.grid = element_blank())
p.tx.ages



## --------------------------------------------------------
## survival; age group comparisons
## FG vs. TX
## --------------------------------------------------------

cancer.kidney.first.model <- rbind(
  cancer.kidney.first.tx %>% filter(BL_AGE<400) %>% 
    dplyr::select(c(FINNGENID, BL_AGE, BL_YEAR, HASTX, EVENT_AGE.kidney,
                    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                    EVENT_OR_FU, EVENT, SEX, BMI, SMOKE2)),
  cancer.kidney.first %>% filter(BL_AGE<400) %>% 
    dplyr::select(c(FINNGENID, BL_AGE, BL_YEAR, HASTX, EVENT_AGE.kidney,
                    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                    EVENT_OR_FU, EVENT, SEX, BMI, SMOKE2))
)

# check
# x11();cancer.kidney.first.model %>% filter(HASTX==1) %>% .$EVENT_AGE.kidney %>% hist
cancer.kidney.first.model <- dplyr::select(cancer.kidney.first.model, -EVENT_AGE.kidney) 

# join with other covars
cancer.kidney.first.model <- inner_join(cancer.kidney.first.model, risk.scores[, c('IID', 'sum')], 
                                        by=c('FINNGENID'='IID'))
# family history data
cancer.kidney.first.model <- left_join(cancer.kidney.first.model, adj.phenos[, c(1, 5:6)])
cancer.kidney.first.model <- cancer.kidney.first.model %>% dplyr::select(-FINNGENID)
# NA = 0 in fam history data
cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM[
  is.na(cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM)] <- 0

# numeric
cancer.kidney.first.model <- model.matrix( ~. , data = cancer.kidney.first.model) %>% 
  .[, -1] %>% data.frame()
cancer.kidney.first.model %>% dim
cancer.kidney.first.model %>% filter(HASTX==1) %>% dim
cancer.kidney.first.model %>% filter(HASTX==0) %>% dim

# cox model
res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., data = cancer.kidney.first.model)
res.surv %>% summary %>% coef

# KM plot
p.txfg.ages <- survfit(Surv(EVENT_OR_FU, EVENT)  ~ HASTX, 
                       data = cancer.kidney.first.model %>% na.omit) %>% 
  ggsurvplot(legend.labs = c('KT (n = 688)', 
                             'non-KT (n = 279,967)') %>% rev,
             legend.title = 'FinnGen cohorts adjusted for\nage, diabetes smoking and BMI',
             risk.table=F, size=0.8, censor = F,
             palette = c('#E7B800', 'cornflowerblue')) %>% .[[1]] #  '#2E9FDF'

p.txfg.ages <- p.txfg.ages + xlab('Age (years)') +
  # annotate(geom = 'text', 
  #          label = paste0(substitute(paste('coxph ', italic('p'), ' = 0.00175'))), 
  #          size = 3.0, x = 20, y = 0.3) +
  annotate(geom = 'text', 
           label = paste0('coxph \n p-value = 0.00175 \n', 
                          'HR = 1.27\n95% CI = 1.09-1.47'), 
           size = 3.0, x = 10, y = 0.3, hjust = 0) +
  ylab('Cancer-free survival probability') +
  theme_minimal() +
  theme(legend.direction = 'vertical',
        legend.position = 'top',
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8),#element_blank(),
        legend.key.size = unit(0.72, 'lines'),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        axis.ticks = element_line(linewidth = 0.3),
        axis.line = element_line(linewidth = 0.3),
        panel.spacing = unit(0, 'lines'),
        panel.grid = element_blank())


p.txfg.ages



## --------------------------------------------------------
## survival; genetic risks in tx
## loop over ages and adjustments
## --------------------------------------------------------

# compute survival HRs when adding older KTs
tx.ages.surv <- map(18:80, function(i) {
  cancer.kidney.first.model <- cancer.kidney.first.tx %>% filter(EVENT_AGE.kidney < i) %>% 
    dplyr::select(c(FINNGENID, BL_AGE, BL_YEAR, EVENT_AGE.kidney,
                    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                    EVENT_OR_FU, EVENT, 
                    SEX))#, SMOKE2, BMI))
  # join genetic risk score data
  cancer.kidney.first.model <- inner_join(cancer.kidney.first.model, risk.scores[, c('IID', 'sum')], 
                                          by=c('FINNGENID'='IID'))
  # family history data
  cancer.kidney.first.model <- left_join(cancer.kidney.first.model, adj.phenos[, c(1, 5:6)])
  cancer.kidney.first.model <- cancer.kidney.first.model %>% dplyr::select(-FINNGENID)
  # NA = 0 in fam history data
  cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM[
    is.na(cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM)] <- 0
  
  # format to numeric
  cancer.kidney.first.model <- model.matrix( ~. , data = cancer.kidney.first.model) %>% 
    .[, -1] %>% data.frame()
  
  # cox model
  # cancer.kidney.first.model <- dplyr::select(cancer.kidney.first.model, -c(EVENT_AGE.kidney))
  res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., data = cancer.kidney.first.model) %>% 
    summary %>% coef %>% data.frame
  data.frame(Var = rownames(res.surv),
             Age = i,
             Adjusted = c('no BMI & smoking'),
             res.surv) %>% return()
})

# compute survival HRs starting from young age group and proceeding to older by 10 year interval
tx.ages.surv.2 <- map(seq(8, 68, by = 10), function(i) {
  # i <- 80
  # filter data subset
  # don't include BL_AGE  BL_YEAR, EVENT_AGE.kidney, as it causes coxph to fail
  cancer.kidney.first.model <- cancer.kidney.first.tx %>% 
    filter(EVENT_AGE.kidney >= i, EVENT_AGE.kidney < i+10) %>% 
    dplyr::select(c(FINNGENID,
                    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                    EVENT_OR_FU, EVENT, 
                    SEX))#, SMOKE2, BMI))
  # join genetic risk score data
  print(i)
  cancer.kidney.first.model <- inner_join(cancer.kidney.first.model, risk.scores[, c('IID', 'sum')], 
                                          by=c('FINNGENID'='IID'))
  # family history data
  cancer.kidney.first.model <- left_join(cancer.kidney.first.model, adj.phenos[, c(1, 6)])
  cancer.kidney.first.model <- cancer.kidney.first.model %>% dplyr::select(-FINNGENID)
  # NA = 0 in fam history data
  # cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM[
  #   is.na(cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM)] <- 0
  
  # format to numeric
  cancer.kidney.first.model <- model.matrix( ~. , data = cancer.kidney.first.model) %>% 
    .[, -1] %>% data.frame()
  
  # cox model
  res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., data = cancer.kidney.first.model,
                    control = coxph.control(iter.max = 200)) %>% 
    summary %>% coef %>% data.frame
  data.frame(Var = rownames(res.surv),
             N = nrow(cancer.kidney.first.model),
             Age = paste0(i, '-', i+9),
             Adjusted = c('no BMI, smoking, Fam'),
             res.surv,
             Upper = exp(res.surv[, 'coef']+(1.96*res.surv[, 'se.coef.'])),
             Lower = exp(res.surv[, 'coef']-(1.96*res.surv[, 'se.coef.']))) %>% return()
}) %>% map_dfr(function(x) filter(x, Var=='sum'))
tx.ages.surv.2$Age <- factor(tx.ages.surv.2$Age, levels = tx.ages.surv.2$Age)

# plot
p.tx.ages.surv.2.2 <- tx.ages.surv.2 %>% ggplot(aes(Age, exp.coef., ymin=Lower, ymax=Upper)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype='dashed', linewidth=0.3) +
  #geom_smooth(aes(group = 1), alpha = 0.2, linewidth = 0.5, se = F) +
  xlab('Age at KT') + ylab('PRS cancer HR') +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.72, 'lines'),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        axis.ticks = element_blank(),
        axis.line = element_line(linewidth=0.3),
        panel.spacing = unit(0, 'lines'),
        plot.margin = margin(20, 20, 2, 20))

p.tx.ages.surv.2.2


# by sample size bins, n = 100
tx.ages.surv.3 <- map(cancer.kidney.first.tx %>% arrange(EVENT_AGE.kidney) %>% .$EVENT_AGE.kidney %>% cut_number(15) %>% levels, function(x) {
  # start.end <- c(0.227,25.6)
  # filter data subset
  # don't include BL_AGE  BL_YEAR, EVENT_AGE.kidney, as it causes coxph to fail
  print(x)
  start.end <- str_split_fixed(x, '\\)|\\(|\\[|\\]|,', 4) %>% as.numeric %>% na.omit
  print(start.end)
  cancer.kidney.first.model <- cancer.kidney.first.tx %>% 
    filter(EVENT_AGE.kidney >= start.end[1], EVENT_AGE.kidney < start.end[2]) %>% 
    dplyr::select(c(FINNGENID,
                    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                    EVENT_OR_FU, EVENT, 
                    SEX))#, SMOKE2, BMI))
  # join genetic risk score data
  cancer.kidney.first.model <- inner_join(cancer.kidney.first.model, risk.scores[, c('IID', 'sum')], 
                                          by=c('FINNGENID'='IID'))
  # family history data
  cancer.kidney.first.model <- left_join(cancer.kidney.first.model, adj.phenos[, c(1, 6)])
  cancer.kidney.first.model <- cancer.kidney.first.model %>% dplyr::select(-FINNGENID)
  
  # format to numeric
  cancer.kidney.first.model <- model.matrix( ~. , data = cancer.kidney.first.model) %>% 
    .[, -1] %>% data.frame()
  
  # cox model
  res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., data = cancer.kidney.first.model,
                    control = coxph.control(iter.max = 200)) %>% 
    summary %>% coef %>% data.frame
  data.frame(Var = rownames(res.surv),
             N = nrow(cancer.kidney.first.model),
             Age = paste0(start.end[1], ' - ', start.end[2]),
             AgeMean = (start.end %>% mean)/max(cancer.kidney.first.tx$EVENT_AGE.kidney),
             Adjusted = c('no BMI, smoking, Fam'),
             res.surv,
             Upper = exp(res.surv[, 'coef']+(1.96*res.surv[, 'se.coef.'])),
             Lower = exp(res.surv[, 'coef']-(1.96*res.surv[, 'se.coef.']))) %>% return()
}) %>% map_dfr(function(x) filter(x, Var=='sum'))

# plot
p.tx.ages.surv.3.1 <- tx.ages.surv.3 %>% ggplot(aes(Age, exp.coef., ymin=Lower, ymax=Upper)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype='dashed', linewidth=0.3) +
  geom_smooth(aes(group = 1), alpha = 0.2, linewidth = 0.5) +
  xlab('Age at KT') + ylab('PRS cancer HR') +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.72, 'lines'),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.text.x = element_text(size = 8, color = 'black', angle = 33, vjust = 1, hjust = 1),
        axis.title = element_text(size = 10, color = 'black'),
        axis.ticks = element_line(linewidth = 0.3),
        axis.line = element_line(linewidth=0.3),
        panel.spacing = unit(0, 'lines'),
        plot.margin = margin(2, 20, 2, 20))

p.tx.ages.surv.3.1


# plot of HRs over cumulative KT ages starting at 18y
p.tx.ages.surv.2 <- rbind(
  tx.ages.surv %>% map_dfr(function(x) filter(x, Var=='sum'))) %>% 
  ggplot(aes(Age, exp.coef.)) +
  geom_point(size=1.2) +
  xlab('Maximum age at KT') + ylab('PRS cancer HR') +
  geom_smooth(method = 'lm', se = F) +
  geom_hline(yintercept = 1, linetype='dashed', linewidth=0.3) +
  stat_cor(p.accuracy = 1e-80, size = 3, 
           label.x.npc = 0.5,
           label.y.npc = 0.95) +
  scale_color_manual(values = c('#65233E')) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.72, 'lines'),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        axis.ticks = element_blank(),
        axis.line = element_line(linewidth=0.3),
        panel.spacing = unit(0, 'lines'))

p.tx.ages.surv.2


## composite plot

jpeg('results/kidney_cancer/Tx_ages_scores.jpg', width=9, height=7, res=1000, units = 'in')
(p.txfg.ages + p.tx.ages + p.tx.scores) / 
  ((p.tx.ages.surv.2 + p.tx.ages.surv.2.2) + plot_layout(widths = c(1, 0.6))) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = 'bold'))
dev.off()




## --------------------------------------------------------
## cancer diagnoses in risk groups
## --------------------------------------------------------

cancer.kidney.first.model %>% dim
cancer.kidney.first.model.diag %>% dim
# remove redundant
cancer.kidney.first.model.diag <- dplyr::select(cancer.kidney.first.model.diag, 
                                                -c('CD2_NEOPLASM', 'C3_CANCER_WIDE', 'C3_CANCER', 
                                                   'C3_SKIN', 'C3_OTHER_SKIN', 'CD2_INSITU', 
                                                   'CD2_INSITU_SKIN',
                                                   "HLP_C3_NSCLC_SQUAM",                         
                                                   "HLP_C3_NSCLC_ADENO",                        
                                                   "HLP_C3_SCLC",                                
                                                   "HLP_C3_ALL",                                 
                                                   "HLP_C3_CLL",                                 
                                                   "HLP_C3_AML",                                 
                                                   "HLP_C3_DLBCL",                               
                                                   "HLP_C3_BURKITT"))

#
cancer.kidney.first.model.diag.long <- cancer.kidney.first.model.diag %>% pivot_longer(21:269)

cancer.kidney.first.model.diag.long %>% filter(grepl('MYELOMA|LYMPHOMA|LEUKEMIA|LEUKAEMIA', name)) %>% 
  .$FINNGENID %>% unique %>% length
cancer.kidney.first.model.diag.long %>% filter(grepl('CARCINOMA|MELANOMA', name)) %>% 
  .$FINNGENID %>% unique %>% length



(grepl('MYELOMA|LYMPHOMA|LEUKEMIA|LEUKAEMIA', cancer.kidney.first.model.diag.long$name))
(grepl('CARCINOMA|MELANOMA', cancer.kidney.first.model.diag.long$name))


# count cancer types in high and low PRS risk groups 
# high group, n=60
# age <40
cancer.kidney.first.model.diag.h <- cancer.kidney.first.model.diag[
  cancer.kidney.first.model$HIGHRISK == 1 &
    !is.na(cancer.kidney.first.model$RiskGroup), 
  21:ncol(cancer.kidney.first.model.diag)]
cancer.kidney.first.model.diag.h %>% apply(., 1, function(x) {
  x[!is.na(x) & x!=0] %>% names
})

# low group, n=64
# age <40
cancer.kidney.first.model.diag.l <- cancer.kidney.first.model.diag[
  cancer.kidney.first.model$HIGHRISK == -1 &
    !is.na(cancer.kidney.first.model$RiskGroup), 
  21:ncol(cancer.kidney.first.model.diag)]
cancer.kidney.first.model.diag.l %>% apply(., 1, function(x) {
  x[!is.na(x) & x!=0] %>% names
})

map(1:nrow(cancer.kidney.first.model.diag.h), function(i) {
  x <- cancer.kidney.first.model.diag.h[i, ]
  z <- x[1:length(x)]
  z <- names(z)[!is.na(z) & z != 0]
  any(grepl('_ALL|_CLL_|_AML|LYMPH|LYMPHOMA|MYELOMA|LEUKAEMIA|LEUKEMIA|
                   HODGKIN|IMMUNO|DLBCL|BURKITT|LYMPHOID', z))
}) %>% unlist
map(1:nrow(cancer.kidney.first.model.diag.l), function(i) {
  x <- cancer.kidney.first.model.diag.l[i, ]
  z <- x[1:length(x)]
  z <- names(z)[!is.na(z) & z != 0]
  any(grepl('_ALL|_CLL_|_AML|LYMPH|LYMPHOMA|MYELOMA|LEUKAEMIA|LEUKEMIA|
                   HODGKIN|IMMUNO|DLBCL|BURKITT|LYMPHOID', z))
}) %>% unlist



# permute samples 100 times for cumulative sum calculation
permuteCumSum <- function(x, summary = T) {
  res <- map_dfc(1:100, function(i) {
    apply(x[sample(1:nrow(x), nrow(x), replace = F), ], 1, function(x) {
      res <- x %>% unlist %>% na.omit()
      res[res>0] %>% length
    }) %>% cumsum
  })
  if(summary) {
    data.frame(Mean = res %>% rowMeans(),
               Upper = apply(res, 1, function(k) quantile(k, 0.95)),
               Lower = apply(res, 1, function(k) quantile(k, 0.05)) %>% return()
    )
  } else {
    res %>% return()
  }
  
}
# permuteCumSum(cancer.kidney.first.model.diag.h)
# permuteCumSum(cancer.kidney.first.model.diag.h, summary = F)

# cumulative sum of indiv diagnoses
diags.cumsum <- rbind(
  data.frame(
    Group = 'High PRS',
    permuteCumSum(cancer.kidney.first.model.diag.h),
    Ind = 1:nrow(cancer.kidney.first.model.diag.h)
  ),
  data.frame(
    Group = 'Low PRS',
    permuteCumSum(cancer.kidney.first.model.diag.l),
    Ind = 1:nrow(cancer.kidney.first.model.diag.l)
  ))
diags.cumsum.2 <- rbind(
  data.frame(
    Group = 'High PRS',
    permuteCumSum(cancer.kidney.first.model.diag.h, F),
    Ind = 1:nrow(cancer.kidney.first.model.diag.h)
  ),
  data.frame(
    Group = 'Low PRS',
    permuteCumSum(cancer.kidney.first.model.diag.l, F),
    Ind = 1:nrow(cancer.kidney.first.model.diag.l)
  ))

# plot cumulative diagnoses in high and low PRS risk groups
p.diags.cumsum <- diags.cumsum %>% ggplot(
  aes(Ind, Mean, ymin = Lower, ymax = Upper, color = Group, fill = Group)) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(alpha = 0.2, linewidth = 0) +
  xlab('Sample index') + ylab('Cumulative sum of\ncancer diagnoses') +
  scale_color_manual(values = c('#00BA38', '#F8766D') %>% rev, name = '') +
  scale_fill_manual(values = c('#00BA38', '#F8766D') %>% rev, name = '') +
  theme_minimal() +
  theme(legend.direction = 'horizontal',
        legend.title.position = 'top',
        legend.position = 'top',
        legend.text = element_text(size = 7.5, hjust = 0.42),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.72, 'lines'),
        axis.text.y = element_text(size = 8, color='black'),
        axis.text.x = element_text(size = 8, color='black'),
        axis.title.x = element_text(size = 10, color='black'),
        axis.title.y = element_text(size = 10, color='black'),
        axis.ticks = element_line(linewidth = 0.3),
        axis.line = element_line(linewidth = 0.3),
        panel.grid = element_blank(),
        plot.margin = margin(1,1,0,1))
p.diags.cumsum


# same names in high and low?
(cancer.kidney.first.model.diag.l %>% colSums(na.rm = T) %>% names)==
  (cancer.kidney.first.model.diag.h %>% colSums(na.rm = T) %>% names)

# make a df of hi-lo diffrences
cancer.kidney.first.model.diag.hml <- 
  data.frame(Cancer = cancer.kidney.first.model.diag.h %>% colSums(na.rm = T) %>% names,
             HighmLow = cancer.kidney.first.model.diag.h %>% colSums(na.rm = T) - 
               cancer.kidney.first.model.diag.l %>% colSums(na.rm = T)) %>% 
  filter(HighmLow != 0) %>% arrange(HighmLow)
cancer.kidney.first.model.diag.hml$Cancer <- factor(cancer.kidney.first.model.diag.hml$Cancer,
                                                    levels = cancer.kidney.first.model.diag.hml$Cancer)
cancer.kidney.first.model.diag.hml <- left_join(cancer.kidney.first.model.diag.hml, fgp, 
                                                by = c('Cancer'='NAME'))
# format
cancer.kidney.first.model.diag.hml$LONGNAME <- gsub('_', '', cancer.kidney.first.model.diag.hml$LONGNAME)
cancer.kidney.first.model.diag.hml <- data.frame(cancer.kidney.first.model.diag.hml, 
                                                 cangroup[nrow(cangroup):1, ])
cancer.kidney.first.model.diag.hml <- cancer.kidney.first.model.diag.hml %>% group_by(new) %>% 
  summarise(SUM=sum(HighmLow)) %>% arrange(SUM)
cancer.kidney.first.model.diag.hml$new <- factor(cancer.kidney.first.model.diag.hml$new, 
                                                 levels = cancer.kidney.first.model.diag.hml$new)

# plot cancer difference in h and l risk <40 KT
p.cancer.kidney.first.model.diag.hml <- ggdotchart(cancer.kidney.first.model.diag.hml %>% data.frame, 
                                                   'new', 'SUM', 
                                                   add = 'segments', rotate = T, 
                                                   sorting = 'desc', 
                                                   dot.size = 2.3,
                                                   color = 'SUM') +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 'dashed', color = 'black') +
  scale_color_viridis_c(labels = c('less', '', '', '', '', 'more'),
                        name = 'Cancer diagnoses\nin high PRS group\n') +
  ylab('Diffrence in counts\n(High PRS - Low PRS)') +
  theme_minimal() +
  theme(legend.direction = 'horizontal',
        legend.title.position = 'left',
        legend.position = 'top',
        # legend.position.inside = c(0.8, 0.1),
        legend.text = element_text(size = 7.5, hjust = 0.42),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.72, 'lines'),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 8, color='black'),
        axis.text.x = element_text(size = 8, color='black'),
        axis.title.x = element_text(size = 10, color='black'),
        axis.title.y = element_blank(),
        axis.ticks = element_line(linewidth = 0.3),
        axis.line = element_line(linewidth = 0.3),
        panel.spacing = unit(0, 'lines'),
        panel.grid = element_blank())

p.cancer.kidney.first.model.diag.hml


## composite plot

jpeg('results/kidney_cancer/Tx_risks_diags.jpg', width=9, height=4.2, res=1000, units = 'in')
p.cancer.kidney.first.model.diag.hml + 
  # (p.diags.cumsum + plot_spacer() + plot_layout(heights = c(1, 0.05))) +
  (p.diags.cumsum) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = 'bold'))
dev.off()




