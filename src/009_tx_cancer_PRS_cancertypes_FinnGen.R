
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
## TX cancer survival analysis F_2023_043
##
## ========================================================


fread2 <- function(x) fread(x, data.table = F)


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


## --------------------------------------------------------
## BMI and smoke info available
## --------------------------------------------------------

cancer.kidney.first.tx %>%  dplyr::select(c(SMOKE2, BMI)) %>% na.omit %>% dim




## --------------------------------------------------------
## Cancer type-specific model in FG
## --------------------------------------------------------

fglongi1 <- fread2('data/tx_cancer_1stCancer.tsv')
fglongi1.n100 <- fglongi1$ENDPOINT_class01 %>% table %>% sort(decreasing = T)
fglongi1.n100 <- fglongi1.n100[fglongi1.n100 > 100] %>% names # n > 100

# loop over cancer endpoints
canspef.cox <- fglongi1.n100 %>% map(function(x) {
  # x <- 'C3_BREAST'
  print(x)
  tmp <- right_join(filter(fglongi1, ENDPOINT_class01 == x),
                    cancer.kidney.first,
                    by = 'FINNGENID')
  tmp <- filter(tmp, !is.na(EVENT))
  
  tmp <- filter(tmp, !(FINNGENID %in% filter(tmp, EVENT==1, is.na(ENDPOINT_class01))$FINNGENID))
  
  cancer.kidney.first.model <- tmp %>%  
    dplyr::select(c(FINNGENID, BL_AGE, BL_YEAR, 
                    PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                    EVENT_OR_FU, EVENT, SEX))#, BMI, SMOKE2))
  
  # join with other covars
  cancer.kidney.first.model <- inner_join(cancer.kidney.first.model, risk.scores[, c('IID', 'sum')], 
                                          by=c('FINNGENID'='IID'))
  
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
  cancer.kidney.first.model$EVENT %>% table
  
  # cox model
  # stage local vs. metastasized
  res.surv <- coxph(Surv(EVENT_OR_FU, EVENT) ~ ., 
                    data = cancer.kidney.first.model %>% dplyr::select(-SEXmale))
  data.frame(CAN = x,
             res.surv %>% summary %>% coef %>% .[13, ] %>% t)
}) %>% do.call(rbind, .)
canspef.cox$FDR <- canspef.cox$Pr...z.. %>% p.adjust(method = 'BH')

canspef.cox <- left_join(canspef.cox, 
                         fglongi1$ENDPOINT_class01 %>% table %>% data.frame,
                         by = c('CAN'='.'))
canspef.cox$Controls <- 380238
colnames(canspef.cox) <- c("CAN", "coef", "exp.coef.", "se.coef.",  "z-value", "p-value",  
                           "FDR", "Cases", "Controls")
fwrite(canspef.cox, 'results/kidney_cancer/Cancer_spedific_sum.tsv', sep = '\t')

formatCoxRes <- function(x) {
  x$HRCI95Upper <- exp(x$coef + 1.96*x[, "se.coef."])
  x$HRCI95Lower <- exp(x$coef - 1.96*x[, "se.coef."])
  x <- rename(x, 'HR' = 'exp.coef.')
  return(x)
}

canspef.cox <- fread2('results/kidney_cancer/Cancer_spedific_sum.tsv')
fg.phe <- read_xlsx(
  '/home/ivm/Documents/R12/data/phenotypes/FINNGEN_ENDPOINTS_DF12_Final_2023-05-17_public.xlsx')
canspef.cox <- left_join(canspef.cox, fg.phe[, c('NAME', 'LONGNAME')], by = c('CAN'='NAME'))

canspef.cox %>% filter(FDR<0.05) %>% formatCoxRes 
canspef.cox.p <- canspef.cox %>% formatCoxRes %>% filter(HRCI95Lower>1) %>% arrange(HR)
canspef.cox.p$LONGNAME <- factor(canspef.cox.p$LONGNAME, levels = canspef.cox.p$LONGNAME)


jpeg('results/kidney_cancer/Cancer_types_PRS_FG.jpg', width = 7, height = 5, res = 1000, units = 'in')
# forest plot of HRs
ggplot(canspef.cox.p, aes(HR, LONGNAME, xmin = HRCI95Lower, xmax = HRCI95Upper)) +
  geom_pointrange() +
  ylab('Cancer endpoint') +
  xlim(c(0.975, 1.11)) +
  geom_vline(xintercept = 1, linetype = 'dashed', linewidth = 0.3) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) +
# barplot of n_cases  
ggplot(canspef.cox.p, aes(Cases, LONGNAME)) +
  geom_col(alpha = 0.5) +
  xlab('# Cases') +
  scale_x_continuous(transform = 'log2', breaks = c(10, 100, 5000)) +
  theme_light() +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  plot_layout(widths = c(1, 0.5))
dev.off()

