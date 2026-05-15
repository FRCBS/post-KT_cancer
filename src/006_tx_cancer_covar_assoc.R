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
## covariance associations
##
## ========================================================


## functions
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
cancer.kidney.first.model <- cancer.kidney.first.tx %>% filter(EVENT_AGE.kidney<100) %>% 
  dplyr::select(c(FINNGENID, BL_AGE, BL_YEAR, EVENT_AGE.kidney,
                  PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                  EVENT_OR_FU, EVENT, # TIME_TX_FU
                  SEX, SMOKE2, BMI))

# join genetic risk score data
cancer.kidney.first.model <- inner_join(cancer.kidney.first.model, risk.scores[, c('IID', 'sum')], 
                                        by=c('FINNGENID'='IID'))
# family history data
cancer.kidney.first.model <- left_join(cancer.kidney.first.model, adj.phenos[, c(1, 5:6)])
# diag data
cancer.kidney.first.model.diag <- left_join(cancer.kidney.first.model, phe)
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

# cox model for all ages with BMI, SMOKE
res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., 
                  data = cancer.kidney.first.model %>% dplyr::select(-c(sum)))
res.surv %>% summary %>% coef

fwrite(tbl_regression(res.surv, exp=T)$table_body %>% data.frame() %>% 
         dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                         statistic, conf.low, conf.high, p.value)) %>% signif2,
       'results/kidney_cancer/cox_all.tsv', sep = '\t'
)


# cox model for under 40 with BMI, SMOKE
res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., 
                  data = cancer.kidney.first.model %>% dplyr::select(-c(sum)) %>% 
                    filter(EVENT_AGE.kidney < 40))
res.surv %>% summary %>% coef

fwrite(tbl_regression(res.surv, exp=T)$table_body %>% data.frame() %>% 
         dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                         statistic, conf.low, conf.high, p.value)) %>% signif2,
       'results/kidney_cancer/cox_under40tx.tsv', sep = '\t'
)



# plot
p.tx.covars.cox <- res %>% ggplot(aes(estimate, variable, xmax=conf.high, xmin=conf.low)) +
  geom_pointrange(position=position_dodge(width=0.7), 
                  fatten=4.0, alpha=0.9, shape=18, linewidth=0.4) + 
  geom_vline(xintercept=1, linetype='dashed', linewidth=0.3) +
  # coord_cartesian(xlim=c(0.99, 1.074), clip='off') +
  # annotate(geom='text', #parse=F, fontface='italic',
  #          label = substitute(paste(italic('p'), ' < 1e-20')),
  #          x=1.066, y=(nrow(tmp)/2), size=3.0) +
  xlab('Hazard ratio') + ylab('Cancer risk score') +
  xlim(c(0.4, 15)) +
  # scale_color_manual(values=c('#e84c4f', '#03396c') %>% rev, name = 'Adjusted for \nBMI and smoking',
  #                    labels = c('no' ,'yes') %>% rev) + # #e84c4f #379AAF #21908CFF #FFBF00 #440154FF
  guides(colour = guide_legend(override.aes = list(size=.2))) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        legend.position='inside',
        legend.position.inside = c(0.75, 0.3),
        legend.text=element_text(size=7.5),
        legend.title=element_text(size=7.5),
        legend.key.size=unit(0.72, 'lines'),
        legend.box.background = element_rect(fill=NULL, linewidth = 0.2),
        strip.text=element_blank(),
        strip.background=element_blank(),
        axis.text.y=element_text(size=8, color='black'),
        axis.text.x=element_text(size=8, color='black'),
        axis.title=element_text(size=10, color='black'),
        axis.ticks=element_blank(),
        axis.line=element_line(linewidth=0.3),
        panel.spacing=unit(0, 'lines'))
p.tx.covars.cox



## --------------------------------------------------------
## ## plot PC
## --------------------------------------------------------

rbind(
  data.frame(PC1 = cancer.kidney.first$PC1,
             PC2 = cancer.kidney.first$PC2,
             Group = 'NoTx'),
  data.frame(PC1 = cancer.kidney.first.tx$PC1,
             PC2 = cancer.kidney.first.tx$PC2,
             Group = 'KT')
) %>% ggplot(aes(PC1, PC2, color = Group)) +
  geom_density2d(bins = 15) +
  theme_minimal()

rbind(
  data.frame(PC1 = cancer.kidney.first$PC1,
             PC3 = cancer.kidney.first$PC3,
             Group = 'NoTx'),
  data.frame(PC1 = cancer.kidney.first.tx$PC1,
             PC3 = cancer.kidney.first.tx$PC3,
             Group = 'KT')
) %>% ggplot(aes(PC1, PC3, color = Group)) +
  geom_density2d(bins = 15) +
  theme_minimal()


## --------------------------------------------------------
## PRS target SNP data
## --------------------------------------------------------

panc <- fread('data/pancancer_filtered.tsv')
fg.panc <-  fread('data/FG_pancancer.tsv')
tmp <- inner_join(panc, fg.panc, by = c('ID_hg38'='Var'))

panc %>% group_by(CAN) %>% summarise(SUM=n())
tmp %>% group_by(CAN) %>% summarise(SUM=n())
snp.nums <- data.frame(panc %>% group_by(CAN) %>% summarise(SUM=n()),
                       tmp %>% group_by(CAN) %>% summarise(SUM=n()) %>% .[, 2])
colnames(snp.nums) <- c('Cancer', 'Original', 'FinnGen')
fwrite(snp.nums, 'results/kidney_cancer/Cancer_type_target_SNP_numbers.tsv', sep='\t')

# SNP data
panc <- fread('data/pancancer_filtered.tsv') %>% 
  dplyr::select(-c(LiftIn, V1, POS_hg38, ID_hg19))
fg.panc <- fread('data/FG_pancancer.tsv')
fg.panc.snps <- left_join(fg.panc, panc[, c('ID_hg38', 'ID', 'RiskAllele', 'CAN')], 
                          by = c('Var'='ID_hg38'))
fg.panc.snps <- fg.panc.snps[, c(1, 4, 9:11)]
colnames(fg.panc.snps) <- c('Chr', 'Position_hg38', 'rsID', 'Risk_allele', 'Cancer_type')
fg.panc.snps <- fg.panc.snps %>% arrange(Cancer_type)
fg.panc.snps$rsID[!grepl('^rs', fg.panc.snps$rsID)] <- NA
# write FG target SNP list
fwrite(fg.panc.snps, 'results/kidney_cancer/FinnGen_PRS_SNP_list.tsv', 
       quote = F, sep = '\t', na = 'NA')



## --------------------------------------------------------
## Fam hist of malignancy survival
## --------------------------------------------------------

cancer.kidney.first.model <- cancer.kidney.first %>% 
  dplyr::select(c(FINNGENID, BL_AGE, BL_YEAR, EVENT_AGE.kidney,
                  PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                  EVENT_OR_FU, EVENT, # TIME_TX_FU
                  SEX, SMOKE2, BMI))

# join genetic risk score data
cancer.kidney.first.model <- inner_join(cancer.kidney.first.model, risk.scores[, c('IID', 'sum')], 
                                        by=c('FINNGENID'='IID'))
# family history data
cancer.kidney.first.model <- left_join(cancer.kidney.first.model, adj.phenos[, c(1, 5:6)])
# diag data
cancer.kidney.first.model.diag <- left_join(cancer.kidney.first.model, phe)
# remove ID from orig
cancer.kidney.first.model <- cancer.kidney.first.model %>% dplyr::select(-c(FINNGENID, EVENT_AGE.kidney))

# format variable matrix to numeric
cancer.kidney.first.model <- model.matrix( ~. , data = cancer.kidney.first.model) %>% 
  .[, -1] %>% data.frame()
cancer.kidney.first.model %>% dim
cancer.kidney.first.model$Z21_FAMILY_HISTORY_MALIG_NEOPLASM %>% table
cancer.kidney.first.model$E4_DIABETES %>% table

# cox model
res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., 
                  data = cancer.kidney.first.model)
res.surv

# plot KM
p.fg.z21 <- list(
  PRS = survfit(Surv(EVENT_OR_FU, EVENT)  ~ Z21_FAMILY_HISTORY_MALIG_NEOPLASM, 
                data = cancer.kidney.first.model)
) %>% ggsurvplot_combine(data = cancer.kidney.first.model,
                         legend.labs = c('Family history of malignancy (n = 2,026)', 
                                         'Control (n = 238,691)') %>% rev,
                         legend.title = 'FinnGen nontransplantation cohorts', 
                         palette = c('#F8766D', 'grey30') %>% rev,
                         risk.table=F, size=0.8, censor = F) %>% .[[1]] # censor.size=3.0

p.fg.z21 <- p.fg.z21 + 
  xlab('Time (years)') +
  annotate(geom = 'text', label = substitute(paste(italic('p'), ' < 2e-16')), 
           size = 3.0, x = 15, y = 0.3) +
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
        axis.ticks = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        panel.spacing = unit(0, 'lines'),
        panel.grid = element_blank())

p.fg.z21


p.fg.diab <- list(
  PRS = survfit(Surv(EVENT_OR_FU, EVENT)  ~ E4_DIABETES, 
                data = cancer.kidney.first.model)
) %>% ggsurvplot_combine(data = cancer.kidney.first.model,
                         legend.labs = c('Diabetes (n = 38,997)', 
                                         'Control (n = 201,720)') %>% rev,
                         legend.title = 'FinnGen nontransplantation cohorts', 
                         palette = c('#F8766D', 'grey30') %>% rev,
                         risk.table=F, size=0.8, censor = F) %>% .[[1]] # censor.size=3.0

p.fg.diab <- p.fg.diab + 
  xlab('Time (years)') +
  annotate(geom = 'text', label = substitute(paste(italic('p'), ' = 0.0115')), 
           size = 3.0, x = 15, y = 0.3) +
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
        axis.ticks = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        panel.spacing = unit(0, 'lines'),
        panel.grid = element_blank())

p.fg.diab


jpeg('results/kidney_cancer/FG_covars_KM.jpg', width=8, height=5, res=1000, units = 'in')
(p.fg.z21 + p.fg.diab) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = 'bold'))
dev.off()





