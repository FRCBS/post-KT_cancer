
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
## FinnGen survival analysis for cancer SNP score
##
## ========================================================


## --------------------------------------------------------
## data read-in
## --------------------------------------------------------

# longitudinal event data
cancer.kidney.first    <- fread('data/phenotypes/cancer.kidney.first')
cancer.kidney.first.tx <- fread('data/phenotypes/cancer.kidney.first.tx')

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
## BMI and smoke info available
## --------------------------------------------------------

cancer.kidney.first %>%  dplyr::select(c(SMOKE2, BMI)) %>% na.omit %>% dim


## --------------------------------------------------------
## FinnGen without Tx: 
## validation of cancer scores
## --------------------------------------------------------

## cancer score cox survival models
## All adjustments, including 
## Z21_FAMILY_HISTORY_MALIG_NEOPLASM & DIAB

# prepare data for cox
# must have BMI and smoking, family history of malignancy
cancer.kidney.first.model <- cancer.kidney.first %>% 
  dplyr::select(c(FINNGENID, BL_YEAR, BL_AGE, cohort,
                  PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                  EVENT_OR_FU, EVENT,
                  SEX, BMI, SMOKE2)) %>% na.omit
cancer.kidney.first.model <- left_join(cancer.kidney.first.model, adj.phenos)
cancer.kidney.first.model <- left_join(cancer.kidney.first.model, risk.scores,
                                       by=c('FINNGENID'='IID')) %>% na.omit
# get IDs
cancer.kidney.first.model.ID <- cancer.kidney.first.model$FINNGENID
# numeric matrix
cancer.kidney.first.model <- model.matrix( ~. ,
                                           data = cancer.kidney.first.model %>% 
                                             dplyr::select(-FINNGENID)) %>% .[, -1] %>% data.frame()
# add IDs
cancer.kidney.first.model <- data.frame(FINNGENID=cancer.kidney.first.model.ID, 
                                        cancer.kidney.first.model)

# separate risk scores and covariates
risk.scores.2 <- cancer.kidney.first.model[, c(1, 54:69)] # extract risk scores
cancer.kidney.first.model <- cancer.kidney.first.model %>% dplyr::select(-c(49:51, 54:69))

# All risk scores:  
# fit cox on each risk score separately
fg.can.score <- map_dfr(2:ncol(risk.scores.2), function(i) {
  # i <- 2
  tmp <- inner_join(cancer.kidney.first.model, 
                    risk.scores.2[, c(1, i)]) %>% 
    dplyr::select(-FINNGENID)
  tmp[, ncol(tmp)] <- scale(tmp[, ncol(tmp)])
  res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., data = tmp)
  res <- tbl_regression(res.surv, exp=T)$table_body %>% data.frame() %>% 
    dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                    statistic, conf.low, conf.high, ci, p.value))
  
  data.frame(ScoreCan=colnames(risk.scores.2)[i], res) %>% return()
})

# save
fwrite(fg.can.score, 'results/kidney_cancer/FG_score_surv_res.tsv', sep='\t')

# save as html table
fg.can.score %>% kable() %>% kable_styling('striped', 'bordered') %>% 
  save_kable('results/kidney_cancer/fg.can.scores.html')

# check table in console
fread('results/kidney_cancer/FG_score_surv_res.tsv') %>% filter(variable %in% ScoreCan) %>% 
  kable(digits=30, format='simple')

# data for plotting
tmp <- fread('results/kidney_cancer/FG_score_surv_res.tsv') %>% filter(variable %in% ScoreCan) %>% 
  arrange(estimate)
tmp$variable <- factor(tmp$variable, levels=tmp$variable %>% unique)

# plot 
tmp %>% ggplot(aes(estimate, variable, xmax=conf.high, xmin=conf.low)) +
  geom_pointrange(shape=15) + 
  geom_vline(xintercept=1, linetype='dashed', linewidth=0.3) +
  coord_cartesian(xlim=c(0.99, 1.08), clip='off') +
  annotate(geom='text', parse=F, fontface='italic',
           label = paste('p =', filter(tmp, variable=='sum')$p.value %>% signif(3)),
           x=1.07, y=nrow(tmp)) +
  xlab('Hazard ratio') + ylab('Cancer risk score') +
  theme_survminer() +
  theme(panel.grid=element_blank(),
        legend.position='right',
        strip.text=element_blank(),
        strip.background=element_blank(),
        axis.text.y=element_text(size=6.5),
        axis.text.x=element_text(size=8),
        axis.title=element_text(size=10),
        axis.ticks=element_blank(),
        plot.margin=margin(8, 2, 10.5, 2),
        axis.line=element_line(linewidth=0.3),
        panel.spacing=unit(0, 'lines'))



## --------------------------------------------------------
## FinnGen without Tx: 
## validation of cancer scores
## --------------------------------------------------------

## cancer score cox survival models
## without BMI and SMOKE 
## Z21_FAMILY_HISTORY_MALIG_NEOPLASM & DIAB

# prepare data for cox
# no BMI and smoking, has family history of malignancy
cancer.kidney.first.model <- cancer.kidney.first %>% 
  dplyr::select(c(FINNGENID, BL_YEAR, BL_AGE, cohort,
                  PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                  EVENT_OR_FU, EVENT,
                  SEX)) %>% na.omit
cancer.kidney.first.model <- left_join(cancer.kidney.first.model, risk.scores,
                                       by=c('FINNGENID'='IID')) %>% na.omit

cancer.kidney.first.model.ID <- cancer.kidney.first.model$FINNGENID
cancer.kidney.first.model <- model.matrix( ~. ,
                                           data = cancer.kidney.first.model %>% 
                                             dplyr::select(-FINNGENID)) %>% .[, -1] %>% data.frame()
cancer.kidney.first.model <- data.frame(FINNGENID=cancer.kidney.first.model.ID, 
                                        cancer.kidney.first.model)

# separate risk scores and covariates
risk.scores.2 <- cancer.kidney.first.model[, c(1, 49:64)] # extract risk scores
cancer.kidney.first.model <- cancer.kidney.first.model %>% dplyr::select(-c(49:64))

# fit cox on each risk score separately
fg.can.score <- map_dfr(2:ncol(risk.scores.2), function(i) {
  tmp <- inner_join(cancer.kidney.first.model, 
                    risk.scores.2[, c(1, i)]) %>% 
    dplyr::select(-FINNGENID)
  tmp[, ncol(tmp)] <- scale(tmp[, ncol(tmp)])
  res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., data = tmp)
  res <- tbl_regression(res.surv, exp=T)$table_body %>% data.frame() %>% 
    dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                    statistic, conf.low, conf.high, ci, p.value))
  
  data.frame(ScoreCan=colnames(risk.scores.2)[i], res) %>% return()
})
# save
fwrite(fg.can.score, 'results/kidney_cancer/FG_score_surv_res_noadj.tsv', sep='\t')

# save as html table
fg.can.score %>% kable() %>% kable_styling('striped', 'bordered') %>% 
  save_kable('results/kidney_cancer/fg.can.scores_noadj.html')

# check table in console
fread('results/kidney_cancer/FG_score_surv_res_noadj.tsv') %>% kable(digits=30, format='simple')


## --------------------------------------------------------
## Forest plot FG score coefficients from 
## adjusted/unadjusted models
## --------------------------------------------------------

# read results
tmp <- rbind(fread('results/kidney_cancer/FG_score_surv_res.tsv') %>% 
               filter(variable %in% ScoreCan) %>% 
               data.frame(Model='Adjusted', .), 
             fread('results/kidney_cancer/FG_score_surv_res_noadj.tsv') %>% 
               filter(variable %in% ScoreCan) %>% 
               data.frame(Model='Unadjusted', .))

# arrange by coefficient
tmp$variable <- factor(tmp$variable, 
                       levels=tmp %>% group_by(variable) %>% 
                         summarise(MM=mean(estimate)) %>% arrange(MM) %>% .$variable)

# plot
p.fg.tx.cores <- tmp %>% ggplot(aes(estimate, variable, xmax=conf.high, xmin=conf.low, color=Model)) +
  geom_pointrange(position=position_dodge(width=0.7), 
                  fatten=4.0, alpha=0.9, shape=18, linewidth=0.4) + 
  geom_vline(xintercept=1, linetype='dashed', linewidth=0.3) +
  coord_cartesian(xlim=c(0.99, 1.074), clip='off') +
  # annotate(geom='text', parse=F, fontface='italic',
  #          label = paste('p =', filter(tmp, variable=='sum', Model=='Adjusted')$p.value %>% signif(3)),
  #          x=1.07, y=(nrow(tmp)/2)-0.50, size=2.3) +
  # annotate(geom='text', parse=F, fontface='italic',
  #          label = paste('p =', filter(tmp, variable=='sum', Model=='Unadjusted')$p.value %>% signif(3)),
  #          x=1.07, y=(nrow(tmp)/2)+0.50, size=2.3) +
  annotate(geom='text', #parse=F, fontface='italic',
           label = substitute(paste(italic('p'), ' < 1e-20')),
           x=1.066, y=(nrow(tmp)/2), size=3.0) +
  xlab('Cancer hazard ratio') + ylab('Cancer risk score type') +
  scale_color_manual(values=c('#e84c4f', '#03396c') %>% rev, name = 'Adjusted for \nBMI and smoking',
                     labels = c('no' ,'yes') %>% rev) + # #e84c4f #379AAF #21908CFF #FFBF00 #440154FF
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
p.fg.tx.cores



## --------------------------------------------------------
## Survival in score extremes
## FinnGen without Tx: cancer score survival model
## All adjustments, 
## including Z21_FAMILY_HISTORY_MALIG_NEOPLASM & DIAB
## --------------------------------------------------------

# prepare data for cox
# must have BMI and smoking, family history of malignancy
cancer.kidney.first.model <- cancer.kidney.first %>% 
  dplyr::select(c(FINNGENID, BL_YEAR, BL_AGE, cohort,
                  PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                  EVENT_OR_FU, EVENT, SEX, BMI, SMOKE2)) %>% na.omit
cancer.kidney.first.model <- left_join(cancer.kidney.first.model, adj.phenos[, c(1, 5:6)]) %>% na.omit
cancer.kidney.first.model.ID <- cancer.kidney.first.model$FINNGENID
cancer.kidney.first.model <- model.matrix( ~. ,
                                           data = cancer.kidney.first.model %>% 
                                             dplyr::select(-FINNGENID)) %>% .[, -1] %>% data.frame()
cancer.kidney.first.model <- data.frame(FINNGENID=cancer.kidney.first.model.ID, 
                                        cancer.kidney.first.model)
cancer.kidney.first.model <- inner_join(cancer.kidney.first.model, risk.scores[, c('IID', 'sum')],
                                        by=c('FINNGENID'='IID')) %>% na.omit
cancer.kidney.first.model <- cancer.kidney.first.model %>% dplyr::select(-FINNGENID)

ind <- which(colnames(cancer.kidney.first.model) %in% c('EVENT_OR_FU', 'EVENT'))
cancer.kidney.first.model[, -c(ind)] <- scale(cancer.kidney.first.model[, -c(ind)])

# fit cox
res.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., 
                  data = cancer.kidney.first.model)
res.surv %>% summary %>% coef
res.surv.tbl <- tbl_regression(res.surv, exp=T)$table_body %>% data.frame() %>% 
  dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                  statistic, conf.low, conf.high, ci, p.value))

# selected covariates for plotting
res.surv.tbl.2 <- res.surv.tbl %>% filter(variable %in% 
                                            c('SEXmale', 'BMI', 'SMOKE2yes', 'E4_DIABETES',
                                              'Z21_FAMILY_HISTORY_MALIG_NEOPLASM')) 
res.surv.tbl.2$variable[res.surv.tbl.2$variable=='SEXmale'] <- 'sex (male)'
res.surv.tbl.2$variable[res.surv.tbl.2$variable=='E4_DIABETES'] <- 'diabetes'
res.surv.tbl.2$variable[res.surv.tbl.2$variable=='SMOKE2yes'] <- 'smoking'
res.surv.tbl.2$variable[res.surv.tbl.2$variable=='Z21_FAMILY_HISTORY_MALIG_NEOPLASM'] <- 
  'family history of\nmalignant neoplasm'
res.surv.tbl.2 <- res.surv.tbl.2 %>% arrange(estimate)
res.surv.tbl.2$variable <- factor(res.surv.tbl.2$variable, 
                                  levels=res.surv.tbl.2$variable)

# covariate forest plot
p.res.surv.tbl.2 <- res.surv.tbl.2 %>% 
  ggplot(aes(estimate, variable, xmax=conf.high, xmin=conf.low)) +
  geom_pointrange(position=position_dodge(width=0.7), 
                  fatten=4.0, alpha=0.9, shape=18, linewidth=0.4) + 
  geom_vline(xintercept=1, linetype='dashed', linewidth=0.3) +
  xlab('Cancer hazard ratio') + ylab('Covariates') +
  guides(colour = guide_legend(override.aes = list(size=.2))) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.75, 0.3),
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.72, 'lines'),
        legend.box.background = element_rect(fill=NULL, linewidth = 0.2),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        axis.ticks = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        panel.spacing = unit(0, 'lines'))

pdf('results/kidney_cancer/FG_score_covariates.pdf', width=3.9, height=2.1)
p.res.surv.tbl.2
dev.off()


# fit sum risk score extremes 1-3 SD cox
tmp <- cancer.kidney.first.model #%>% dplyr::select(-FINNGENID)
# scale sum score
tmp[, ncol(tmp)] <- scale(tmp[, ncol(tmp)])

# plot density of sum score
p.fg.sum.histo <- ggplot(tmp, aes(sum)) +
  geom_histogram(aes(y = ..density..), color = 1, fill = 'white', linewidth = 0.3) +
  geom_density(color = 4, fill = 4, alpha = 0.25, linewidth = 0.3) +
  geom_density(aes(sum, ifelse(x > 2, after_stat(density), NA)), color = 'red', fill = 'red', alpha = 0.5) +
  geom_density(aes(sum, ifelse(x < -2, after_stat(density), NA)), color = 'forestgreen', 
               fill = 'forestgreen', alpha = 0.5) +
  coord_cartesian(xlim = c(-4, 4)) +
  geom_vline(xintercept = 0, linewidth = 0.3, linetype = 'dashed') +
  xlab('PRS') + ylab('Density') +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = 'inside',
        legend.position.inside = c(0.75, 0.3),
        legend.text = element_text(size = 7.5),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.72, 'lines'),
        legend.box.background = element_rect(fill = NULL, linewidth = 0.2),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.title = element_text(size = 10, color = 'black'),
        axis.ticks = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        panel.spacing = unit(0, 'lines'))
p.fg.sum.histo

# SD bin score cox model
cancer.kidney.first.model.bin <- map(seq(-3, 3, by=0.5), function(x) {
  tt <- tmp
  tt$RiskGroup <- rep(NA, nrow(tt))
  
  if(x >= 0) {
    tt$RiskGroup[tt[, ncol(tt)-1] >  x] <- 1
    tt$RiskGroup[tt[, ncol(tt)-1] < -x] <- 0
  } else {
    tt$RiskGroup[tt[, ncol(tt)-1] <  x] <- 1
    tt$RiskGroup[tt[, ncol(tt)-1] > -x] <- 0
  }
  
  tt <- dplyr::select(tt, -sum)
  #tt$RiskGroup <- scale(tt$RiskGroup)
  ind <- which(colnames(tt) %in% c('EVENT_OR_FU', 'EVENT'))
  tt[, -c(ind)] <- scale(tt[, -c(ind)])
  tt.surv <- coxph(Surv(EVENT_OR_FU, EVENT)  ~ ., data = tt)
  tt.surv %>% tbl_regression(., exp=T) %>% .$table_body %>% data.frame() %>% 
    dplyr::select(c(variable, n_obs, n_event, estimate, std.error, 
                    statistic, conf.low, conf.high, ci, p.value))
})

cancer.kidney.first.model.bin.res <- map_dfr(cancer.kidney.first.model.bin, 
                                            function(x) data.frame(x[nrow(x), ])) %>% data.frame
cancer.kidney.first.model.bin.res$variable <- paste0('', seq(-3, 3, by=0.5),  ' SD')
cancer.kidney.first.model.bin.res$variable <- factor(cancer.kidney.first.model.bin.res$variable,
                                                     levels = cancer.kidney.first.model.bin.res$variable)
cancer.kidney.first.model.bin.res$col <- ifelse(cancer.kidney.first.model.bin.res$statistic>0, 'risk', 'notrisk')

# covariate forest plot
p.res.surv.tbl.3 <- cancer.kidney.first.model.bin.res %>% filter(n_event>1000, n_event<50000) %>% 
  ggplot(aes(estimate, variable, xmax=conf.high, xmin=conf.low, color = col)) +
  geom_pointrange(position=position_dodge(width=0.7), 
                  fatten=4.0, alpha=0.9, shape=18, linewidth=0.4) + 
  scale_color_manual(values = c('red', 'forestgreen') %>% rev) +
  coord_cartesian(xlim = c(0.8, 1.20)) +
  geom_vline(xintercept=1, linetype='dashed', linewidth=0.3) +
  xlab('Cancer hazard ratio') + ylab('Binary PRS') +
  guides(colour = guide_legend(override.aes = list(size=.2))) +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        legend.position='none',
        #legend.position.inside = c(0.75, 0.3),
        #legend.text=element_text(size=7.5),
        #legend.title=element_text(size=7.5),
        #legend.key.size=unit(0.72, 'lines'),
        legend.box.background = element_rect(fill=NULL, linewidth = 0.2),
        strip.text=element_blank(),
        strip.background=element_blank(),
        axis.text.y=element_text(size=8, color='black'),
        axis.text.x=element_text(size=8, color='black'),
        axis.title=element_text(size=10, color='black'),
        axis.ticks=element_blank(),
        axis.line=element_line(linewidth=0.3),
        panel.spacing=unit(0, 'lines'))
p.res.surv.tbl.3



## --------------------------------------------------------
## Plot composite fig
## --------------------------------------------------------

jpeg('results/kidney_cancer/FG_scores.jpg', width=7.8, height=4.5, res=1000, units = 'in')
p.fg.tx.cores + (p.fg.sum.histo / p.res.surv.tbl.3) +
  plot_layout(widths = c(1, 0.65)) + plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(face = 'bold'))
dev.off()







