
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



## --------------------------------------------------------
## tx patients' cancer diagnoses after KT
## --------------------------------------------------------

# cancers in tx patients
phe <- fread('data/phenotypes/TX_EVENT_cancers.tsv', data.table=F)

# finngen phenos
fgp <- read_xlsx('data/phenotypes/FINNGEN_ENDPOINTS_DF12_Final_2023-05-17_public.xlsx')[, 1:4]

# counts of different cancer types
phe <- data.frame(Endpoint = phe[, grepl('CD2_|C3_', colnames(phe))] %>% colnames,
                  Count = phe[, grepl('CD2_|C3_', colnames(phe))] %>% colSums(na.rm=T)) %>% 
  filter(Count>4, !(Endpoint %in% c('CD2_NEOPLASM', 'C3_CANCER_WIDE', 'C3_CANCER', 'C3_SKIN', 
                                    'C3_OTHER_SKIN', 'CD2_INSITU'))) %>% 
  arrange(Count)

phe <- left_join(phe, fgp, by = c('Endpoint'='NAME'))
# phe$Endpoint <- factor(phe$Endpoint, levels = phe$Endpoint)
phe$LONGNAME <- factor(phe$LONGNAME, levels = phe$LONGNAME)

p.tx.diag <- ggplot(phe, aes(Count, LONGNAME)) +
  geom_bar(stat='identity') +
  ylab('Cancer type') +
  xlab('Count (min 5)') +
  coord_cartesian(expand = 0) +
  theme_minimal() +
  theme(axis.text.y=element_text(size=7),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())   
p.tx.diag

jpeg('results/kidney_cancer/Tx_diagnoses.jpg', width=7, height=7, res=600, units='in')
#pdf('results/kidney_cancer/Tx_cancer_diagnoses.pdf', width=6, height=7)
p.tx.diag
dev.off()



## --------------------------------------------------------
## tx patients' primary diagnoses prior to tx
## --------------------------------------------------------

# FG endpoints definitions and chapters
eps <- fread('data/phenotypes/FinnGen_def.tsv')
eps$code <- gsub('#', '',  eps$code, fixed=T)

# primary diagnoses before tx
tx.longi.primdiag <- fread('data/phenotypes/TX_LONG_phe_primary.tsv')

# primary diags
tx.longi.primdiag$Primary <- rep(NA, nrow(tx.longi.primdiag)) 
tx.longi.primdiag$Primary[
  grepl('N14_GLOMER_NEPHRITIS', tx.longi.primdiag$ENDPOINT)] <- 
  'Glomerulonephritis'
tx.longi.primdiag$Primary[
  grepl('T1D|D1_ICD8|E4_DM1|T2D|DM2|DIAB', tx.longi.primdiag$ENDPOINT)] <- 
  'Diabetes'
tx.longi.primdiag$Primary[
  grepl('Q17_CYSTIC_KIDNEY_DISEA', tx.longi.primdiag$ENDPOINT)] <- 
  'Cystic kidney disease'
tx.longi.primdiag$Primary[
  grepl('HYPTENS|HYPERTENSION', tx.longi.primdiag$ENDPOINT)] <- 
  'Hypertension'
tx.longi.primdiag$Primary[
  grepl('N14_RENALTUB', tx.longi.primdiag$ENDPOINT)] <-
  'Renal tubulo-interstitial disease'
tx.longi.primdiag$Primary[
  grepl('N14_CHRONKIDNEYDIS', tx.longi.primdiag$ENDPOINT)] <- 
  'Chronic kidney disease'
tx.longi.primdiag$Primary[
  grepl('Q17_CONGEN_MALFO_DEFORMAT_CHROMOSOMAL_ABNORMALITI|Q17_OTHER_CONGEN_MALFO_KIDNEY|
        Q17_CONGEN_MALFO_URINARY_SYSTEM|Q17_OTHER_CONGEN_MALFO_URINARY_SYSTEM', 
        tx.longi.primdiag$ENDPOINT)] <- 'Congenital malformation of the urinary system'


# add event age
tx.longi.primdiag <- left_join(tx.longi.primdiag, 
                               filter(tx.longi.primdiag, 
                                      ENDPOINT=='Z21_KIDNEY_TRANSPLANT_STATUS') %>% 
                                 .[, c('FINNGENID', 'EVENT_AGE')] %>% 
                                 rename(., 'TX_AGE'='EVENT_AGE'))
tx.longi.primdiag$TX_AGE_2 <- round(tx.longi.primdiag$TX_AGE/5)*5

# write
fwrite(tx.longi.primdiag, 'data/phenotypes/TX_LONG_phe_primary_2.tsv', sep='\t')

# read tx longitudinal data
tx.longi.primdiag <- fread('data/phenotypes/TX_LONG_phe_primary_2.tsv')
tx.longi.primdiag$Primary[tx.longi.primdiag$Primary==''] <- NA 

tx.longi.primdiag %>% 
  filter(!is.na(Primary), TimeDiffY_fromTX<15, TimeDiffY_fromTX>0) %>% 
  .$Primary %>% table

# plot
p.tx.primdiag.year.2 <- tx.longi.primdiag %>% 
  filter(!is.na(Primary), TimeDiffY_fromTX<15, TimeDiffY_fromTX>0) %>%
  mutate(Years=round(TimeDiffY_fromTX/0.5)*0.5) %>% 
  ggplot(aes(Years, fill=Primary)) +
  geom_bar(color=NA, alpha=0.9) +
  #scale_fill_d3() +
  scale_fill_manual(values=pal_d3()(7)[c(2,1,6,3,5,4,7)], name='Endpoint') +
  #scale_y_continuous(labels=c('0', '100000', '200000', '300000', '400000')) +
  scale_x_reverse() +
  xlab('Years before transplantation') +
  ylab('Registry entries count') +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        legend.position='right',
        legend.text=element_text(size=7.3),
        legend.title=element_blank(),
        legend.key.size=unit(0.72, 'lines'))

p.tx.primdiag.2.legend <- get_legend(p.tx.primdiag.year.2) %>% as_ggplot()

p.tx.primdiag.year.2 <- tx.longi.primdiag %>% 
  filter(!is.na(Primary), TimeDiffY_fromTX<15, TimeDiffY_fromTX>0) %>%
  mutate(Years=round(TimeDiffY_fromTX/0.5)*0.5) %>% 
  ggplot(aes(Years, fill=Primary)) +
  geom_bar(color=NA, alpha=0.9) +
  #scale_fill_d3() +
  scale_fill_manual(values=pal_d3()(7)[c(2,1,6,3,5,4,7)], name='Endpoint') +
  #scale_y_continuous(labels=c('0', '100000', '200000', '300000', '400000')) +
  scale_x_reverse() +
  xlab('Years before KT') +
  ylab('# Registry entries') +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        panel.border=element_blank(),
        legend.position='none',
        axis.title=element_text(size=10),
        axis.text=element_text(size=8),
        axis.line=element_line(linewidth=0.3))


p.tx.primdiag.years.ages.2 <- rbind(
  data.frame(tx.longi.primdiag %>% filter(TX_AGE<35, !is.na(Primary), TimeDiffY_fromTX<10), Age=' <35'),
  data.frame(tx.longi.primdiag %>% filter(TX_AGE>55, !is.na(Primary), TimeDiffY_fromTX<10), Age=' >55')
) %>% 
  ggplot(aes(TimeDiffY_fromTX, Age, fill=Primary)) +
  geom_boxplot(outlier.shape=NA, coef=0, linewidth=0.2, color='white', width=1.1, alpha=0.9, varwidth=T, ) +
  #geom_violin(scale='count', linewidth=0.2, alpha=0.9) +
  stat_summary(fun='mean', geom='point', shape=23, size=1.7, stroke=0.2) +
  facet_wrap(vars(Primary), nrow=7) +
  scale_fill_manual(values=pal_d3()(7)[c(2,1,6,3,5,4,7)]) +
  scale_x_reverse() +
  coord_cartesian(xlim=c(4.2, 0)) +
  xlab('Years before KT') + ylab('KT age') +
  theme_minimal() +
  theme(panel.grid=element_blank(),
        legend.position='none',
        strip.text=element_blank(),
        strip.background=element_blank(),
        axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=8),
        axis.title=element_text(size=10),
        plot.margin=margin(8, 2, 10.5, 2),
        axis.line=element_line(linewidth=0.3),
        panel.spacing=unit(0, 'lines'))


jpeg('results/kidney_cancer/Tx_primary_diag_2.jpg', width=8, height=2.7, res=600, units='in') 
# pdf('results/kidney_cancer/Tx_primary_diag_2.pdf', width=8, height=2.7) 
(p.tx.primdiag.year.2 + inset_element(p.tx.primdiag.2.legend, left=0.17, bottom=0.5, right=0.58, top=0.95)) %>% 
  annotate_figure(fig.lab.pos='top.left', fig.lab='a', fig.lab.face='bold') + 
  # (p.tx.primdiag.years.ages.2 / plot_spacer() + plot_layout(heights=c(1, 0.001))) %>% 
  p.tx.primdiag.years.ages.2 %>% 
  annotate_figure(fig.lab.pos='top.left', fig.lab='b', fig.lab.face='bold') + 
  plot_layout(widths=c(1.9, 1), ncol=2) 
dev.off()

jpeg('results/kidney_cancer/Tx_primary_diag.jpg', width=6, height=2.7, res=600, units='in') 
# pdf('results/kidney_cancer/Tx_primary_diag_1.pdf', width=6, height=3) 
(p.tx.primdiag.year.2 + inset_element(p.tx.primdiag.2.legend, left=0.17, bottom=0.5, right=0.58, top=0.95))
dev.off()
# jpeg('results/kidney_cancer/Tx_primary_diag_2.jpg', width=8, height=2.7, res=600, units='in') 
# ggarrange(
#   (p.tx.primdiag.year.2 + inset_element(p.tx.primdiag.2.legend, left=0.2, bottom=0.5, right=0.6, top=1)) %>% 
#   annotate_figure(fig.lab.pos='top.left', fig.lab='A'),
# 
#   p.tx.primdiag.years.ages.2 %>% 
#   annotate_figure(fig.lab.pos='top.left', fig.lab='B'),
#   
#   ncol=2, align='hv', widths=c(1.9, 1)  
# )  
# dev.off()





