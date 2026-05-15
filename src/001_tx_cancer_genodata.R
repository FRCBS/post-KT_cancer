
library(tidyverse)
library(data.table)


## ========================================================
##
## FinnGen kidney tx cancer genotype data preparation
##
## ========================================================


## --------------------------------------------------------
## Target SNP data extraction
## --------------------------------------------------------

# FG variants
bim     <- fread('/finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12.bim')
bim     <- data.frame(bim, str_split_fixed(bim$V2, '_', 3))
bim     <- unite(bim, Var, X1:X2, sep='_')
bim$Var <- gsub('chr', '', bim$Var)

# pancancer variants
panc <- fread('data/pancancer_filtered.tsv')
filter(bim, Var %in% panc$ID_hg38) %>% dim
fg.panc <- filter(bim, Var %in% panc$ID_hg38)
fwrite(fg.panc, 'data/FG_pancancer.tsv', sep='\t')
write(fg.panc$V2, 'data/pancSNP.list', ncolumns=1)

# run plink to extract pancancer variants
system(paste0("plink2 --bfile /finngen/library-red/finngen_R12/genotype_plink_2.0/data/finngen_R12 ",
              "--extract data/pancSNP.list ",
              "--make-bed --out data/genotypes/tx_pancancer"))

# split into 2 and create dosages
system(paste0("plink2 --bfile data/genotypes/tx_pancancer ",
              "--keep data/FG_split_sample.list ",
              "--recode A --out data/genotypes/tx_pancancer_01"))
system(paste0("plink2 --bfile data/genotypes/tx_pancancer ",
              "--remove data/FG_split_sample.list ",
              "--recode A --out data/genotypes/tx_pancancer_02"))

# keep kidney tx samples
system(paste0("plink2 --bfile data/genotypes/tx_pancancer ",
              "--keep data/kidney_tx.list ",
              "--make-bed --out data/genotypes/tx_kidney_pancancer"))

# convert to dosage
system(paste0("plink2 --bfile data/genotypes/tx_kidney_pancancer ",
              "--extract data/pancSNP.list ",
              "--recode A --out data/genotypes/tx_kidney_pancancer"))



## --------------------------------------------------------
## SNP dosage data and risk alleles
## --------------------------------------------------------

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
## Compute risk scores in two parts
## --------------------------------------------------------

# parts 01 and 02, one after the other
# dos <- fread('data/genotypes/tx_pancancer_01.raw')[, -c(1, 3:6)]
dos <- fread('data/genotypes/tx_pancancer_02.raw')[, -c(1, 3:6)]

# continue
dos.vars <- colnames(dos)[-1]
dos.vars <- data.frame(DosVars=dos.vars, str_split_fixed(dos.vars, '_', 5))
dos.vars$X1 <- gsub('chr', '', dos.vars$X1)
dos.vars <- unite(dos.vars, 'DosID', X1, X2, sep='_', remove=F)
dos.vars <- left_join(dos.vars, fg.panc, by=c('DosID'='Var')) %>% 
  dplyr::select(-c(X1, X2, X4, X5))
dos.vars <- right_join(dos.vars, panc, by=c('DosID'='ID_hg38')) %>% 
  dplyr::select(-c(X3.y, V1, V3, V4))
colnames(dos.vars)[3] <- 'DosAllele'
head(dos.vars)
dim(dos.vars)

# number of variants per cancer type
can.n <- dos.vars %>% group_by(CAN) %>% summarise(n=n()) %>% filter(n>1)
fwrite(can.n, 'results/kidney_cancer/FG_cancer_SNP_num.tsv', sep = '\t', quote = F)
dos.vars <- filter(dos.vars, CAN %in% can.n$CAN)


## orient dosage data according to risk allele

isDosRisk <- map(unique(dos.vars$CAN), function(x) {
  # x <- 'kidney'
  print(x)
  o <- filter(dos.vars, CAN==x)
  col.keep <- match(o$DosVars, colnames(dos)) %>% na.omit
  o <- filter(o, DosVars %in% colnames(dos))
  col.keep <- match(o$DosVars, colnames(dos)) %>% na.omit
  
  if(length(col.keep) == 1) {
    d <- data.frame(dos)[, c(col.keep, col.keep)]
    d[, (o$DosAllele != o$RiskAllele)] <- 2-d[, (o$DosAllele != o$RiskAllele)]
    return(data.frame(IID=dos$IID, d[, 1]))
  } else {
    d <- data.frame(dos)[, col.keep]
    d[, (o$DosAllele != o$RiskAllele)] <- 2-d[, (o$DosAllele != o$RiskAllele)]
    return(data.frame(IID=dos$IID, d))
  }
  
})
names(isDosRisk) <- unique(dos.vars$CAN)

# check
filter(dos.vars, CAN=='bladder') %>% head
isDosRisk[['bladder']][1:5, 1:6]
filter(dos.vars, CAN=='kidney') %>% head


# compute risk scores for each indiv by averaging risk alleles per cancer type
risk.scores <- map2(1:length(isDosRisk), names(isDosRisk), function(i, y) {
  # i <- 13
  # y <- 'kidney'
  rowmeans2 <- function(x) {
    if(is.null(dim(x)) & length(x) > 1) x else rowMeans(x)
  }
    
  o <- data.frame(IID = isDosRisk[[i]]$IID, 
                  C = isDosRisk[[i]][, -1] %>% rowmeans2)
  colnames(o)[2] <- y
  return(o)                
}) %>% Reduce(inner_join, .)

risk.scores %>% head

# highest per individual
risk.scores$max <- risk.scores %>% dplyr::select(one_of(dos.vars$CAN %>% unique)) %>% scale %>% apply(., 1, max)
risk.scores$min <- risk.scores %>% dplyr::select(one_of(dos.vars$CAN %>% unique)) %>% scale %>% apply(., 1, min)
risk.scores$sum <- risk.scores %>% dplyr::select(one_of(dos.vars$CAN %>% unique)) %>% scale %>% apply(., 1, sum)
risk.scores$pos <- risk.scores %>% dplyr::select(one_of(dos.vars$CAN %>% unique)) %>% scale %>% apply(., 1, function(x) {
  sum(x[x>0])
})

# fwrite(risk.scores, 'data/FG_pancancer_scores_01.tsv', sep='\t')
fwrite(risk.scores, 'data/FG_pancancer_scores_02.tsv', sep='\t')

