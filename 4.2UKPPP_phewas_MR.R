



rm(list = ls())
library(TwoSampleMR)
library(tidyverse)
library(data.table)
setwd("M:/pQTL_BC_MR/")
load('UKBPPP_as_expo_0.1_100kb/ALPI_P09923_OID30497_v1_Inflammation_II.rdata')

setwd("H:/2data/GWAS data/FinnR9_all_outcome_rdata/")
outfile <- list.files(pattern = 'rdata')
AllOR <- data.frame()
for (i in outfile) {
  setwd("H:/2data/GWAS data/FinnR9_all_outcome_rdata")

  print(which(outfile==i))
  load(i)
  trait <-sapply(strsplit(i,split = '_finndata'), '[',1)  
  
  out <- finndata %>% dplyr::filter(SNP %in% exp$SNP)
  if (is.null(out)) {
    next              }
  else{
    repeat{dat <- try(TwoSampleMR::harmonise_data(exposure_dat = exp, outcome_dat = out,
                                                  action = 2))
    if(!('try-error' %in% class(dat))){
      break }}
    dat <- subset(dat, dat$mr_keep == TRUE)
    if (nrow(dat) == 0) {
      next            }
    else{
      res <- GagnonMR::primary_MR_analysis(dat = dat)
      res$exposure <- 'ALPI'
      res$outcome <- trait
      OR <- generate_odds_ratios(res)
      
      if (OR$method == 'Inverse variance weighted') {
        het <- mr_heterogeneity(dat,method_list = 'mr_ivw')
        pleio <- mr_pleiotropy_test(dat)
        direct <- directionality_test(dat)
        OR <- cbind(OR,het_Qval=het$Q_pval,pleio_pval=pleio$pval,
                    direction=direct$correct_causal_direction)
        AllOR <- rbind(AllOR,OR)
      }
      else{
        direct <- directionality_test(dat)
        OR <- cbind(OR,het_Qval=NA,pleio_pval=NA,
                    direction=direct$correct_causal_direction)
        AllOR <- rbind(AllOR,OR)
      }
    }
  }
}


setwd("M:/pQTL_BC_MR/")

AllOR$BFR <- p.adjust(AllOR$pval,method = "bonferroni")
save(AllOR,file = '4.2UKPPP_ALPI_phewas_MR.rdata')
write.csv(AllOR,file = '4.2UKPPP_ALPI_phewas_MR.csv')

AllOR_NOhete1 <- subset(AllOR,BFR<0.05 & het_Qval>0.05 & pleio_pval>0.05 & direction==TRUE) 
AllOR_NOhete2 <- subset(AllOR,BFR<0.05 & is.na(het_Qval) & direction==TRUE)
AllOR_NOhete <- rbind(AllOR_NOhete1,AllOR_NOhete2)
save(AllOR_NOhete,file = '4.2UKPPP_ALPI_phewas_MR_BFR_NOhete.rdata')
