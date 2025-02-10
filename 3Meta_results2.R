

rm(list = ls())
library(meta)
library(tidyverse)
library(export)
library(data.table)
setwd("M:/pQTL_BC_MR/")


out_plots <- function(filename,pic_width=5,pic_height=7){
  graph2png(file=filename,width=pic_width,height=pic_height)
  graph2pdf(file=filename,width=pic_width,height=pic_height)
}

load("M:/pQTL_BC_MR/all_pQTL_mrres_Finn_NOhete_FDR.rdata")
AllOR_Finn <- AllOR_NOhete_FDR

load("M:/pQTL_BC_MR/all_pQTL_mrres_BCAC_NOhete_FDR.rdata")
AllOR_BCAC <- AllOR_NOhete_FDR

load("M:/pQTL_BC_MR/all_pQTL_mrres_BCAC_iCOGS_NOhete_FDR.rdata")
AllOR_BCAC_iCOGS <- AllOR_NOhete_FDR

AllOR_Finn$short_expo <- sapply(strsplit(AllOR_Finn$exposure,split = '_'),'[',3) 
AllOR_BCAC$short_expo <- sapply(strsplit(AllOR_BCAC$exposure,split = '_'),'[',3) 
AllOR_BCAC_iCOGS$short_expo <- sapply(strsplit(AllOR_BCAC_iCOGS$exposure,split = '_'),'[',3) 

same_exp <- Reduce(intersect, list(AllOR_Finn$short_expo, AllOR_BCAC$short_expo,AllOR_BCAC_iCOGS$short_expo))

setwd('M:/pQTL_BC_MR/3Meta_results2/')
same_exp_metasig <- as.character()
for (i in same_exp) {

  sig_Finn1 <- AllOR_Finn[AllOR_Finn$short_expo==i,]
  sig_Finn1$study <- 'FinnGen'
  sig_Finn1 <- sig_Finn1 %>% dplyr::select(short_expo,study,b,se)
  
  sig_BCAC1 <- AllOR_BCAC[AllOR_BCAC$short_expo==i,]
  sig_BCAC1$study <- 'BCAC_OncoArray'
  sig_BCAC1 <- sig_BCAC1 %>% dplyr::select(short_expo,study,b,se)
  
  sig_BCAC1_iCOGS <- AllOR_BCAC_iCOGS[AllOR_BCAC_iCOGS$short_expo==i,]
  sig_BCAC1_iCOGS$study <- 'BCAC_iCOGS'
  sig_BCAC1_iCOGS <- sig_BCAC1_iCOGS %>% dplyr::select(short_expo,study,b,se)
  
  sig_merge <- rbind(sig_Finn1,sig_BCAC1,sig_BCAC1_iCOGS)

  sig_meta<- meta::metagen(TE = b, seTE = se, data = sig_merge, studlab = study, sm = "OR",
                            common = T,random = T)
  
  settings.meta("JAMA")
  meta::forest(sig_meta)
  out_plots(filename=paste0('Meta_BC_',i) ,pic_width=8,pic_height=6)

  inf <- meta::metainf(sig_meta)
  forest(inf)
  out_plots(filename=paste0('Inf_BC_',i),pic_width=8,pic_height=6)
  
  if (sig_meta$pval.Q>0.05 & sig_meta$I2<0.5) {
    if (sig_meta$pval.common<0.05) {
      same_exp_metasig <- c(same_exp_metasig,i)
    }
  }
  else{
    if (sig_meta$pval.random<0.05) {
      same_exp_metasig <- c(same_exp_metasig,i)
    }
  }
}
setwd("M:/pQTL_BC_MR/")
save(same_exp,same_exp_metasig,file = '3same_exp_metasig.rdata')

