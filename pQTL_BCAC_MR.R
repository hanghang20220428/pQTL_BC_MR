

rm(list = ls())
library(TwoSampleMR)
library(tidyverse)
library(data.table)
load("M:/pQTL_BC_MR/0BC_data_BCAC/BCAC_sus_2020_onco.rdata")

AllOR <- data.frame()

setwd("M:/pQTL_BC_MR/deCODE_as_expo_0.1_100kb/")
expfile <- list.files(pattern = 'rdata')
for (i in expfile) {
  setwd("M:/pQTL_BC_MR/deCODE_as_expo_0.1_100kb/")

  trait <-sapply(strsplit(i,split = '.rda'), '[',1)  
    
  load(i)
  
  out <- BCAC_sus_2020_onco %>% dplyr::filter(SNP %in% exp$SNP)
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
      res$exposure <- trait
      res$outcome <- "BCAC_BC"
      OR <- generate_odds_ratios(res)
      if (res$pval>0.05 | res$pval=='NaN') {
        OR <- cbind(OR,het_Qval=NA,pleio_pval=NA,
                    direction=NA)
        AllOR <- rbind(AllOR,OR)
        next }
      else{
        setwd("M:/pQTL_BC_MR/")
        if (OR$method == 'Inverse variance weighted') {
          het <- mr_heterogeneity(dat,method_list = 'mr_ivw')
          pleio <- mr_pleiotropy_test(dat)
          direct <- directionality_test(dat)
          OR <- cbind(OR,het_Qval=het$Q_pval,pleio_pval=pleio$pval,
                      direction=direct$correct_causal_direction)
          save(exp,out,dat,res,OR,file = paste0('./pQTL_BCAC_MR/',i,
                                                "_mrResult.rdata"))
          AllOR <- rbind(AllOR,OR)

          single <- mr_leaveoneout(dat)
          mr_leaveoneout_plot(single)
          ggsave(filename = paste0('./pQTL_BCAC_MR/',trait,'_leaveoneout_plot.pdf'),width = 6,height = 5)
          
          mr_scatter_plot(res,dat)
          ggsave(filename = paste0('./pQTL_BCAC_MR/',trait,'_mr_scatter_plot.pdf'),width = 6,height = 5)
          
          res_single <- mr_singlesnp(dat,all_method = c("mr_ivw")) 
          mr_forest_plot(res_single)
          ggsave(filename = paste0('./pQTL_BCAC_MR/',trait,'_mr_forest_plot.pdf'),width = 5,height = 6)
          
          mr_funnel_plot(res_single)
          ggsave(filename = paste0('./pQTL_BCAC_MR/',trait,'_mr_funnel_plot2.pdf'),width = 6,height = 4.5)
          
        }
        else{
          direct <- directionality_test(dat)
          OR <- cbind(OR,het_Qval=NA,pleio_pval=NA,
                      direction=direct$correct_causal_direction)
          save(exp,out,dat,res,OR,file = paste0('./pQTL_BCAC_MR/',trait,
                                                "_mrResult.rdata"))
          AllOR <- rbind(AllOR,OR)
        }
      }
    }
  }
}

setwd("M:/pQTL_BC_MR/")

#FDR校正
AllOR$FDR <- p.adjust(AllOR$pval,method = "BH")

save(AllOR,file = 'all_pQTL_mrres_BCAC.rdata')
write.csv(AllOR,file = 'all_pQTL_mrres_BCAC.csv')

AllOR_NOhete1 <- subset(AllOR,FDR<0.05 & het_Qval>0.05 & pleio_pval>0.05 & direction==TRUE) 
AllOR_NOhete2 <- subset(AllOR,FDR<0.05 & is.na(het_Qval) & direction==TRUE)
AllOR_NOhete_FDR <- rbind(AllOR_NOhete1,AllOR_NOhete2)
save(AllOR_NOhete_FDR,file = 'all_pQTL_mrres_BCAC_NOhete_FDR.rdata')
write.csv(AllOR_NOhete_FDR,file = 'all_pQTL_mrres_BCAC_NOhete_FDR.csv')

AllOR$label = ifelse(AllOR$BFR < 0.05,AllOR$short_expo, NA)

ggplot(AllOR,aes(x = b, y = -log10(FDR))) +
  geom_point(aes(size = abs(b)), alpha = 0.6, color = "#5E4FA2") +
  geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 6.5))+
  labs(x = "Beta (effect size)", 
       y = parse(text = "-log[10]*(FDR)"),
       title = 'All MR results between proteins and breast cancer (BCAC study)')+
  ggrepel::geom_label_repel(aes(label = label),size = 2)
ggsave(filename = 'volcano_plot_BCAC.pdf',width = 6.5,height = 4.5)                















    