

rm(list = ls())
library(TwoSampleMR)
library(tidyverse)
library(data.table)
setwd("M:/pQTL_BC_MR/")
load('0BC_data_finngen_R9/finndataR9_BC_as_outcome.rdata')

setwd('F:/deCODE_pQTL/')
allp <- list.files(pattern = '*.gz')

AllOR <- data.frame()
Wrongfiles <- data.frame()

for (i in allp) {

  shorti <- sapply(strsplit(i,'.txt'),'[',1)

  setwd('F:/deCODE_pQTL/')


  data0 <- try(fread(i,data.table = F, integer64 = "numeric"))
  
  if('try-error' %in% class(data0)){
    Wrongfiles <- rbind(Wrongfiles,i)
    next
  }else{
    #整理数据
    data1 <- data0 %>% dplyr::select("rsids","effectAllele","otherAllele","Beta",
                                     "SE","ImpMAF","Pval","N","Chrom","Pos")
    colnames(data1) <- c("SNP","effect_allele","other_allele","beta","se",
                         "eaf",'p','samplesize',"chr","pos")
    data2 <- na.omit(data1)
    MHC <- data2 %>% filter(chr=='chr6' & pos>28000000 & pos<34000000)
    data3 <- data2[!data2$SNP %in% MHC$SNP,]
    
    exp_dat1 <- data3 %>% dplyr::select("SNP","p")
    
    colnames(exp_dat1)[1] <- 'rsid'
    colnames(exp_dat1)[2] <- 'pval'

    exp_dat2<- try(ieugwasr::ld_clump_local(exp_dat1,clump_p=5e-08,clump_r2=0.1,clump_kb=100,
                                        bfile="D:/1kg.v3/EUR",
                                        plink_bin="D:/Program files/plink_win64_20230116/plink.exe"))
    if('try-error' %in% class(exp_dat2)){
      next
     }else{
    if (nrow(exp_dat2)>0) {

      sameSNP <- intersect(data2$SNP,exp_dat2$rsid)
      data <- data2[data2$SNP %in% sameSNP,] %>% 
        dplyr::arrange(p) %>% 
        dplyr::distinct(SNP,.keep_all = T)
      
      exp <- TwoSampleMR::format_data(data,type = "exposure",  snp_col = "SNP",
                                          beta_col = "beta",se_col = "se",
                                          eaf_col = "eaf",effect_allele_col = "effect_allele",
                                          other_allele_col = "other_allele",
                                          pval_col = "p",samplesize_col = "samplesize",
                                          chr_col = "chr",pos_col = "pos")

      exp$R2<-exp$beta.exposure*exp$beta.exposure*2*(exp$eaf.exposure)*(1-exp$eaf.exposure)
      exp$f<-(exp$samplesize.exposure-2)*exp$R2/(1-exp$R2)
      

      exp <- exp[exp$f>10,]
      
      if (is.null(exp)) {
        next            }
      else{
        setwd('F:/deCODE_as_expo_0.1_100kb/')
        save(exp,file = paste0(shorti,'.rdata'))
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
            res$exposure <- shorti
            res$outcome <- "Finn_BC"
            OR <- generate_odds_ratios(res)
            if (res$pval>0.05 | res$pval=='NaN') {
              OR <- cbind(OR,het_Qval=NA,pleio_pval=NA,
                          direction=NA)
              AllOR <- rbind(AllOR,OR)
              next           }
            else{
              
              setwd("M:/pQTL_BC_MR/")
              if (OR$method == 'Inverse variance weighted') {
                het <- mr_heterogeneity(dat,method_list = 'mr_ivw')
                pleio <- mr_pleiotropy_test(dat)
                direct <- directionality_test(dat)
                OR <- cbind(OR,het_Qval=het$Q_pval,pleio_pval=pleio$pval,
                            direction=direct$correct_causal_direction)
                save(exp,out,dat,res,OR,file = paste0('./pQTL_FinnR9_MR/',shorti,
                                                      "_mrResult.rdata"))
                AllOR <- rbind(AllOR,OR)
                single <- mr_leaveoneout(dat)
                mr_leaveoneout_plot(single)
                ggsave(filename = paste0('./pQTL_FinnR9_MR/',shorti,'_leaveoneout_plot.pdf'),width = 6,height = 5)
                
                mr_scatter_plot(res,dat)
                ggsave(filename = paste0('./pQTL_FinnR9_MR/',shorti,'_mr_scatter_plot.pdf'),width = 6,height = 5)
                
                res_single <- mr_singlesnp(dat,all_method = c("mr_ivw")) 
                mr_forest_plot(res_single)
                ggsave(filename = paste0('./pQTL_FinnR9_MR/',shorti,'_mr_forest_plot.pdf'),width = 5,height = 6)
                
                mr_funnel_plot(res_single)
                ggsave(filename = paste0('./pQTL_FinnR9_MR/',shorti,'_mr_funnel_plot2.pdf'),width = 6,height = 4.5)
                
              }
              else{
                direct <- directionality_test(dat)
                OR <- cbind(OR,het_Qval=NA,pleio_pval=NA,
                            direction=direct$correct_causal_direction)
                save(exp,out,dat,res,OR,file = paste0('./pQTL_FinnR9_MR/',shorti,
                                                      "_mrResult.rdata"))
                AllOR <- rbind(AllOR,OR)
              }
            }
          }
        }
      }
    }
  }
  }
}

setwd("M:/pQTL_BC_MR/")

AllOR$FDR <- p.adjust(AllOR$pval,method = "BH")

save(AllOR,file = 'all_pQTL_mrres_Finn.rdata')
write.csv(AllOR,file = 'all_pQTL_mrres_Finn.csv')


AllOR_NOhete1 <- subset(AllOR,FDR<0.05 & het_Qval>0.05 & pleio_pval>0.05 & direction==TRUE) 
AllOR_NOhete2 <- subset(AllOR,FDR<0.05 & is.na(het_Qval) & direction==TRUE)
AllOR_NOhete_FDR <- rbind(AllOR_NOhete1,AllOR_NOhete2)
save(AllOR_NOhete_FDR,file = 'all_pQTL_mrres_Finn_NOhete_FDR.rdata')
write.csv(AllOR_NOhete_FDR,file = 'all_pQTL_mrres_Finn_NOhete_FDR.csv')

colnames(AllOR)
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
       title = 'All MR results between proteins and breast cancer (FinnGen study)')+
  ggrepel::geom_label_repel(aes(label = label),size = 2)
ggsave(filename = 'volcano_plot_Finn.pdf',width = 6.5,height = 4.5)                


    