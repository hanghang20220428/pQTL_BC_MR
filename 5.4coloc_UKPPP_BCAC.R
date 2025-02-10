


rm(list = ls())
library(tidyverse)
library(data.table)
library(coloc)
setwd("M:/pQTL_BC_MR")
load("M:/pQTL_BC_MR/3same_exp_metasig.rdata")
same_exp_metasig <- same_exp_metasig[2:4]
load('UKBPPP_SNP_ref.rdata')

all_coloc_results <- data.frame()
for (i in same_exp_metasig) {
  load("M:/0BC_data_BCAC/2BCAC_2020_suscep_as_coloc.rdata")

  setwd('F:/2940UK_pQTL data')
  allp <- list.files(pattern = '.rdata$')
  load(allp[grep(paste0("^",i,'_'),allp)])

  data$SNP <- UKBPPP_SNP$rsid[match(data$ID,UKBPPP_SNP$ID)]
  data$Pval <- 10^(-data$LOG10P)
  data$samplesize <- 54219
  
  head(data)
  data1 <- data %>% dplyr::select("SNP","CHROM","GENPOS","ALLELE1","ALLELE0",
                                  "A1FREQ","BETA","SE","Pval","samplesize")
  data2 <- na.omit(data1)
  colnames(data2) <- c('SNP','chrom',"Pos",'effect_allele','other_allele',
                       "eaf","beta","se","P","samplesize")
  data3 <- as.data.frame(data2)
  data3$varbeta <- data3$se^2
  data3$MAF <- data3$eaf
  data3$z = data3$beta/data3$se
  lead <- data3 %>% dplyr::arrange(P)
  leadchr <- lead$chrom[1]
  leadPos <- lead$Pos[1]
  QTLdata <- data3[data3$chrom==leadchr & data3$Pos>leadPos-1000000 & data3$Pos<leadPos+1000000,] %>% 
    na.omit() %>% dplyr::distinct(SNP,.keep_all = T)
  
  GWASdata <- GWASdata[GWASdata$SNP %in% QTLdata$SNP, ]
  
  
  coloc_data <- list(dataset1=list(snp=QTLdata$SNP,beta=QTLdata$beta,varbeta=QTLdata$varbeta,
                                   N=QTLdata$samplesize,MAF=QTLdata$MAF,z = QTLdata$z,
                                   pvalues=QTLdata$P,type="quant"), 
                     dataset2=list(snp=GWASdata$SNP,beta=GWASdata$beta,varbeta=GWASdata$varbeta,
                                   N=GWASdata$samplesize,MAF=GWASdata$MAF,z = GWASdata$z,
                                   pvalues=GWASdata$P,type="cc",s=GWASdata$s))
  
  result <- coloc.abf(dataset1=coloc_data$dataset1, dataset2=coloc_data$dataset2)
  
  result$results %>% dplyr::arrange(desc(SNP.PP.H4))
  SNPresult <- result$results %>% dplyr::arrange(desc(SNP.PP.H4))
  
  co_result <-data.frame(SNPresult$SNP.PP.H4[1]) 

  rownames(co_result) <- i 
  all_coloc_results <- rbind(all_coloc_results,co_result)
  
  setwd("M:/pQTL_BC_MR/7.4coloc_UKPPP_BCAC/")
  save(result,file =paste0(i,'_BCAC_coloc.rdata'))
  
  QTL_coloc<- QTLdata %>%  dplyr::select(SNP,P) %>% dplyr::rename(rsid=SNP,pval=P) 
  GWAS_coloc <- GWASdata %>%  dplyr::select(SNP,P)%>% dplyr::rename(rsid=SNP,pval=P)
  fwrite(QTL_coloc,file = 'QTL_coloc.tsv',sep = "\t",row.names = F,col.names = T)
  fwrite(GWAS_coloc,file = 'GWAS_coloc.tsv',sep ="\t",row.names = F,col.names = T )
  
  p <- locuscomparer::locuscompare(in_fn1 ='GWAS_coloc.tsv',  in_fn2 ='QTL_coloc.tsv',
                                   title1 = "Breast cancer risk (BCAC)",title2 = paste0(i,"(pQTL)"),
                                   genome = "hg38",legend_position = 'topleft')
  
  pdf(file =paste0(i,"_BCAC_coloc.pdf"),width = 10,height = 8)
  print(p)
  dev.off()
  
  file.remove('QTL_coloc.tsv')
  file.remove('GWAS_coloc.tsv')
}

setwd("M:/pQTL_BC_MR")
save(all_coloc_results,file = '5.4coloc_UKPPP_BCAC.rdata')
fwrite(all_coloc_results,file ='5.4coloc_UKPPP_BCAC.csv',row.names = T) 






