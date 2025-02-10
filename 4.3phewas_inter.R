



rm(list = ls())
library(tidyverse)
library(data.table)
setwd("M:/pQTL_BC_MR/")

load("M:/pQTL_BC_MR/4.1DeCODE_ALPI_phewas_MR_BFR_NOhete.rdata")
DeCODE_ALPI_phewas <- AllOR_NOhete
colnames(DeCODE_ALPI_phewas) <- paste0(colnames(DeCODE_ALPI_phewas),'_DeCODE')

load("M:/pQTL_BC_MR/4.2UKPPP_ALPI_phewas_MR_BFR_NOhete.rdata")
UKPPP_ALPI_phewas <- AllOR_NOhete
colnames(UKPPP_ALPI_phewas) <- paste0(colnames(UKPPP_ALPI_phewas),'_UKPPP')

same_outcome <- Reduce(intersect, list(DeCODE_ALPI_phewas$outcome_DeCODE, UKPPP_ALPI_phewas$outcome_UKPPP))

DeCODE_ALPI_phewas <- DeCODE_ALPI_phewas %>% dplyr::filter(outcome_DeCODE %in% same_outcome)
UKPPP_ALPI_phewas <- UKPPP_ALPI_phewas %>% dplyr::filter(outcome_UKPPP %in% same_outcome)

ALPI_phewas_inter <- cbind(DeCODE_ALPI_phewas,UKPPP_ALPI_phewas)
save(ALPI_phewas_inter,file = '4.3ALPI_phewas_inter.rdata')
fwrite(ALPI_phewas_inter,file ='4.3ALPI_phewas_inter.csv')




