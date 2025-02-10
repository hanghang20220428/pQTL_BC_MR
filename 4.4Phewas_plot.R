

rm(list = ls())
library(tidyverse)
library(data.table)
library(tibble)
library(circlize)
library(ComplexHeatmap)
setwd("M:/pQTL_BC_MR/")

load("M:/pQTL_BC_MR/4.3ALPI_phewas_inter.rdata")
finninfo <- fread('Finn_R9_data.csv',data.table = F)
head(ALPI_phewas_inter)
head(finninfo)
data0 <- merge(ALPI_phewas_inter,finninfo,by.x = 'outcome_DeCODE',by.y = 'phenocode')
head(data0)

data1 <- data0 %>% dplyr::select(name,exposure_DeCODE,b_DeCODE,category) %>% 
  dplyr::rename(exposure=exposure_DeCODE,beta=b_DeCODE)
data1$category <- sapply(strsplit(data1$category,split = '[(]'),'[',1 )
data1$category <- as.factor(data1$category)
data1$exposure <- paste0(data1$exposure,'_DeCODE')

data2 <- data0 %>% dplyr::select(name,exposure_UKPPP,b_UKPPP,category)%>% 
  dplyr::rename(exposure=exposure_UKPPP,beta=b_UKPPP)
data2$category <- sapply(strsplit(data2$category,split = '[(]'),'[',1 )
data2$category <- as.factor(data2$category)
data2$exposure<- paste0(data2$exposure,'_UKPPP')

data3 <- rbind(data1,data2)

data <- data3 %>% dplyr::arrange('category') %>% 
  pivot_wider(names_from=exposure,values_from = beta) %>% 
  column_to_rownames('name') %>% 
  dplyr::select(-category,category) %>% 
  mutate(dplyr::across(.cols=c('ALPI_DeCODE','ALPI_UKPPP'),.fns=as.character)) %>% 
  mutate(dplyr::across(.cols=c('ALPI_DeCODE','ALPI_UKPPP'),.fns=as.numeric))

levels(data$category)
col_n <- c("#a6cee3", "#1f78b4", "#b2df8a","#33a02c", "#fb9a99",
           "#e31a1c","#fdbf6f", "#ff7f00", "#cab2d6", "#62B197",
           "#E16E6D","#9392BE","#D0E7ED","#D5E4A8","#6a3d9a",
           "#4197d8","#E79397", "#f8c120","#CC88B0","#413496",
           '#083356','#ce1554',"#223e9c","#b12b23","#aebea6",
           "#edae11","#0f6657")
names(col_n) <- unique(data$category)
                    
col_fun = list(col_1 = colorRamp2(c(-1, 0, 1), 
                                  c("#0775b4", "white", "#ff7f00")),
               col_2 = colorRamp2(c(-1, 0, 1), 
                                  c("#6a3d9a", "white", "#e31a1c")),
               col_3 = col_n
)

#画图
if(T){
pdf("4.4Phewas_plot.pdf", height = 20, width = 20)
circos.par$gap.degree <- 60
circos.par$start.degree <- 30
circos.par$track.margin <- c(0.001, 0.001)

for (i in 1:3) {

  data_tmp <- as.matrix(data[,i])
  if (i == 1) {
    rownames(data_tmp) <- rownames(data)
  }
  colnames(data_tmp) <- colnames(data)[i]
  if (i<3) {
    circos.heatmap(data_tmp, 
                   col = col_fun[[i]],
                   rownames.side = "outside", 
                   rownames.cex = 0.8,cluster = F,
                   cell.border = "white",show.sector.labels = T,
                   track.height = 0.02)
  } else {
    circos.heatmap(data_tmp, 
                   col = col_fun[[i]],
                   rownames.side = "outside", 
                   cluster = T,
                   track.height = 0.02)
  }
}

lgd1 <- Legend(title = "DeCODE", border = "black", grid_height = unit(5, "mm"),
               legend_width = unit(20, "mm"),
               at = c(-1, 0, 1), title_position = "topcenter",
               col_fun = col_fun[[1]], direction = "horizontal")

lgd2 <- Legend(title = "UKPPP", border = "black", grid_height = unit(5, "mm"),
               legend_width = unit(20, "mm"),
               at = c(-1, 0, 1), title_position = "topcenter",
               col_fun = col_fun[[2]], direction = "horizontal")

pd <- packLegend(lgd1, lgd2, row_gap = unit(2, "mm"))
draw(pd, x = unit(0.6, "npc"), y = unit(0.7, "npc"))

lgd = Legend(labels = names(col_n),
             legend_gp = gpar(fill = col_n), title = "Category", 
             ncol = 1, row_gap = unit(1, "mm"))

draw(lgd, x = unit(0.515, "npc"), y = unit(0.5, "npc"))

dev.off()
circos.clear()
} 
                  
   
                    
