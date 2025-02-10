

rm(list = ls())
library(tidyverse)
library(data.table)
setwd("M:/pQTL_BC_MR")

load("M:/pQTL_BC_MR/5.1coloc_DeCODE_Finn.rdata")
DeCODE_Finn <- all_coloc_results
colnames(DeCODE_Finn) <- 'DeCODE_FinnGen'

load("M:/pQTL_BC_MR/5.2coloc_DeCODE_BCAC.rdata")
DeCODE_BCAC <- all_coloc_results
colnames(DeCODE_BCAC ) <- 'DeCODE_BCAC'

load("M:/pQTL_BC_MR/5.3coloc_UKPPP_Finn.rdata")
UKPPP_Finn <- all_coloc_results
colnames(UKPPP_Finn) <- 'UKPPP_FinnGen'

load("M:/pQTL_BC_MR/5.4coloc_UKPPP_BCAC.rdata")
UKPPP_BCAC <- all_coloc_results
colnames(UKPPP_BCAC) <- 'UKPPP_BCAC'

data <- cbind(DeCODE_Finn,DeCODE_BCAC,UKPPP_Finn,UKPPP_BCAC)
save(data,file = '5.5all_coloc_results.rdata')
write.csv(data,file = '5.5all_coloc_results.csv')

#宽转长
data1 <- data %>% rownames_to_column(var = "Protein") %>% pivot_longer(!Protein,
                      names_to ='Group',values_to = 'PP.H4' )

data1 %>% ggplot(aes(x = Protein, y = Group, size=PP.H4)) +
  geom_point(color="#1F78B4") + 
  labs(x = NULL,y=NULL) +
  ggtitle("All PP.H4 results of colocalization") +
  theme_bw() +
  theme(panel.grid =element_blank(),
        panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        axis.text = element_text(size = 13),
        title = element_text(size = 16))
ggsave(filename = '5.5Dotplot_coloc.pdf',width = 6,height = 3.2)  



