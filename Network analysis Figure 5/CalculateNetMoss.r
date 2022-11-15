load("NetMoss_sourceData.RData")
# abundance tables were organized separately for each disease state
# net edges were organized into correlation matrix

load("microbeAbb.RData")

library(NetMoss2)
library(Hmisc)
library(dplyr)

# NetMoss ====================================

nodes_result_healthy.vs.preCOPD = 
  NetMoss(case_dir = abund_preCOPD,
          control_dir = abund_Health,
          net_case_dir = net_preCOPD,
          net_control_dir = net_Health)



nodes_result_preCOPD.vs.COPD = 
  NetMoss(case_dir = abund_COPD,
          control_dir = abund_preCOPD,
          net_case_dir = net_COPD,
          net_control_dir = net_preCOPD)


NMSS_1.2 = nodes_result_healthy.vs.preCOPD[[1]] %>% 
  mutate(taxon=sapply(taxon_names, function(x) microbeAbb$V2[which(microbeAbb$abb == x)])) %>%
  mutate(`1-2` = NetMoss_Score)


NMSS_2.3 = nodes_result_preCOPD.vs.COPD[[1]]  %>% 
  mutate(taxon=sapply(taxon_names, function(x) microbeAbb$V2[which(microbeAbb$abb == x)])) %>%
  mutate(`2-3` = NetMoss_Score)


res <- merge(NMSS_1.2, 
             NMSS_2.3, 
             by=c("taxon_names","taxon")) %>%
  select(taxon_names, taxon, `1-2`, `2-3`)
colnames(res)[1] <- "microbeAbb"

write.table(res, file = "NetMoss.txt", sep = "\t", quote = F, row.names = F)


