load("NetMoss_sourceData.RData")

library(dplyr)
library(Hmisc)
library(NetMoss2)

all(abund_COPD$X == abund_Health$X)
all(abund_Health$X == abund_preCOPD$X)

abund_all <- cbind.data.frame(abund_COPD, abund_Health %>% select(-X), abund_preCOPD %>% select(-X))

# calculate noise NetMoss score by randomly shuffle sample IDs 100 times -------------------------
AllShuff_Res <- NULL
for(shuf in 1:100){
  
  # shuf = 1
  
  shuffled  <- abund_all %>% tibble::column_to_rownames("X")
  
  set.seed(shuf)
  shuffled <- shuffled[,sample(1:699)]
  
  abund_COPD <- shuffled[,1:63]
  abund_Health <- shuffled[,64:559]
  abund_preCOPD <- shuffled[,560:699]
  
  # remove all 0 features
  abund_COPD <- abund_COPD[!rownames(abund_COPD) %in% names(which(apply(abund_preCOPD, 1, function(x) sum(x != 0)) <=1)), ]
  abund_Health <- abund_Health[!rownames(abund_Health) %in% names(which(apply(abund_Health, 1,function(x) sum(x != 0)) <=1)),]
  abund_preCOPD <- abund_preCOPD[!rownames(abund_preCOPD) %in% names(which(apply(abund_preCOPD, 1, function(x) sum(x != 0)) <=1)),]
  
  # COPD
  Corr <- rcorr(abund_COPD %>% t(), type="spearman")
  net_COPD <- Corr$r
  p_COPD <- Corr$P
  
  net_COPD_revised <- net_COPD
  net_COPD_revised[which(p_COPD > 0.05)] = 0
  net_COPD_revised[which(abs(net_COPD_revised) < 0.4)] = 0
  
  # preCOPD 
  Corr <- rcorr(abund_preCOPD %>% t(), type="spearman")
  net_preCOPD <- Corr$r
  p_preCOPD <- Corr$P
  
  net_preCOPD_revised <- net_preCOPD
  net_preCOPD_revised[which(p_preCOPD > 0.05)] = 0
  net_preCOPD_revised[which(abs(net_preCOPD_revised) < 0.4)] = 0
  
  # Health
  Corr <- rcorr(abund_Health %>% t(), type="spearman")
  net_Health <- Corr$r
  p_Health <- Corr$P
  
  net_Health_revised <- net_Health
  net_Health_revised[which(p_Health > 0.05)] = 0
  net_Health_revised[which(abs(net_Health_revised) < 0.4)] = 0
  
  # netMoss ----
  # Health vs preCOPD
  netMossRes_Heal.vs.pre = try(NetMoss(case_dir = abund_preCOPD,
                                           control_dir = abund_Health ,
                                           net_case_dir = net_preCOPD,
                                           net_control_dir = net_Health))
  
  netMossRes_Heal.vs.pre_revisedNet = try(NetMoss(case_dir = abund_preCOPD,
                                           control_dir = abund_Health ,
                                           net_case_dir = net_preCOPD_revised,
                                           net_control_dir = net_Health_revised))
  
  # Health vs COPD
  netMossRes_Heal.vs.copd = try(NetMoss(case_dir = abund_COPD,
                                           control_dir = abund_Health ,
                                           net_case_dir = net_COPD,
                                           net_control_dir = net_Health))
  
  netMossRes_Heal.vs.copd_revisedNet = try(NetMoss(case_dir = abund_COPD,
                                           control_dir = abund_Health ,
                                           net_case_dir = net_COPD_revised,
                                           net_control_dir = net_Health_revised))
  
  # preCOPD vs COPD
  netMossRes_pre.vs.copd = try(NetMoss(case_dir = abund_COPD,
                                           control_dir = abund_preCOPD,
                                           net_case_dir = net_COPD,
                                           net_control_dir = net_preCOPD))
  
  netMossRes_pre.vs.copd_revisedNet = try(NetMoss(case_dir = abund_COPD,
                                           control_dir = abund_preCOPD ,
                                           net_case_dir = net_COPD_revised,
                                           net_control_dir = net_preCOPD_revised))
  
  # summarise results
  if('try-error' %in% class(netMossRes_Heal.vs.pre) ){
    NNMS.1_2 <- cbind.data.frame(taxon_names = abund_all$X,
                                 NMSS.1_2 = NA, 
                                 stringsAsFactors=F)
  }else{
    NNMS.1_2 <- netMossRes_Heal.vs.pre[[1]] %>% mutate(NMSS.1_2 = NetMoss_Score) %>% select(taxon_names, NMSS.1_2)
  }
  
  if('try-error' %in% class(netMossRes_Heal.vs.pre_revisedNet) ){
    NNMS.1_2.revisedNet <- cbind.data.frame(taxon_names = abund_all$X,
                                            NMSS.1_2.revised = NA, 
                                            stringsAsFactors=F)
  }else{
    NNMS.1_2.revisedNet <- netMossRes_Heal.vs.pre_revisedNet[[1]] %>% 
      mutate(NMSS.1_2.revised = NetMoss_Score) %>% 
      select(taxon_names, NMSS.1_2.revised)
  }
  
  if('try-error' %in% class(netMossRes_Heal.vs.copd) ){
    NNMS.1_3 <- cbind.data.frame(taxon_names = abund_all$X,
                                 NMSS.1_3 = NA, 
                                 stringsAsFactors=F)
  }else{
    NNMS.1_3 <-  netMossRes_Heal.vs.copd[[1]] %>% mutate(NMSS.1_3 = NetMoss_Score) %>% select(taxon_names, NMSS.1_3)
  }
  
  if('try-error' %in% class(netMossRes_Heal.vs.copd_revisedNet) ){
    NNMS.1_3.revisedNet <- cbind.data.frame(taxon_names = abund_all$X,
                                            NMSS.1_3.revised = NA, 
                                            stringsAsFactors=F)
  }else{
    NNMS.1_3.revisedNet <-  netMossRes_Heal.vs.copd_revisedNet[[1]] %>% mutate(NMSS.1_3.revised = NetMoss_Score) %>% select(taxon_names, NMSS.1_3.revised)

  }
  
 
  if('try-error' %in% class(netMossRes_pre.vs.copd) ){
    NNMS.2_3 <- cbind.data.frame(taxon_names = abund_all$X,
                                 NMSS.2_3 = NA, 
                                 stringsAsFactors=F)
  }else{
    NNMS.2_3 <- netMossRes_pre.vs.copd[[1]] %>% mutate(NMSS.2_3 = NetMoss_Score) %>% select(taxon_names, NMSS.2_3)
    
  }
  

  if('try-error' %in% class(netMossRes_pre.vs.copd_revisedNet) ){
    NNMS.2_3.revisedNet <- cbind.data.frame(taxon_names = abund_all$X,
                                            NMSS.2_3.revised = NA, 
                                            stringsAsFactors=F)
  }else{
    NNMS.2_3.revisedNet <- netMossRes_pre.vs.copd_revisedNet[[1]] %>% mutate(NMSS.2_3.revised = NetMoss_Score) %>% select(taxon_names, NMSS.2_3.revised)
    
  }
 
  
  
  shuffRes <- merge(merge(merge(NNMS.1_2, NNMS.1_2.revisedNet, by = "taxon_names", all=T),
              merge(NNMS.1_3, NNMS.1_3.revisedNet, by="taxon_names", all=T),
              by="taxon_names", all=T),
        merge(NNMS.2_3, NNMS.2_3.revisedNet, by="taxon_names", all=T),
        by="taxon_names", all=T) %>%
    mutate(shuffle = shuf)
  
  AllShuff_Res <- bind_rows(AllShuff_Res, shuffRes)
}


# one sample t test to calculate whether the NetMoss score for each node are larger than noise NetMoss score  -------------------------
library(dplyr)

dat.shuf <- AllShuff_Res


dat.1sp <- data.table::fread("NetMoss.txt", data.table = F)

results <- NULL
for(nd in dat.1sp$microbeAbb){
  # nd=dat.1sp$microbeAbb[1]
  
  # 1_2 ============================
  value.1_2 = dat.1sp$`1-2`[which(dat.1sp$microbeAbb == nd)]
  shuff.1_2 = dat.shuf$NMSS.1_2.revised[which(dat.shuf$taxon_names == nd)]
  wilcoxP.1_2 <- wilcox.test(shuff.1_2, mu = value.1_2, alternative = "less")$p.value
  ttestP.1_2 <- t.test(shuff.1_2, mu = value.1_2, alternative = "less")$p.value
  
  
  # 2-3 ============================
  value.2_3 = dat.1sp$`2-3`[which(dat.1sp$microbeAbb == nd)]
  shuff.2_3 = dat.shuf$NMSS.2_3.revised[which(dat.shuf$taxon_names == nd)]
  wilcoxP.2_3 <- wilcox.test(shuff.2_3, mu = value.2_3, alternative = "less")$p.value
  ttestP.2_3 <- t.test(shuff.2_3, mu = value.2_3, alternative = "less")$p.value
  
  res_c <- c("node"=nd,
             "ttestPval.1-2"=ttestP.1_2, "wilcoxPval.1-2"=wilcoxP.1_2, 
             "ttestPval.2-3"=ttestP.2_3, "wilcoxPval.2-3" = wilcoxP.2_3)
  results <- bind_rows(results, res_c)
  
}

results$taxon = sapply(results$node, function(x) dat.1sp$taxon[which(dat.1sp$microbeAbb == x)])
results$`ttestPval.1-2` <- as.numeric(results$`ttestPval.1-2`)
results$`wilcoxPval.1-2` <- as.numeric(results$`wilcoxPval.1-2`)
results$`ttestPval.2-3` <- as.numeric(results$`ttestPval.2-3`)
results$`wilcoxPval.2-3` <- as.numeric(results$`wilcoxPval.2-3`)

write.table(results, file = 'NetMoss_oneSampleTest.w.100shuffle.txt', quote = F, sep = "\t", row.names = F)
