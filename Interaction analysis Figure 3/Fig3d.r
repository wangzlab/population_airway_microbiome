library(data.table)
library(dplyr)


load("../_data/microbialFeatures_relAbund_rpkm_notStandardized.RData")
grouping_df<- fread("occupation_PC1_meta.txt", data.table = F)
for (df.name in c("ARG_rel","ARG_rpkm","bacteria_rel","fungi_rel","KO_rel",
                  "MetagMod","MetagTaxa_rel","VF_rel","VF_rpkm" )) {
  
  # df.name = "MetagMod"
  
  mf_dat <- eval(parse(text = df.name))
  
  mf_dat.l <- mf_dat %>%
    tibble::rownames_to_column("Sample") %>%
    mutate(Occupation =  sapply(Sample,
                                function(x){
                                  if( x %in% grouping_df$`#NAME`){
                                    grouping_df$Occupation[which(grouping_df$`#NAME` == x)]
                                  }else NA
                                } ) ) %>%
    mutate(NeissGroup = sapply(Sample,
                               function(x){
                                 if( x %in% grouping_df$`#NAME`){
                                   grouping_df$Group[which(grouping_df$`#NAME` == x)]
                                 }else NA
                               }) ) %>%
    reshape2::melt(id.vars=c("Sample","Occupation","NeissGroup"), 
                   variable.name = "MicrobialFeature")
  
  FC <- mf_dat.l %>%
    group_by(MicrobialFeature, NeissGroup) %>%
    summarise(mean.Occup_Y = mean(value[which(Occupation == "Y")]),
              mean.Occup_N = mean(value[which(Occupation == "N")])) %>%
    filter(!is.na(NeissGroup))
  
  
  if(df.name == "MetagMod") {
    FC <- FC %>% 
      mutate(FC = mean.Occup_Y - mean.Occup_N) %>%
      mutate(logFC = FC)  # to keep consistent with other microbial features, use logFC as the colname, but only need the sign (+ or -)
  }else{
    FC <- FC %>% 
      mutate(FC = (mean.Occup_Y + 0.00001)/(mean.Occup_N + 0.00001)) %>%
      mutate(logFC = log2(FC))
  }
  
  assign(x=paste0("FC_", df.name), value = FC)
  write.table(FC, file = paste0("FC_output/FC_", df.name, ".txt"), sep = "\t", row.names = F, quote = F )
  
}

# plotting =================================================
library(data.table)
library(dplyr)

data.n <- "bacteria_rel" # manually change: ARG_rel, bacteria_rel, fungi_rel, MetagMod, VF_rel
Pval.f <- "microbial_occupation_association/lmRes_metagMod.features.st_occupation.txt" # mannual change

FC.f <- list.files(path = "FC_output/", pattern = data.n, full.names = T)
FC <- fread(FC.f,data.table = F) 

PVal <- fread(Pval.f, data.table = F) 

tmp <- merge(
  PVal %>% select(taxon, group, Pvalue),
  FC %>% select( MicrobialFeature, NeissGroup, logFC) ,
  by.x=c("taxon", "group"), by.y=c("MicrobialFeature","NeissGroup")
)

res <- tmp %>%
  mutate(log10P = -log10(Pvalue)) %>%
  mutate(value = sign(logFC) * log10P)

plotDat <- res %>%
  reshape2::dcast(taxon~group, value.var = "value")

plotDat[is.na(plotDat)] = 0

plotDat$Significance <- sapply(1:nrow(plotDat),
                               function(i){
                                 if(abs(plotDat$Neis_high[i]) < -log10(0.05)){
                                   if(abs(plotDat$Neis_low[i]) < -log10(0.05)){
                                     "insignificant"
                                   }else{
                                     "sig_Neis_Low"
                                   }
                                 }else{
                                   if(abs(plotDat$Neis_low[i]) < -log10(0.05)){
                                     "sig_Neis_high"
                                   }else{
                                     "sig_Neiss_both"
                                   }
                                 }
                               })


for(i in 1:nrow(plotDat)){
  # i=1
  
  if(abs(plotDat$Neis_high[i]) < -log10(0.05)){
    if(abs(plotDat$Neis_low[i]) < -log10(0.05)){
      plotDat$Significance[i] <- "insignificant"
    }else{
      plotDat$Significance[i] <- "sig_Neis_Low"
    }
  }else{
    if(abs(plotDat$Neis_low[i]) < -log10(0.05)){
      plotDat$Significance[i] <- "sig_Neis_high"
    }else{
      plotDat$Significance[i] <- "sig_Neiss_both"
    }
  }
}

relAbund_df <- eval(parse(text = data.n))
relAbund_vec <- sapply(relAbund_df, mean) 

plotDat$relAbund <- sapply(plotDat$taxon, function(x) relAbund_vec[names(relAbund_vec) == x])

# plotting ================================
library(ggplot2)
library(ggrepel)

p<-ggplot(plotDat, aes(x= Neis_high, y=Neis_low) ) +
  geom_point(aes(color=Significance, size=relAbund)) +
  scale_color_manual(values = c("gray","indianred","skyblue3","gold3"))+
  scale_size(range = c(1,3)) +
  # geom_text_repel(aes(label=taxon)) +
  theme_bw() + theme(panel.grid = element_blank()) +
  geom_vline(xintercept = c(-1.30103,1.30103), linetype="dashed", color="darkgray") +
  geom_hline(yintercept = c(-1.30103,1.30103), linetype="dashed", color="darkgray")

p

ggsave(p, filename = paste0("Fig3d_",data.n,".pdf"), device = "pdf", width = 5.5, height = 4)

