# calculate and plot the networks in three datasets

# ===============================================================================
# Figure 5a    
# ===============================================================================


# calculate networks ===================
library(dplyr)
library(data.table)

meta <- fread("../_data/meta_COPD.txt")
meta$Disease[meta$Disease == "Pre-COPD"] <- "preCOPD"
unique(meta$Disease)

load("../_data/microbialFeatures_relAbund_rpkm_notStandardized.RData")
load("microbeAbb.RData")

colnames(bacteria_rel) <- 
  sapply(colnames(bacteria_rel),
         function(x) if(x %in% microbeAbb$V2) microbeAbb$abb[which(microbeAbb$V2 == x)] else x)

colnames(fungi_rel) <- 
  sapply(colnames(fungi_rel),
         function(x) if(x %in% microbeAbb$V2) microbeAbb$abb[which(microbeAbb$V2 == x)] else x)

combined.dat <- 
  merge(bacteria_rel %>% select(-bact.shannon),
        fungi_rel %>% select(-fung.shannon), 
        by = 0) %>%
  filter(Row.names %in% meta$`#NAME`) %>%
  mutate(Group = sapply(Row.names, function(x) meta$Disease[which(meta$`#NAME` == x)]) )


library(Hmisc)


edges_3Diseases <- NULL
for(dss in unique(meta$Disease)){
  
  dat.tmp <- combined.dat %>% 
    filter(Group == dss) %>% 
    tibble::column_to_rownames("Row.names") %>%
    select(-Group)
  
  Corr <- rcorr(as.matrix(dat.tmp) , type="spearman")
  occor.r <- Corr$r
  occor.p <- Corr$P
  
  # Hide lower triangle for r (so each correlation is only going to appear once)
  lower <- occor.r
  lower[lower.tri(occor.r)]<-NA
  
  # edges
  r_df.l <- lower %>% reshape2::melt(value.name = "r")
  p_df.l <- occor.p %>% reshape2::melt(value.name = "p")
  
  if(all(r_df.l$Var1 == p_df.l$Var1) & all(r_df.l$Var2 == p_df.l$Var2)){
    occor.rp_l <- cbind.data.frame(r_df.l, p=p_df.l$p, stringsAsFactors=F) # quicker
  }else{
    occor.rp_l <- merge(occor.r %>% reshape2::melt(value.name = "r"), 
                        occor.p %>% reshape2::melt(value.name = "p"),
                        by=c("Var1", "Var2"))  # too slowwwwwwwww
  }
  
  
  edges <- occor.rp_l %>% 
    filter(!is.na(r)) %>% 
    filter(abs(r) < 1) %>% 
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) 
  
  edges$edge = paste(edges$Var1, edges$Var2, sep = "|")
  edges$disease = dss
  edges$edgeType <-
    paste0(substr(edges$edge,1,1), 
           substr(sapply(strsplit(edges$edge, "|", fixed = T), "[[", 2), 1, 1))
  
  edges_3Diseases <- bind_rows(edges_3Diseases, edges)
  
}

edges_3Diseases$r <- as.numeric(edges_3Diseases$r)
edges_3Diseases$p <- as.numeric(edges_3Diseases$p)

edges_3Diseases$source <- sapply(edges_3Diseases$Var1,
                                 function(x) microbeAbb$V2[which(microbeAbb$abb == x)])
edges_3Diseases$target <- sapply(edges_3Diseases$Var2,
                                 function(x) microbeAbb$V2[which(microbeAbb$abb == x)])

for(dss in c( "COPD","Health","preCOPD")){
  
  edges <- edges_3Diseases %>% 
    filter(disease == dss) %>%
    mutate(abs.r =  abs(r)) %>%
    filter(abs.r > 0.4) %>%
    filter(p < 0.05) %>%
    relocate(source, target) 
  
  write.csv(edges, file = paste0("edges/",dss,"_edges.csv"),quote = F, row.names = F)
}

# for each edge file, calculate net modules using Gephi
# Manually identify overlapped nodes in modules between healthy and preCOPD nets, and between preCOPD and COPD nets
# 
# nodes and corresponding module colors were saved in : Module_colors.RData

remove(list = ls())
load("NodeModuleColors.RData")
load("microbeAbb.RData")

dss = "COPD" # manually change: Health, preCOPD, COPD
edges <- fread(paste0("edges/",dss,"_edges.csv"), data.table = F) 
nodeColors <- eval(parse(text = paste0("nodesColors_",dss)))

extraNodes <- microbeAbb$V2[!microbeAbb$V2 %in% as.character(nodeColors$Id )]
selfEdges <- cbind.data.frame(source = extraNodes, 
                               target = extraNodes,
                               edgeType = "selfLink",
                               stringsAsFactors=F) 

edges <- bind_rows(edges, selfEdges)

nodes <- 
  data.frame(table(c(edges$source, edges$target)) ) %>%
  mutate(id = Var1) %>% select(-Var1) %>% 
  mutate(abb = sapply(id, function(x) microbeAbb$abb[which(microbeAbb$V2 == x)])) %>%
  mutate(type = sub("_\\d+","", abb)) %>%
  mutate(phylum =  sapply(id, function(x) microbeAbb$phylum[which(microbeAbb$V2 == x)])) %>%
  mutate(phylum6 = sapply(phylum, 
                          function(x){
                            if(x %in% c("Bacteroidetes","Ascomycota","Basidiomycota","Proteobacteria","Firmicutes")) x else "Others"
                          } )) %>%
  relocate(id)  



# plot with igraph ====================
library(igraph)
net <- graph_from_data_frame(d=edges, vertices=nodes, directed=F) 
net


V(net)$module.color <- 
  sapply(names(V(net)),
         function(x){
           if(x %in% nodeColors$Id){
             nodeColors$color[which(nodeColors$Id == x)]
           }else "gray"
         }) 

df <- cbind.data.frame(
  V(net)$module.color,
  V(net)$phylum6,
  V(net)$type
)
assign(paste0(dss,"_df"), df, envir = .GlobalEnv)

V(net)$type.color <-
  sapply(1:length(V(net)),
         function(i){
           if(V(net)$type[i] == "bacteria") "white" else V(net)$module.color[i]
         })


E(net)$source.module.color  <- 
  sapply(E(net)$Var1,
         function(x) {
           if(is.na(x) | !x %in% nodeColors$Id) "gray" else nodeColors$color[which(nodeColors$Id == x)]
         })


net.simp <- net - E(net)[E(net)$edgeType=="selfLink"]  #remove self links, only keep the nodes


# Fig5a:
pdf(paste("NetPlots/", dss,".pdf",sep = ""), 
    width = 4 , height = 4) 
#储存的图片大小（对应圆的面积）跟node个数成正比，所以长宽与node个数的平方根成正比
par(mar=c(0.1,0.1,0.1,0.1)) 
set.seed(100)

plot(net.simp, vertex.label=NA, 
     vertex.size = 4,
     vertex.color = V(net)$type.color,
     vertex.frame.color= V(net)$module.color,
     edge.size = 1, 
     edge.color = E(net)$source.module.color, 
     layout=layout_with_fr(net)
)

dev.off()


# ===============================================================================
# Figure 5b    
# ===============================================================================


tmp <- 
  rbind.data.frame(
    Health_df %>% mutate(Disease = "Health"),
    preCOPD_df %>% mutate(Disease = "preCOPD"),
    COPD_df %>% mutate(Disease = "COPD"),
    stringsAsFactors = F
  ) %>% 
  mutate(module.color=`V(net)$module.color`,
         phylum = `V(net)$phylum6`,
         type = `V(net)$type`) %>%
  mutate(module =  sapply(module.color,
                          function(x){
                            if(x %in% c("#43997A", "#47A265")){
                              "green"
                            }else if(x %in% c("#7F5477", "#91569B")){
                              "purple"
                            }else if(x == "#5D6795"){
                              "blue"
                            }else if(x == "#CB6651"){
                              "brown"
                            }else if(x == "#E41A1C"){
                              "red"
                            }else if(x == "#FFD422"){
                              "yellow"
                            }else x
                          }))
plotDat.phylumPie <- tmp %>%
  group_by(module, phylum) %>%
  summarise(n=n()) %>%
  #calculate relative proportions
  group_by(module) %>%
  mutate(perc=n/sum(n))

library(ggplot2)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

Fig5b.phylumPie <- ggplot(plotDat.phylumPie %>% filter(module != "gray"), aes(x="", y=perc, fill=phylum))+
  facet_grid(.~module) + 
  geom_bar(width = 1, stat = "identity", color="white") +
  coord_polar("y", start=0) +
  blank_theme +
  theme(axis.text.x=element_blank())+
  scale_fill_manual(values = c('#84AADA','#65A484',"#D2BDDB","#E1E0AB","#DAD7D7","#D56565"))



# plot module size ---------------------
plotDat.moduleSize <- tmp %>%
  filter(module != "gray") %>%
  group_by(Disease, module) %>%
  summarise(n=n())


plotDat.moduleSize$Disease <- factor(plotDat.moduleSize$Disease, levels = c("Health","preCOPD","COPD"))
Fig5b.moduleSize <- ggplot(plotDat.moduleSize ) +
  geom_point(aes(x=Disease, y=module, size=n, color=module), shape = 21) +
  scale_color_manual(values = c("#5D6795","#CB6651","#43997A","pink2","#7F5477","#E41A1C","#FFD422")) +
  scale_size(range = c(2, 12)) +
  theme_bw() + theme(panel.grid = element_blank())
Fig5b.moduleSize


# plot module type percentage (bacteria vs fungi ) --------------
plotDat.moduleTypePerc <- tmp %>%
  filter(module != "gray") %>%
  group_by(Disease, module, type) %>%
  summarise(n=n()) %>%
  group_by(Disease, module) %>%
  mutate(perc = n/sum(n)) 

plotDat.moduleTypePerc$fillType <- 
  sapply(1:nrow(plotDat.moduleTypePerc),
         function(i){
           if(plotDat.moduleTypePerc$type[i] == "bacteria"){
             "white"
           }else{
             if(plotDat.moduleTypePerc$module[i] == "blue"){
               "#5D6795"
             }else if(plotDat.moduleTypePerc$module[i] == "brown"){
               "#CB6651"
             }else if(plotDat.moduleTypePerc$module[i] == "green"){
               "#43997A"
             }else if(plotDat.moduleTypePerc$module[i] == "pink2"){
               "pink2"
             }else if(plotDat.moduleTypePerc$module[i] == "purple"){
               "#7F5477"
             }else if(plotDat.moduleTypePerc$module[i] == "red"){
               "#E41A1C"
             }else if(plotDat.moduleTypePerc$module[i] == "yellow"){
               "#FFD422"
             }
           }
         })

plotDat.moduleTypePerc$Disease <- factor(plotDat.moduleTypePerc$Disease, levels = c("Health","preCOPD","COPD"))
Fig5b.moduleTypePerc <- ggplot(plotDat.moduleTypePerc, aes(x="", y=perc, fill=fillType)) + 
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  facet_grid(module ~ Disease) +
  scale_fill_manual(values =  unique(plotDat.moduleTypePerc$fillType)[order(unique(plotDat.moduleTypePerc$fillType))])+
  #scale_color_manual(values = unique(plotDat.moduleTypePerc$module)[order(unique(plotDat.moduleTypePerc$module))]) +
  theme_void() # remove background, grid, numeric labels


# manually merge Fig5b.moduleTypePerc and Fig5b.moduleSize



# ===============================================================================
# Figure 5c    
# ===============================================================================
remove(list = ls())
dat <- fread("NetMoss.txt", data.table = F)
head(dat)

# color by phylum and node significance
load("microbeAbb.RData")
dat$phylum = sapply(dat$taxon,
                    function(x){
                      microbeAbb$phylum[which(microbeAbb$V2 == x)]
                    })

# see "shuffleSample_NetMoss_noise.r" for how nodes significance were calcualted 
noise <- fread("NetMoss_oneSampleTest.w.100shuffle.txt",data.table = F)
insig.nodes <- noise$node[which(noise$`ttestPval.1-2` > 0.05 & noise$`ttestPval.2-3` > 0.05) ]


dat$color = sapply(1:nrow(dat),
                   function(i){
                     if(dat$microbeAbb[i] %in% insig.nodes){
                       "gray"
                     }else if(dat$phylum[i] == "Bacteroidetes"){
                       "#E1E0AB"
                     }else if(dat$phylum[i] == "Firmicutes"){
                       "#65A484"
                     }else if(dat$phylum[i] == "Proteobacteria"){
                       "#D56565"
                     }else if(dat$phylum[i] == "Ascomycota"){
                       "#84AADA"
                     }else if(dat$phylum[i] == "Basidiomycota"){
                       "#D2BDDB"
                     }else{
                       "#C6BA9A"
                     }
                   })


# mean relative abundance in all samples as circle size
load("../_data/microbialFeatures_relAbund_rpkm_notStandardized.RData")

tmp1 <- c(sapply(bacteria_rel, mean),sapply(fungi_rel, mean))

dat$relAbund = sapply(dat$taxon,
                      function(x){
                        tmp1[which(names(tmp1) == x)]
                      })


colorlvls <- unique(dat$color)[order(unique(dat$color))]




Fig5c <- ggplot(dat) +
  geom_point(aes(x=`1-2`, y=`2-3`, color=color, size=relAbund)) +
  scale_size(range = c(3,4)) +
  scale_color_manual(values = colorlvls) +
  theme_bw() +theme(panel.grid = element_blank())

