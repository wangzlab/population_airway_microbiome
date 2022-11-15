library("migest")
library(data.table)
#library(dplyr)
library(tidyverse)


dat <- fread("circos_20220818.txt") 
colnames(dat) <- c("Exposure","MicroType","Microbiome","Pval","log10P")

# adjust the log10P value to a range 
library(scales)
dat$flow <-  rescale(dat$log10P, to = c(1.5, 5))


dat <- dat %>%
  mutate(Exposure_MicroType = paste(Exposure, MicroType,sep = "_"))

dat.d <- dat %>% 
  mutate(orig = Exposure_MicroType,
         #   flow=log10P,
         dest = Microbiome ) %>%
  select(orig, dest, flow, Pval)


# define the colors  ------

config_type <- cbind.data.frame(
  order = 1:10,
  names = c(unique(dat$Exposure), unique(dat$MicroType)),
  labels = c(unique(dat$Exposure), unique(dat$MicroType)),
  col = c("#FF0000", "#D2960C", "#7DAF00", "#7500FF", "#A0007D", # five colors for exposures
          "#95C6C6", "#49B1DD", "#7EABCA", "#E98DAF", "#E9CF91"),  # five colors for Microtype
  stringsAsFactors=F
)

config_element <- cbind.data.frame(
  order = 1:length(c(unique(dat$Exposure_MicroType), unique(dat$Microbiome))),
  names = c(unique(dat$Exposure_MicroType), unique(dat$Microbiome)),
  labels = c(rep(NA, length(unique(dat$Exposure_MicroType))), unique(dat$Microbiome)), 
  stringsAsFactors=F
) %>% 
  mutate(type = sapply(names, function(x){
    if(x %in% dat$Exposure_MicroType) {
      strsplit(x,"_",fixed = T)[[1]][1]
    } else {
      dat$MicroType[which(dat$Microbiome == x)[1]]
    }
  }))%>%
  mutate(typeLabel = sapply(type,
                            function(x) if(x %in% c("ko","arg","vf")) toupper(x) else x )) %>%
  mutate(col = sapply(type, function(x) config_type$col[which(config_type$names == x)])) %>%
  mutate(type_end = sapply(names, function(x){
    if(x %in% dat$Exposure_MicroType) {
      strsplit(x,"_",fixed = T)[[1]][2]
    } else {
      dat$MicroType[which(dat$Microbiome == x)[1]]
    }
  })) %>%
  mutate(col_end = sapply(type_end, function(x) config_type$col[which(config_type$names == x)]))


config_element$type <- factor(config_element$type, levels = c(config_type$labels))

config_element <- config_element %>%
  arrange(names) %>% arrange(type)

# find appropriate label position for bend2
unique(config_element$type)
i <- which(config_element$type == "Module")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""


i <- which(config_element$type == "Bacteria")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] <- ""


i <- which(config_element$type == "Fungi")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] <- ""


i <- which(config_element$type == "ARG")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] <- ""

i <- which(config_element$type == "VF")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] <- ""


i <- which(config_element$type == "Biofuel")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""

i <- which(config_element$type == "Occupation")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""

i <- which(config_element$type == "PM25")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""


i <- which(config_element$type == "SHS")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""


i <- which(config_element$type == "Smoking")
i.l <- round(median(i),0)
config_element$typeLabel[i[i != i.l] ] = ""



# label with ARG and VF names ---------------
ARG.abb <- fread("arg_name_mapping.txt", header = F, col.names = c("abb","ARG"), data.table = F)
VF.abb <- fread("VFID.names.txt", data.table = F, col.names = c("abb","VF"))

config_element$labels_full <-
  sapply(config_element$labels,
         function(x){
           if(is.na(x)){
             NA
           }else if(startsWith(x,"ARG")){
             ARG.abb$ARG[which(ARG.abb$abb == x)] 
           }else if(startsWith(x,"VF")){
             VF.abb$VF[which(VF.abb$abb == x)]
           }else x
           
         })


# plot with the names of element types 
pdf("1_4.circular_typeNames.pdf")
mig_chord(x = dat.d, 
          order = config_element$names,
          grid.col = config_element %>% 
            select(names, col_end) %>%
            deframe(), 
          #transparency = dat.d$Pval,
          
          # lab = config_element %>% select(names, labels_full) %>% deframe(), 
          lab_bend2 = config_element %>% select(names, typeLabel) %>% deframe(), 
          
          gap.degree = 1, 
          label_size = 1,
          no_axis = T
          #axis_breaks = 100
)
dev.off()

# plot with the names of elements 
pdf("1_4.circular_elementFullnames.pdf")
mig_chord(x = dat.d, 
          order = config_element$names,
          grid.col = config_element %>% 
            select(names, col_end) %>%
            deframe(), 
          #transparency = dat.d$Pval,
          
          lab = config_element %>% select(names, labels_full) %>% deframe(), 
          # lab_bend2 = config_element %>% select(names, typeLabel) %>% deframe(), 
          
          gap.degree = 1, 
          label_size = 0.5,
          no_axis = T
          #axis_breaks = 100
)
dev.off()



# plot the zscore from lm analysis as an out circle heatmap =============
# remove(list = ls())


load("plotSourceData/integrated_microbial.exposure.associations.RData") 
combined_dat$Zscore <- -qnorm(combined_dat$Pvalue/2)
combined_dat$Zscore <- combined_dat$Zscore * sign(combined_dat$effect.size)



all(dat$Microbiome %in% combined_dat$taxon)
dat$Exposure %>% unique()
combined_dat$phenotype %>% unique()

dat$Exposure[dat$Exposure == "Biofuel"] <- "Biofuel_exposure"
dat$Exposure[dat$Exposure == "Occupation"] <- "Occupational_pollution"
dat$Exposure[dat$Exposure == "PM25"] <- "year2pm25"
dat$Exposure[dat$Exposure == "SHS"] <- "SHS_binary"
dat$Exposure[dat$Exposure == "Smoking"] <- "Smoking_binary"


plotDat_hm <- merge(dat, combined_dat %>% select(phenotype, taxon, Zscore),
                    by.x = c("Exposure","Microbiome"), by.y = c("phenotype","taxon"))

plotDat_hm$mf.label <- sapply(plotDat_hm$Microbiome, 
                              function(x) config_element$labels_full[which(config_element$names == x)])

mf.lvl <- config_element$labels_full[!is.na(config_element$labels_full)]
all(plotDat_hm$mf.label %in% mf.lvl) 


plotDat_hm$mf.label <- factor(plotDat_hm$mf.label,
                              levels = mf.lvl)
head(plotDat_hm)

max.limit = 5
plotDat_hm$PlotValue <- sapply(1:nrow(plotDat_hm),
                               function(i){
                                 if(plotDat_hm$Pval[i] > 1){
                                   NA
                                 }else{
                                   #limits: 
                                   if( plotDat_hm$Zscore[i] > max.limit) {
                                     max.limit
                                   }else if(plotDat_hm$Zscore[i] < (-1 * max.limit)){
                                     -1 * max.limit
                                   }else{
                                     plotDat_hm$Zscore[i]
                                   }
                                 }
                               })




library(ggplot2)
pdf("exposure_microbiome_lm_hm.pdf", width = 4,height = 10)
ggplot(plotDat_hm) +
  geom_tile(aes(x=Exposure, y=mf.label, fill=PlotValue),color="#EFEFEF")+
  scale_fill_gradient2(low = '#2F7E77', high = '#9C6625', mid = '#F5F5F4', midpoint = 0,na.value="white") +
  #scale_fill_continuous_divergingx(palette="PuOr", mid=0, na.value="white") + 
  #  geom_text(aes(x=Feature, y=Taxa, label=Text, color=Text)) +
  #  scale_y_discrete(labels = Taxa.labels$Taxa.label) +
  scale_color_manual(values = c("black","white")) +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text.x = element_text(angle = 90))
dev.off()



