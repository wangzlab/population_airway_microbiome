library(data.table)
library(dplyr)
library(ggplot2)

# ===============================================================================
# Figure 4b   
# ===============================================================================

AMHI <- fread("AMHI.txt",data.table = F)

Health <- fread("../_data/meta_airwayHealth_overallHealth.txt")
Health <- cbind.data.frame(
  SampleID = Health$`#NAME`,
  sapply(Health[,2:3], function(x) as.logical(x))
) %>% tibble::column_to_rownames("SampleID") %>%
  mutate(AirwayHealth = !Airway_disease,
         GeneralHealth = !All_disease) %>%
  select(AirwayHealth, GeneralHealth)


load("../_data/meta_phenotypes.RData")
colnames(meta.cat) <- sub("\\W","_", colnames(meta.cat))
airwayDiseases <- c("Allergic_rhinitis","Asthma","Chronic_Bronchitis","COPD",
                    "Rhinitis","Tuberculosis"  )

n.air.diseases <- c("Anemia","Arrhythmia","CHD",
                    'Diabetes',"Gastroesophageal_reflux","Hypertension","Osteoporosis")

RespiratorySymptom <- c("Cough","Dyspnea","Phlegm","Wheeze")

# switch Y/N to logical
meta.cat[2:ncol(meta.cat)] <- 
  sapply(meta.cat[2:ncol(meta.cat)],
       function(x) as.logical(sub("N",F,sub("Y",T, x))))
meta.cat <- meta.cat %>% tibble::column_to_rownames("SampleID")

#calculate No symptom samples:
meta.symptoms <- meta.cat %>% tibble::column_to_rownames("SampleID") %>% select(all_of(RespiratorySymptom))
meta.noSymp <- data.frame("NoSymptom" = apply(meta.symptoms, 1, function(x) !any(x)))


# CAT score 
meta.CAT <- meta.num %>% 
  mutate(largeCAT = cut(CAT_score, breaks=c(-Inf, 9, Inf), labels=c(F, T))) %>%
  mutate(smallCAT = cut(CAT_score, breaks=c(-Inf, 9, Inf), labels=c(T, F))) %>%
  select(SampleID, largeCAT, smallCAT) %>%
  tibble::column_to_rownames("SampleID")
head(meta.CAT)

meta.CAT <- cbind.data.frame(
  SampleID = rownames(meta.CAT),
  sapply(meta.CAT, function(x) as.logical(x))
) %>% tibble::column_to_rownames("SampleID")


common.sps <- intersect(intersect(rownames(Health), rownames(meta.cat)),
                        intersect(rownames(meta.noSymp), rownames(meta.CAT)))

Meta <- bind_cols(Health[common.sps, ],
                  meta.cat[common.sps, c(airwayDiseases, n.air.diseases)],
                  meta.noSymp[common.sps, , drop=F],
                  meta.symptoms[common.sps, ],
                  meta.CAT[common.sps, ])

tmp <- Meta %>%  
  tibble::rownames_to_column("SampleID") %>%
  reshape2::melt(id.vars ="SampleID") %>%
  filter(value) 

plotDat <- merge(tmp, AMHI, by="SampleID")
plotDat$variable <- factor(plotDat$variable, levels = colnames(Meta))

colorRampPalette(brewer.pal(11, "Set3"))(22)

Fig4c <-ggplot(plotDat, aes(x=variable, y=AMHI))  +
  geom_hline(yintercept = 0, linetype="dashed", color="gray", lwd=1) +
  geom_violin(aes(x=variable, y=AMHI, fill=variable), alpha=0.6, trim = F) +
  geom_boxplot(aes(x=variable, y=AMHI), width=0.1, outlier.shape = NA) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(11, "Set3"))(22)) +
  theme_bw() + theme(panel.grid = element_blank(),
                     legend.position = "none",
                     axis.text.x = element_text(angle = 90))  +      
  stat_compare_means(label = "p.signif", method =  "wilcox.test",
                     ref.group = "AirwayHealth",  # manually change to "GeneralHealth" and "NoSymptom" and revise stars accordingly
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                        symbols = c("***", "**", "*", "+", "ns")))               

Fig4c

# n for each condition to add into axis x text:
plotDat %>% group_by(variable) %>% dplyr::summarise(n=n())


# ===============================================================================
# Figure 4d  
# ===============================================================================
colnames(meta.cat)
Exposures <- c("Smoking_binary","SHS_binary","Biofuel_exposure","Occupational_pollution")

conditions <- append(
  
  append(
  as.list(as.data.frame(combn(Exposures,m=1))),
  as.list(as.data.frame(combn(Exposures,m=2)))),
  append(
    as.list(as.data.frame(combn(Exposures,m=3))),
    as.list(as.data.frame(combn(Exposures,m=4))))
)
names(conditions) <- paste0("condition", seq(1:length(conditions)))


Meta.exposures <- NULL
# the first column is no exposure:
tmp <- meta.cat %>% select(all_of(Exposures))
meta.noExpo <- data.frame("noExposure"=apply(tmp, 1, function(x) !any(x)))
Meta.exposures <- bind_rows(Meta.exposures, meta.noExpo)
# the other columns are combinations of exposures
for(i in 1:length(conditions)){
  # i=1
   cond = conditions[[i]]
   
   tmp <- meta.cat %>% select(all_of(cond))
   remains <- meta.cat %>% select(all_of(Exposures)) %>% select(-all_of(cond))
   
   meta.cond <- data.frame(apply(tmp, 1, all) & apply(remains, 1, function(x) !any(x)))
     
   Meta.exposures <- bind_cols(Meta.exposures, meta.cond)
   colnames(Meta.exposures)[ncol(Meta.exposures)] <- paste(cond, collapse = ";")
  
}


tmp <- Meta.exposures %>%  
  tibble::rownames_to_column("SampleID") %>%
  reshape2::melt(id.vars ="SampleID") %>%
  filter(value) 

plotDat.exps <- merge(tmp, AMHI, by="SampleID")
plotDat.exps$variable <- factor(plotDat.exps$variable,
                                levels = colnames(Meta.exposures))

plotDat.exps$numExpos <-
  as.character(sapply(as.character(plotDat.exps$variable), 
                      function(x) length(strsplit(x,";",fixed = T)[[1]])) )
plotDat.exps$numExpos[plotDat.exps$variable == "noExposure"] <- 0

library(ggpubr)

Fig4d <- ggplot(plotDat.exps,aes(x=variable, y=AMHI)) +
  geom_hline(yintercept = 0, linetype="dashed", color="gray",lwd=1)+
  geom_violin(aes(x=variable, y=AMHI, fill=numExpos), alpha=0.6, trim = F) +
  geom_boxplot(aes(x=variable, y=AMHI), width=0.1, outlier.shape = NA) +
  scale_fill_manual(values = c("#EDEDEC","#EBCBD1","#D997A2","#C66475","#B7364B")) +
  theme_bw() + theme(panel.grid = element_blank(),
                     legend.position = "none",
                     axis.text.x = element_text(angle = 90)) +      
  stat_compare_means(label = "p.signif", method =  "wilcox.test",
                     ref.group = "noExposure",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                                        symbols = c("***", "**", "*", "+", "ns")))               
Fig4d

# n for each condition:
plotDat.exps %>% group_by(variable) %>% dplyr::summarise(n=n())
