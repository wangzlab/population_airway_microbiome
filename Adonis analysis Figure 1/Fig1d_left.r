
# Left panel: top 10 DAFs ============================

library(data.table)
library(dplyr)
load("../_data/meta_phenotypes.RData")
meta <- merge(merge(meta.cat, meta.num, by = "SampleID"),
              meta.geo_name, by = "SampleID")

vars_toAnalyze <- c("District","Smoking_binary","Gender","year2pm25", "Age","A_FEV1_FVC_Post",
                    "Biofuel_exposure","FEV1pred_post","COPD_diagnosis","Hypertension","Medication",
                    "Occupational_pollution","Rhinitis","CAT_score","Allergic_rhinitis","SHS_binary",
                    "Osteoporosis","BMI","Phlegm","Stroke","Bronchiectasis","Cough","Anemia","OSAS","Emphysema",
                    "Wheeze","Asthma","Depression","Tuberculosis","HF","Gastroesophageal_reflux","Dyspnea","CHD",
                    "Chronic_Bronchitis","Diabetes","Arrhythmia","IPF","Other_cancer")

all(vars_toAnalyze %in% colnames(meta))

meta.analyze <- meta %>% select(all_of(c("SampleID", vars_toAnalyze)))


catVars <- names(which(sapply(meta.analyze, class) == "character"))
catVars <- catVars[catVars != "SampleID"]
numVars <- names(which(sapply(meta.analyze, class) == "numeric"))


NumVars_site_statsDf <- NULL
for(v in numVars){
  # v=numVars[1]
  dat.test <- meta.analyze %>% select(all_of(c("SampleID", v)))
  dat.test <- dat.test[complete.cases(dat.test),]
  rownames(dat.test) <- NULL
  dat.test <- dat.test %>% tibble::column_to_rownames("SampleID")
  colnames(dat.test) <- "Var"
  
  sites <- meta.analyze %>% select(SampleID, District)
  sites <- sites[match( rownames(dat.test), sites$SampleID),]
  colnames(sites)[which(colnames(sites) == "District")] <- "Site"
  
  # kruskal
  k = kruskal.test(dat.test$Var ~ Site, data = sites)
  kruskal_c <- c("Variable" = v, "Evaluation"="Kruskal-Wallis",
                 "Value"= unname(k$statistic), "pvalue"=unname(k$p.value))
  
  NumVars_site_statsDf <- bind_rows(NumVars_site_statsDf,  kruskal_c) # 不算adonis
  
}# variables
NumVars_site_statsDf$pvalue <- as.numeric(NumVars_site_statsDf$pvalue)
NumVars_site_statsDf$pajd <- p.adjust(NumVars_site_statsDf$pvalue) 
#write.table(NumVars_site_statsDf, file = "7_1.numVars.Site_kruskal_batch2.txt", quote=F,  row.names = F, sep = "\t")




# categorical variables using Chi-Square test  ----
CatVars_site_statsDf <- NULL
catVars <- catVars[!catVars %in% c("District")]
for(v in catVars){
  
  dat.test <- meta.analyze %>% dplyr::select(SampleID, all_of(c(v, "District") ))
  dat.test <- dat.test[complete.cases(dat.test),]
  
  colnames(dat.test)[colnames(dat.test) == v] <- "Var"
  colnames(dat.test)[colnames(dat.test) == "District"] <- "Site"
  
  c <- chisq.test(dat.test$Var, dat.test$Site, correct = F)
  
  res_c <- c("Variable" = v, "Evaluation"="X-squared",
             "Value"= unname(c$statistic), "pvalue"=unname(c$p.value))
  
  CatVars_site_statsDf <- bind_rows(CatVars_site_statsDf, res_c)
}# variables

CatVars_site_statsDf$pvalue <- as.numeric(CatVars_site_statsDf$pvalue)
CatVars_site_statsDf$pajd <- p.adjust(CatVars_site_statsDf$pvalue) 

#write.table(CatVars_site_statsDf, file = "7_1.catVariables.site_chisq_batch2.txt", sep = '\t', row.names = F, quote = F)


districtAssocation <- bind_rows(NumVars_site_statsDf, CatVars_site_statsDf)



