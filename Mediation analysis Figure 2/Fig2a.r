library(data.table)
library(dplyr)
library(foreach)

# Load data ==================================================
load("../_data/microbialFeatures_standardized.RData")

# read meta data 
load("../_data/meta_phenotypes.RData")
meta <- merge(merge(meta.num, meta.cat, by="SampleID"), meta.geo_name %>% select(SampleID, District), by="SampleID")
vars_toAnalyze <- c("Smoking_binary","year2pm25","Biofuel_exposure","Occupational_pollution","SHS_binary",
                    "A_FEV1_FVC_Post","FEV1pred_post","CAT_score")
all(vars_toAnalyze %in% colnames(meta))



# df of covariables, rowname=sample, colnames=variable 
covariables <- c("Age","BMI", "Gender", "District","Medication")
covar <- meta %>% select(SampleID, all_of(covariables)) %>% tibble::column_to_rownames("SampleID") 
covar$Gender <- as.integer(sub( "F", 1, sub("M",0,covar$Gender)))


# df of phenotypes excluding covariables,  rowname=sample, colnames=phenotypes 
pheno <- meta %>% 
  select(SampleID, all_of(vars_toAnalyze)) %>% 
  tibble::column_to_rownames("SampleID")  
pheno[pheno == "N"] <- 0
pheno[pheno == "Y"] <- 1
colnames(pheno) <- gsub("\\W","_", colnames(pheno))



# lm to calculate association between phenotypes and microbial features =================
# matrix of features：rowname=sample, colnames=features

for(df.name in c("bact.features.st", "fung.features.st", "VF.features.st","ARG.features.st",
                 "metagMod.features.st", "metagTaxa.features.st")){
  # df.name = "metagMod.features.st" 
  ftrs_transformed <- eval(parse(text = df.name))
  
  
  # belowing are codes from DMP paper: https://github.com/GRONINGEN-MICROBIOME-CENTRE/DMP 
  # Run multivariate models, multi-thread implementation
  # ======================================================
  # prep parallelization
  writeLines(paste0("parallel calculation of lm started for:  ",df.name))
  # registerDoSEQ() # debug mode = single threaded
  threads = 4
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  # debug: timer
  # t1 <- Sys.time()
  # loop over all phenotypes
  result_ftrs = foreach(i = 1:ncol(pheno),.combine = rbind) %:%
    
    #result_ftrs = foreach(i = 1:50,.combine = rbind) %:%     # debug/test
    
    # loop over all features of requested type (ftrs, VFs, PWYs or CARDs)
    # ==========================================================================
  # foreach(j = 1:ncol(ftrs_transformed),.combine = rbind) %do% { # single-threaded implementation for debug
  
  foreach(j = 1:ncol(ftrs_transformed),.combine = rbind) %dopar% {  # parallel implementation
    
    #debug output/mode
    #foreach(j = 1:2,.combine = rbind) %dopar% {  #debug/test
    #print(i)
    #print(j)
    predictors = data.frame(covar[!is.na(pheno[,i]),],
                            model.matrix(
                              as.formula(paste0("~ ",colnames(pheno)[i])),data = pheno)[,-1,drop = F])
    
    cleaned_data = predictors[complete.cases(predictors),]
    rn <- rownames(cleaned_data)
    rn <- rn[rn %in% rownames(ftrs_transformed)]
    ftrs.cleaned = ftrs_transformed[rn,]
    cleaned_data = cleaned_data[rn,]
    if (nrow(cleaned_data) > 3) {
      # debug: print model
      #print(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data))))
      
      
      # lm as applied in DMP scripts: use lm because Y is features, always continuous variable with normal distribution
      # if covariates are considered
      if(T){
        # make model
        s1 = lm(
          as.formula(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)))),
          data = cleaned_data
        )
        # debug: print model
        #print(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)[1:ncol(covar)])))
        
        # make model with extra covariates
        s0 = lm(
          as.formula(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)[1:ncol(covar)]))),
          data = cleaned_data
        )
        
        # compare models
        an1 = anova(s1,s0)
        output = data.frame(
          phenotype = colnames(pheno)[i],
          taxon = colnames(ftrs.cleaned)[j],
          Nsamples = nrow(cleaned_data),
          levels = if(class(pheno[,i]) == "factor") paste(collapse=":",levels(pheno[,i])) else "Not Applicable",
          levels_SampleSize = 
            if(class(pheno[,i]) == "factor" | length(table(pheno[,i]))==2) paste(collapse= ":",table(pheno[,i])) else "Not Applicable",
          effect.size =
            if(class(pheno[,i]) == "factor") {
              paste(collapse = ":",c(0,round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))])))
            } else round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))]) ,
          
          R2 = summary(s1)$r.squared - summary(s0)$r.squared,
          F.stat = an1[2,5],
          Pvalue = an1[2,6]
        )
      }
      
      
      #add covariates
      output
    }# condition: nrow(cleaned_data)>3 etc
  }# parallel inside
  # debug
  #t2 <-  Sys.time()
  writeLines("ftrs_done")
  # debug
  #print(t2-t1)
  rownames(result_ftrs) <- NULL
  #result_ftrs$FDR <- p.adjust(result_ftrs$Pvalue,method = "BH")
  #write.table(,file = phenoOut,row.names=F,quote = F,sep="\t")
  #close parallelization
  #registerDoSEQ()
  on.exit(stopCluster(cl))
  
  # return results
  result_ftrs
  
  # export results
  write.table(result_ftrs, file = paste0("lmRes_",df.name,"_withCovariates.txt"), quote = F, row.names = F, sep = "\t" )
  writeLines(paste0(df.name, " has finished and results exported."))
}# loop through data types

# padj for each results data frame ========================================
res <-fread("lmRes_bact.features.st_withCovariates.txt")
test<-
  bind_rows( res %>%
               filter(!taxon %in% c("bact.alpha.st","bact.pco1")) %>%
               group_by(phenotype) %>%
               mutate(padj = p.adjust(Pvalue)),
             res %>% filter(taxon %in% c("bact.alpha.st","bact.pco1"))) %>%
  select(-levels,-levels_SampleSize) %>%
  arrange(Pvalue)
write.table(test, file = "lmRes_bact.features.st_withCovariates_padj.txt", quote = F, row.names = F, sep = "\t")

res <-fread("lmRes_ARG.features.st_withCovariates.txt")
test<-
  bind_rows( res %>%
               filter(!taxon %in% c("ARG.alpha.st","ARG.pco1","ARG.abund.st")) %>%
               group_by(phenotype) %>%
               mutate(padj = p.adjust(Pvalue)),
             res %>% filter(taxon %in% c("ARG.alpha.st","ARG.pco1","ARG.abund.st"))) %>%
  select(-levels,-levels_SampleSize) %>%
  arrange(Pvalue)
write.table(test, file = "lmRes_ARG.features.st_withCovariates_padj.txt", quote = F, row.names = F, sep = "\t")

res <-fread("lmRes_fung.features.st_withCovariates.txt")
test<-
  bind_rows( res %>%
               filter(!taxon %in% c("fung.alpha.st","fung.pco1")) %>%
               group_by(phenotype) %>%
               mutate(padj = p.adjust(Pvalue)),
             res %>% filter(taxon %in% c("fung.alpha.st","fung.pco1"))) %>%
  select(-levels,-levels_SampleSize) %>%
  arrange(Pvalue)
write.table(test, file = "lmRes_fung.features.st_withCovariates_padj.txt", quote = F, row.names = F, sep = "\t")

res <-fread("lmRes_metagMod.features.st_withCovariates.txt")
test<-
  bind_rows( res %>%
               filter(!taxon %in% c("ko.alpha.st","MetagMod.pco1.st")) %>%
               group_by(phenotype) %>%
               mutate(padj = p.adjust(Pvalue)),
             res %>% filter(taxon %in%  c("ko.alpha.st","MetagMod.pco1.st"))) %>%
  select(-levels,-levels_SampleSize) %>%
  arrange(Pvalue)
write.table(test, file = "lmRes_metagMod.features.st_withCovariates_padj.txt", quote = F, row.names = F, sep = "\t")

res <-fread("lmRes_metagTaxa.features.st_withCovariates.txt")
test<-
  bind_rows( res %>%
               filter(!taxon %in% c("metagTaxa.alpha.st","metagTaxa.pco1")) %>%
               group_by(phenotype) %>%
               mutate(padj = p.adjust(Pvalue)),
             res %>% filter(taxon %in% c("metagTaxa.alpha.st","metagTaxa.pco1"))) %>%
  select(-levels,-levels_SampleSize) %>%
  arrange(Pvalue)
write.table(test, file = "lmRes_metagTaxa.features.st_withCovariates_padj.txt", quote = F, row.names = F, sep = "\t")

res <-fread("lmRes_VF.features.st_withCovariates.txt")
test<-
  bind_rows( res %>%
               filter(!taxon %in% c("VF.alpha.st","VF.pco1","VF.abund.st")) %>%
               group_by(phenotype) %>%
               mutate(padj = p.adjust(Pvalue)),
             res %>% filter(taxon %in% c("VF.alpha.st","VF.pco1","VF.abund.st"))) %>%
  select(-levels,-levels_SampleSize) %>%
  arrange(Pvalue)
write.table(test, file = "lmRes_VF.features.st_withCovariates_padj.txt", quote = F, row.names = F, sep = "\t")

# plotting ============================================================

remove(list = ls())

library(data.table)
library(dplyr)
plotFiles <- list.files(path = "plotSourceData/", full.names = T, pattern = "select")
combined_dat <- NULL
for(f in plotFiles) {
  
  target = strsplit( basename(f) ,"_", fixed = T)[[1]][1]
  dat<- fread(f) %>% mutate(Target = target) # %>% filter(!Feature %in% ft.rm)
  combined_dat <- bind_rows(combined_dat, dat)
  
}

diversFiles <- list.files(path = "plotSourceData/", full.names = T, pattern = "diversity")
diversity_dat <- NULL
for(f in diversFiles) {
  
  target = sub("\\.txt", "", basename(f))
  dat<- fread(f) %>% mutate(Target = target)# %>% filter(!Feature %in% ft.rm)
  diversity_dat <- bind_rows(diversity_dat, dat)
  
}

combined_dat <- bind_rows(combined_dat, diversity_dat)

# calculate zscore ==========================================
combined_dat$Zscore <- -qnorm(combined_dat$Pvalue/2)
combined_dat$Zscore <- combined_dat$Zscore * sign(combined_dat$effect.size)


# plot together -------------------
plotDat <- combined_dat

# colors
plotDat$direction <- sign(plotDat$effect.size)


max.limit = 5
plotDat$PlotValue <- sapply(1:nrow(plotDat),
                            function(i){
                              if(plotDat$Pvalue[i] > 0.25){
                                NA
                              }else{
                                #limits: 
                                if( abs(plotDat$Zscore[i]) > max.limit) {
                                  max.limit*plotDat$direction[i]
                                }else{
                                  plotDat$Zscore[i]
                                }
                              }
                            })

colnames(plotDat)[colnames(plotDat) == "padj"] <- "FDR"
plotDat$Text <- sapply(1:nrow(plotDat),
                       function(i){
                         if(grepl("diversity", plotDat$Target[i])){
                           if( plotDat$Pvalue[i] > 0.01){
                             NA
                           } else if(plotDat$direction[[i]] == 1) {
                             "+"
                           }else if(plotDat$direction[[i]] == -1) {
                             "-"
                           }
                           
                         }else{
                           if( plotDat$Pvalue[i] > 0.01){
                             NA
                           } else if(plotDat$direction[[i]] == 1) {
                             "+"
                           }else if(plotDat$direction[[i]] == -1) {
                             "-"
                           }
                         }
                         
                       }) 


# sequences -----

colnames(plotDat)[colnames(plotDat) == "phenotype"] <- "Feature"
colnames(plotDat)[colnames(plotDat) == "taxon"] <- "Taxa"

# cluster as a whole
plotDat.w <- plotDat %>% reshape2::dcast(Feature~Taxa, value.var = "Zscore") %>% tibble::column_to_rownames("Feature")
plotDat.w[is.na(plotDat.w)] <-0
if(T){
  library(ggdendro)
  df <- t(plotDat.w) 
  x <- as.matrix(scale(df))
  dd.col <- as.dendrogram(hclust(dist(x), method = "ward.D"))
  col.ord <- order.dendrogram(dd.col)
  
  dd.row <- as.dendrogram(hclust(dist(t(x)), method = "ward.D"))
  row.ord <- order.dendrogram(dd.row)
  
  xx <- scale(df)[col.ord, row.ord] 
  xx_names <- attr(xx, "dimnames") 
  #df <- as.data.frame(xx)
  ddata_x <- dendro_data(dd.row) 
  ddata_y <- dendro_data(dd.col) 
}
xx_names[[1]]
envLevels.whole <- xx_names[[2]] # taxa

#plotDat$Feature <- factor(plotDat$Feature, levels = envLevels.whole)  #定Feature(environmental variables) levels

# cluster inside each data type
df.w <- plotDat %>% filter(Target == "arg") %>%
  reshape2::dcast(Feature~Taxa, value.var = "Zscore") %>% tibble::column_to_rownames("Feature")
if(T){
  df <- t(df.w) 
  x <- as.matrix(scale(df))
  dd.col <- as.dendrogram(hclust(dist(x), method = "ward.D"))
  col.ord <- order.dendrogram(dd.col)
  
  dd.row <- as.dendrogram(hclust(dist(t(x)), method = "ward.D"))
  row.ord <- order.dendrogram(dd.row)
  
  xx <- scale(df)[col.ord, row.ord] 
  xx_names <- attr(xx, "dimnames") 
  #df <- as.data.frame(xx)
  ddata_x <- dendro_data(dd.row) 
  ddata_y <- dendro_data(dd.col) 
} 
ArgLevels <- xx_names[[1]]

df.w <- plotDat %>% filter(Target == "bac") %>%
  reshape2::dcast(Feature~Taxa, value.var = "Zscore") %>% tibble::column_to_rownames("Feature")
if(T){
  df <- t(df.w) 
  x <- as.matrix(scale(df))
  dd.col <- as.dendrogram(hclust(dist(x), method = "ward.D"))
  col.ord <- order.dendrogram(dd.col)
  
  dd.row <- as.dendrogram(hclust(dist(t(x)), method = "ward.D"))
  row.ord <- order.dendrogram(dd.row)
  
  xx <- scale(df)[col.ord, row.ord] 
  xx_names <- attr(xx, "dimnames") 
  #df <- as.data.frame(xx)
  ddata_x <- dendro_data(dd.row) 
  ddata_y <- dendro_data(dd.col) 
}
BacLevels <- xx_names[[1]]

df.w <- plotDat %>% filter(Target == "fun") %>%
  reshape2::dcast(Feature~Taxa, value.var = "Zscore") %>% tibble::column_to_rownames("Feature")
if(T){
  df <- t(df.w) 
  x <- as.matrix(scale(df))
  dd.col <- as.dendrogram(hclust(dist(x), method = "ward.D"))
  col.ord <- order.dendrogram(dd.col)
  
  dd.row <- as.dendrogram(hclust(dist(t(x)), method = "ward.D"))
  row.ord <- order.dendrogram(dd.row)
  
  xx <- scale(df)[col.ord, row.ord]
  xx_names <- attr(xx, "dimnames") 
  #df <- as.data.frame(xx)
  ddata_x <- dendro_data(dd.row) 
  ddata_y <- dendro_data(dd.col) 
} 
FunLevels <- xx_names[[1]]

df.w <- plotDat %>% filter(Target == "module") %>%
  reshape2::dcast(Feature~Taxa, value.var = "Zscore") %>% tibble::column_to_rownames("Feature")
if(T){
  df <- t(df.w) 
  x <- as.matrix(scale(df))
  dd.col <- as.dendrogram(hclust(dist(x), method = "ward.D"))
  col.ord <- order.dendrogram(dd.col)
  
  dd.row <- as.dendrogram(hclust(dist(t(x)), method = "ward.D"))
  row.ord <- order.dendrogram(dd.row)
  
  xx <- scale(df)[col.ord, row.ord] 
  xx_names <- attr(xx, "dimnames") 
  #df <- as.data.frame(xx)
  ddata_x <- dendro_data(dd.row)
  ddata_y <- dendro_data(dd.col)
} 
KoLevels <- xx_names[[1]]

df.w <- plotDat %>% filter(Target == "vf") %>%
  reshape2::dcast(Feature~Taxa, value.var = "Zscore") %>% tibble::column_to_rownames("Feature")
if(T){
  df <- t(df.w) 
  x <- as.matrix(scale(df))
  dd.col <- as.dendrogram(hclust(dist(x), method = "ward.D"))
  col.ord <- order.dendrogram(dd.col)
  
  dd.row <- as.dendrogram(hclust(dist(t(x)), method = "ward.D"))
  row.ord <- order.dendrogram(dd.row)
  
  xx <- scale(df)[col.ord, row.ord] 
  xx_names <- attr(xx, "dimnames") 
  #df <- as.data.frame(xx)
  ddata_x <- dendro_data(dd.row) 
  ddata_y <- dendro_data(dd.col) 
} 
VfLevels <- xx_names[[1]]

# order by feature type
if(F){
  taxaOrders_df <- plotDat %>% 
    select(Target, Taxa) %>%
    unique() %>%
    arrange(Taxa) %>% arrange(Target)
  
  plotDat$Taxa <- factor(plotDat$Taxa, levels = taxaOrders_df$Taxa)
}

# order by feature cluster as a whole
#plotDat$Taxa <- factor(plotDat$Taxa, levels = xx_names[[1]])

#order by feature cluster in each type
alphaLevels <- c("VF.pco1","VF.alpha.st","VF.abund.st","ARG.pco1","ARG.alpha.st","ARG.abund.st",
                 "MetagMod.pco1.st","ko.alpha.st","fung.pco1","fung.alpha.st", "bact.pco1","bact.alpha.st")
plotDat$Taxa <- factor(plotDat$Taxa, levels = c(VfLevels,ArgLevels,KoLevels,FunLevels,BacLevels,alphaLevels)) #定taxon levels


plotDat$Feature <- factor(plotDat$Feature, levels = c("A_FEV1_FVC_Post","FEV1pred_post", "CAT_score",
                                                      "year2pm25","year2pm10","Smoking_binary","Occupational_pollution",
                                                      "Biofuel_exposure","SHS_binary"))

# change ARG abb back to ARG name
ARG.names <- fread("arg_name_mapping.txt", header = F)
plotDat$Taxa.label <- sapply(as.character(plotDat$Taxa), 
                             function(x){
                               if(x %in% ARG.names$V1) ARG.names$V2[which(ARG.names$V1==x)] else x
                             } )
Taxa.labels <- plotDat %>% select(Taxa, Taxa.label) %>% unique() %>% arrange(Taxa)


library(ggplot2)
pdf("Fig2a.microbial_phenotypes_heatmap.pdf",width = 4.5,height = 15)  # ,height = 10
ggplot(plotDat) +
  geom_tile(aes(x=Feature, y=Taxa, fill=PlotValue),color="#EFEFEF") +
  scale_fill_gradient2(low = '#2F7E77', high = '#9C6625', mid = '#F5F5F4', midpoint = 0,na.value="white") +
  #scale_fill_continuous_divergingx(palette="PuOr", mid=0, na.value="white") + 
  geom_text(aes(x=Feature, y=Taxa, label=Text, color=Text)) +
  scale_y_discrete(labels = Taxa.labels$Taxa.label) +
  scale_color_manual(values = c("black","white")) +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.text.x = element_text(angle = 90))
dev.off()



