# lm(microbialFeature ~ covariates + exposure ) in Neiss_High and Neiss_low groups 
# microbialFeatures : alpha diversity of bacteria, fungi, metag taxon, ko, ARG, VF
# abundance of ARG and VF ( zscore )


library(data.table)
library(dplyr)
library(foreach)

load("../_data/microbialFeatures_standardized.RData")


# read meta data 
load("../_data/meta_phenotypes.RData")
meta <- merge(merge(meta.num, meta.cat, by="SampleID"), meta.geo_name %>% select(SampleID, District), by="SampleID")
vars_toAnalyze <- c("Occupational_pollution")
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


# Neisseria high and low samples
grouping_df <- fread("occupation_PC1_meta.txt", data.table=F)


# matrix of features：rowname=sample, colnames=features

for(df.name in c("bact.features.st", "fung.features.st", "VF.features.st","ARG.features.st",
                 "metagMod.features.st", "metagTaxa.features.st")){
  # df.name = "bact.features.st" 
  ftrs_transformed <- eval(parse(text = df.name))
  
  
  # belowing are codes from DMP paper: https://github.com/GRONINGEN-MICROBIOME-CENTRE/DMP 
  # Run multivariate models, multi-thread implementation
  # ======================================================
  # prep parallelization
  writeLines(paste0("parallel calculation of lm started for:  ",df.name))
  # registerDoSEQ() # debug mode = single threaded
  threads = 12
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
  
  foreach(j = 1:ncol(ftrs_transformed),.combine = rbind) %:%  # parallel implementation
    
    #debug output/mode
    #foreach(j = 1:2,.combine = rbind) %dopar% {  #debug/test
    #print(i)
    #print(j)
    
    foreach(grp = c("Neis_high","Neis_low"),.combine = rbind) %dopar% {
      # grp="Neis_high"
      
      predictors_2grps = data.frame(covar[!is.na(pheno[,i]),],
                                    model.matrix(
                                      as.formula(paste0("~ ",colnames(pheno)[i])),data = pheno)[,-1,drop = F])
      
      
      sps <- grouping_df$`#NAME`[which(grouping_df$Group == grp)]
      predictors <- predictors_2grps[sps,]
      
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
            group = grp,
            Nsamples = nrow(cleaned_data),
            levels = 
              if(class(cleaned_data[,ncol(cleaned_data)]) == "factor" | length(table(cleaned_data[,ncol(cleaned_data)]))==2) {
                paste(collapse=":",names(table(cleaned_data[,ncol(cleaned_data)])))
              } else "Not Applicable",
            levels_SampleSize = 
              if(class(cleaned_data[,ncol(cleaned_data)]) == "factor" | length(table(cleaned_data[,ncol(cleaned_data)]))==2){
                paste(collapse= ":",table(cleaned_data[,ncol(cleaned_data)])) 
              } else "Not Applicable",
            effect.size =
              if(class(cleaned_data[,ncol(cleaned_data)]) == "factor") {
                paste(collapse = ":",c(0,round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))])))
              } else round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))]) ,
            
            R2 = summary(s1)$r.squared - summary(s0)$r.squared,
            F.stat = an1[2,5],
            Pvalue = an1[2,6]
          )
        }
        
        
        # if no covariates are considered
        if(F){
          s1 = lm(
            as.formula(paste("ftrs.cleaned[,j] ~ ",colnames(cleaned_data)[ncol(cleaned_data)])),
            data = cleaned_data
          )
          
          # compare models
          an1 = anova(s1)
          output = data.frame(
            phenotype = colnames(pheno)[i],
            taxon = colnames(ftrs.cleaned)[j],
            group = grp,
            Nsamples = nrow(cleaned_data),
            levels = 
              if(class(cleaned_data[,ncol(cleaned_data)]) == "factor" | length(table(cleaned_data[,ncol(cleaned_data)]))==2) {
                paste(collapse=":",names(table(cleaned_data[,ncol(cleaned_data)])))
              } else "Not Applicable",
            levels_SampleSize = 
              if(class(cleaned_data[,ncol(cleaned_data)]) == "factor" | length(table(cleaned_data[,ncol(cleaned_data)]))==2){
                paste(collapse= ":",table(cleaned_data[,ncol(cleaned_data)])) 
              } else "Not Applicable",
            effect.size =
              if(class(cleaned_data[,ncol(cleaned_data)]) == "factor") {
                paste(collapse = ":",c(0,round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))])))
              } else round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))]) ,
            
            R2 = summary(s1)$r.squared ,
            F.stat = an1[1,4],
            Pvalue = an1[1,5]
          )
        }
        
        
        #add covariates
        output
      }# condition: nrow(cleaned_data)>3 etc
      # }# loop through Neisseria high and low 
      
      
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
  write.table(result_ftrs, file = paste0("microbial_occupation_association/lmRes_",df.name,"_occupation.txt"), quote = F, row.names = F, sep = "\t" )
  
  # export another format: directly compare Neiss_high and low groups
  totalMFeatures <- c(
    result_ftrs$taxon[grepl("^bact",result_ftrs$taxon )] %>% as.character() %>%  unique(),
    result_ftrs$taxon[grepl("^fun",result_ftrs$taxon )] %>% as.character() %>%  unique(),
    result_ftrs$taxon[grepl("^metag",result_ftrs$taxon )] %>% as.character() %>%  unique(),
    result_ftrs$taxon[grepl("^Meta",result_ftrs$taxon )]  %>% as.character() %>%  unique(),
    result_ftrs$taxon[grepl("^ARG_",result_ftrs$taxon )]  %>% as.character() %>%  unique(),
    result_ftrs$taxon[grepl("^VF_",result_ftrs$taxon )]  %>% as.character() %>%  unique()
  )
  
  
  test <- bind_rows( result_ftrs %>%
                       filter(!taxon %in% totalMFeatures) %>%
                       group_by(phenotype, group) %>%
                       mutate(padj = p.adjust(as.numeric(Pvalue))),
                     result_ftrs %>% filter(taxon %in% totalMFeatures)) 
  
  result_ftrs.w <- 
    merge(test %>% reshape2::dcast(phenotype + taxon ~ group, value.var = "Pvalue"),
          test %>% reshape2::dcast(phenotype + taxon ~ group, value.var = "padj"),
          by=c("phenotype","taxon")) 
  
  colnames(result_ftrs.w) <- gsub( "\\.y", "\\.Padj",  gsub("\\.x", "\\.Pvalue",  colnames(result_ftrs.w) ))
  
  write.table(result_ftrs.w, file = paste0("microbial_occupation_association/lmRes_wide_",df.name,"_occupation.txt"), quote = F, row.names = F, sep = "\t" )
  
  writeLines(paste0(df.name, " has finished and results exported."))
}# loop through data types

