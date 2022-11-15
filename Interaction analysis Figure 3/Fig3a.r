library(data.table)
library(dplyr)


load("../_data/microbialFeatures_standardized.RData")

# read meta data 
load("../_data/meta_phenotypes.RData")
meta <- merge(
  merge(meta.num, meta.cat,by="SampleID"),
  meta.geo_name %>% select(SampleID, District), 
  by="SampleID"
) 
meta$District <- as.factor(meta$District)

# Y: lung function
Ys <- c("CAT_score","A_FEV1_FVC_Post")
Y_df <- meta %>% select(SampleID, all_of(Ys))


# Exposure: 
Expos <- c( "Biofuel_exposure","Occupational_pollution","SHS_binary",
            "Smoking_binary","year2pm25")
Expos_df <- meta %>% select(SampleID, all_of(Expos))
Expos_df <- cbind.data.frame("SampleID" = Expos_df$SampleID,
                             sapply(Expos_df[,-1], function(x) as.numeric(sub( "Y", 1, sub("N",0,x))) ),
                             stringsAsFactors=F) 


# df of covariables, rowname=sample, colnames=variable 
covariables <- c("Age","BMI", "Gender","District","Medication")
covar <- meta %>% select(SampleID, all_of(covariables)) %>% tibble::column_to_rownames("SampleID") 
covar$Gender <- as.integer(sub( "F", 1, sub("M",0,covar$Gender)))



# interaction ------------------------------------
if(any(grepl("District", covariables))) i.pval = 12 else i.pval = 7

Pvalue_res <- NULL
for(y in Ys){
  # y=Ys[1]
  writeLines(paste0("lung function:" , y, " , ", which(Ys == y), " out of ", length(Ys)))
  for(ep in Expos){
    # ep = Expos[1]  
    writeLines(paste0("exposure: " , ep, " , ", which(Expos == ep), " out of ", length(Expos)))
    
    meta.tmp <- merge(meta %>% select(SampleID, all_of(y)),
                      Expos_df %>% select(SampleID, all_of(ep)),
                      by="SampleID")
    
    for(df.name in c("bact.features.st", "fung.features.st", "metagMod.features.st", 
                     "metagTaxa.features.st", "ARG.features.st", "VF.features.st")){
      # df.name = "bact.features.st"
      
      microbFtr_df <- eval(parse(text = df.name)) %>% as.data.frame()
      colnames(microbFtr_df) <- gsub("\\W", "_", colnames(microbFtr_df))
      writeLines(paste0("Microbial feature type: " , df.name))
      
      for(i in 1:ncol(microbFtr_df)){
        # i=1
        # writeLines(paste0(i, " out of ", ncol(microbFtr_df), " microbial features."))
        
        mf = colnames(microbFtr_df)[i]
        mf_df.tmp <- microbFtr_df %>% select(mf)
        
        dat <- merge(merge(meta.tmp, mf_df.tmp, by.x = "SampleID", by.y = 0),
                     covar,
                     by.x="SampleID",by.y=0)
        
        dat <- dat[complete.cases(dat),]
        
        fml <- paste0(y, " ~ ", paste(colnames(covar), collapse = " + "), " + ", ep, " + ", mf, " + ",
                      ep, " * ", mf )
        m <- lm(as.formula(fml), data = dat)
        
        #m.res <- summary(m)$coefficients %>% as.data.frame()
        #m.res$`Pr(>|t|)`[i.pval]
        
        an1 = anova(m)
        F.stat = an1[7,4] 
        Pvalue = an1[7,5]  
        
        res_vec <- c("LungFunction" = y, "Exposure" = ep,"MicrobType" = df.name, "MicrobFeature" = mf,
                     "F.value" = F.stat,"P.interaction" = Pvalue)
        Pvalue_res <- bind_rows(Pvalue_res, res_vec)
        
      }
      
    } # loop through microbe features
  }# loop through exposures
}# loop through lung functions


write.csv(Pvalue_res, file = "Interaction_exposure.microbial_lungFunction.csv", quote = F, row.names = F)
