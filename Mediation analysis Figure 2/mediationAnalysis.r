library(data.table)
library(dplyr)
library(mediation)


# load data ============================================
load("../_data/microbialFeatures_standardized.RData") 
load("../_data/meta_phenotypes.RData")
meta <- merge(meta.cat, meta.num, by="SampleID")

# exposures are the treators: 
Exposures <- c("Biofuel_exposure","Occupational_pollution","Smoking_binary","SHS_binary","year2pm25")

# lung functions:
LungFuns <- c("A_FEV1_FVC_Post","FEV1pred_post","CAT_score")


# Covariates:
Covars <- c("Age","BMI","Gender","Medication")
covar_df <- meta %>% dplyr::select(SampleID, all_of(Covars))
covar_df$Gender <- as.integer(sub( "F", 1, sub("M",0,covar_df$Gender)))
covar_df$Medication <- as.integer(sub( "Y", 1, sub("N",0, covar_df$Medication)))


# mediation analysis =================
Mediation.forward <- NULL
Mediation.reverse <- NULL

for(exp in Exposures){
  #exp=Exposures[1]
  meta.expo <- meta %>% dplyr::select(SampleID, all_of(exp))
  meta.expo$Biofuel_exposure <- as.integer(sub("N",0,sub("Y",1,meta.expo$Biofuel_exposure)))
  
  for(lf in LungFuns){
     #lf = LungFuns[1]
     meta.lf = meta %>% dplyr::select(SampleID, all_of(lf))
     
     for(mf.dn in c("bact.features.st", "fung.features.st", "metagMod.features.st", "ARG.features.st", "VF.features.st")){
       #mf.dn = "bact.features.st"
       
       mf.data <- eval(parse(text = mf.dn))
       
       for(i in 1:ncol(mf.data)){
         mf.dat.tmp <- mf.data[,i, drop=F]
         mf = colnames(mf.data)[i]
         
         dat <- merge(merge(merge(meta.expo, meta.lf, by = "SampleID"),
                            covar_df, by="SampleID"),
                      mf.dat.tmp, by.x="SampleID", by.y = 0)
         
         dat <- dat[complete.cases(dat),]
         
         # forward mediation --------------------
         id = paste(exp, mf, lf, sep = "|")
         
         model.m = lm(as.formula( paste(mf, " ~ ", exp, " + ", 
                                        paste(colnames(dat)[!(colnames(dat) %in% c("SampleID",lf, mf, exp))], collapse = " + "),
                                        sep = "") ), data = dat)
         
         model.y = lm(as.formula( paste(lf, " ~ ", exp," + ", mf," + ", 
                                        paste(colnames(dat)[!(colnames(dat) %in% c("SampleID",lf, mf, exp))], collapse = " + "),
                                        sep = "") ), data = dat)
         
         summary = summary(mediate(model.m,model.y,treat=exp,mediator=mf,boot=F,sims=1000))
         res <- capture.output(summary,append=FALSE)
         
         tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
         tmp <- tmp[tmp != "" & tmp!="."]
         tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
         ACME.p <- tmp[length(tmp)]
         
         tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
         tmp <- tmp[tmp != "" & tmp!="."]
         tmp <- tmp[!grepl("*",tmp,fixed = T) ]
         ADE.p <- tmp[length(tmp)]
         
         tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
         tmp <- tmp[tmp != "" ]
         i_str = which(grepl("Mediated", tmp))
         prop.mediated <- tmp[(i_str + 1)]
         
         forw_vec = c(id, ACME.p, ADE.p, prop.mediated)
         names(forw_vec) <- c("Treat_Mediator_Y", "ACME.p", "ADE.p", "prop.mediated")
         
         Mediation.forward <- bind_rows(Mediation.forward, forw_vec)
         
         # reverse mediation --------------------
         id = paste(exp, lf, mf, sep = "|")
         
         model.m = lm(as.formula( paste(lf, " ~ ", exp, " + ", 
                                        paste(colnames(dat)[!(colnames(dat) %in% c("SampleID",lf, mf, exp))], collapse = " + "),
                                        sep = "") ), data = dat)
         
         model.y = lm(as.formula( paste(mf, " ~ ", exp," + ", lf," + ", 
                                        paste(colnames(dat)[!(colnames(dat) %in% c("SampleID",lf, mf, exp))], collapse = " + "),
                                        sep = "") ), data = dat)
         
         summary = summary(mediate(model.m,model.y,treat=exp,mediator=lf,boot=F,sims=1000))
         res <- capture.output(summary,append=FALSE)
         
         #sub( "^()\\s", "\\1", res[7])
         tmp <-  base::strsplit(res[grep("ACME",res)],"\\s")[[1]]
         tmp <- tmp[tmp != "" & tmp!="."]
         tmp <- tmp[!grepl("*",tmp,fixed = T)] # remove stars in case the last element is star
         ACME.p <- tmp[length(tmp)]
         
         tmp <- base::strsplit(res[grep("ADE",res)],"\\s")[[1]]
         tmp <- tmp[tmp != "" & tmp!="."]
         tmp <- tmp[!grepl("*",tmp,fixed = T) ]
         ADE.p <- tmp[length(tmp)]
         
         tmp <- base::strsplit(res[grep("Prop. Mediated", res)],"\\s")[[1]]
         tmp <- tmp[tmp != "" ]
         i_str = which(grepl("Mediated", tmp))
         prop.mediated <- tmp[(i_str + 1)]
         
         revers_vec = c(id, ACME.p, ADE.p, prop.mediated)
         names(revers_vec) <- c("Treat_Mediator_Y", "ACME.p", "ADE.p", "prop.mediated")
         
         Mediation.reverse <- bind_rows(Mediation.reverse, revers_vec)
         
         
         
       }# loop through individual microbial features
     }# loop through microbial data frames
  }# loop through lung functions
}# loop through exposures