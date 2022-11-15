library(dplyr)
library(vegan)


load("../_data/microbialFeatures_relAbund_rpkm_notStandardized.RData")
load("../_data/meta_phenotypes.RData")
meta <- merge(meta.geo_name, 
              merge(meta.num, meta.cat, by="SampleID"),
              by="SampleID")
colnames(meta) <- gsub("\\W","_", colnames(meta))


vars2analyze <- c("Smoking_binary","year2pm25","Gender","A_FEV1_FVC_Post","Age", "BMI",
                  "Biofuel_exposure","FEV1pred_post","Rhinitis","Hypertension",
                  "Occupational_pollution","Gastroesophageal_reflux","Phlegm","Emphysema",
                  "Bronchiectasis","Anemia","BMI","CAT_score","Osteoporosis","OSAS","SHS_binary",
                  "Cough","HF","Wheeze","Allergic_rhinitis","Asthma","Depression",
                  "Chronic_Bronchitis","Diabetes","CHD","Stroke","Pulmonary_hypertension",
                  "IPF","Dyspnea","Other_cancer","Tuberculosis","Arrhythmia","Medication") 


Adonis_res <- NULL
for(df.name in c("ARG_rel", "bacteria_rel", "fungi_rel", "KO_rel", "VF_rel","MetagTaxa_rel")){ 
    # df.name = "ARG_rel"
    
    mf_df <- eval(parse(text = df.name))
    
    # remove total microbial features
    totalMFeatures <- c(
        colnames(mf_df)[grepl("^bact",colnames(mf_df) )] %>% as.character() %>%  unique(),
        colnames(mf_df)[grepl("^fun",colnames(mf_df) )] %>% as.character() %>%  unique(),
        colnames(mf_df)[grepl("^metag",colnames(mf_df) )] %>% as.character() %>%  unique(),
        colnames(mf_df)[grepl("^Meta",colnames(mf_df) )]  %>% as.character() %>%  unique(),
        colnames(mf_df)[grepl("^ARG\\.",colnames(mf_df) )]  %>% as.character() %>%  unique(),
        colnames(mf_df)[grepl("^VF\\.",colnames(mf_df) )]  %>% as.character() %>%  unique()
    )
    
    mf_df <- mf_df %>% select(-all_of(totalMFeatures))
    
    
    for(v in vars2analyze){
        #v=vars2analyze[1]
        
        dat.test <- meta %>% select(all_of(c("SampleID", v)))
        dat.test <- dat.test[complete.cases(dat.test),]
        
        # match sps 
        common.sps <- intersect(rownames(mf_df), dat.test$SampleID)
        
        mf_df.tmp <- mf_df[match(common.sps, rownames(mf_df)),]
        dat.test <- dat.test[match(common.sps, dat.test$SampleID),]
        
        fml <- paste0("mf_df.tmp ~ ", v)
        Adonis <- try(adonis2(as.formula(fml), data = dat.test, permutations = 999, 
                              method="bray", by="margin"))
        
        if('try-error' %in% class(Adonis)){
            R2 <- NULL
            p <- NULL
        }else{
            R2 <- Adonis$R2[1]
            p <- Adonis$`Pr(>F)`[1]
        } 
        
        res <- c("data" = df.name,
                 "variable" = v,
                 "adonis.r2" = R2,
                 "adonis.p" = p)
        
        Adonis_res <- bind_rows(Adonis_res, res)
        
    }# loop through vars2analyze
}# loop through df.name



