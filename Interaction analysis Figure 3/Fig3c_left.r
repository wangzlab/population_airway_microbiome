# Adonis to calculate effect of occupation exposure on bacteria, fungi, ko, ARG, VF 
# in Nesseria high and Nesseria low groups respectively 
# adonis(microbialDataFrame ~ exposure, by="margin")  



library(data.table)
library(dplyr)
library(vegan)

grouping_df <- fread("occupation_PC1_meta.txt", data.table=F)


# load microbial feature data ------------------------
load("../_data/microbialFeatures_relAbund_rpkm_notStandardized.RData")


# exposures --------------------------------------
load("../_data/meta_phenotypes.RData")
meta <-merge(
  merge(meta.num, meta.cat,  by = "SampleID"),
  meta.geo_name,by = "SampleID"
) 


Expo.vars <- c("Occupational_pollution")

meta.simp <- meta %>% select(all_of(c("SampleID", Expo.vars)) )
meta.simp <- meta.simp[complete.cases(meta.simp),]


# Adonis --------------------------------------



Adonis_res <- NULL

for(df.name in c(bacteria_rel,fungi_rel,KO_rel,ARG_rel,VF_rel)){
  # df.name =  "bacteria_rel"
  
  dat <- eval(parse(text = df.name))
  dat <- dat[complete.cases(dat),]
  
  totalMFeatures <- c(
    colnames(dat)[grepl("^bact",colnames(dat))] %>% as.character() %>%  unique(),
    colnames(dat)[grepl("^fung",colnames(dat))] %>% as.character() %>%  unique(),
    colnames(dat)[grepl("^KO",colnames(dat))]  %>% as.character() %>%  unique(),
    colnames(dat)[grepl("^ARG\\.",colnames(dat))]  %>% as.character() %>%  unique(),
    colnames(dat)[grepl("^VF\\.",colnames(dat))]  %>% as.character() %>%  unique()
  )
  
  dat <- dat %>% select(-all_of(totalMFeatures))
  
  # match microbial feature with meta 
  commonSps <- intersect(intersect(rownames(dat), meta.simp$SampleID), grouping_df$`#NAME` )
  
  dat <- dat[match(commonSps, rownames(dat)),]
  vars_df <- meta.simp[match(commonSps, meta.simp$SampleID),]
  grouping_df.match <- grouping_df[match(commonSps, grouping_df$`#NAME`),]
  
  for(expo in Expo.vars){
    # expo = Expo.vars[1]
    
    for(grp in c("Neis_high","Neis_low")){
      # grp = "Neis_high"
      sps <- grouping_df.match$`#NAME`[which(grouping_df.match$Group == grp)]
      dat.tmp <- dat[sps,]
      vars_df.tmp <- vars_df[which(vars_df$SampleID %in% sps),]
      grouping_df.tmp <- grouping_df.match[which(grouping_df.match$`#NAME` %in% sps),]
      
      fml <- paste0("dat.tmp ~ ", expo)
      Adonis <- try(adonis2(as.formula(fml), data = vars_df.tmp, permutations = 999, method="bray", by="margin"))
      
      if('try-error' %in% class(Adonis)){
        R2 <- NULL
        p <- NULL
      }else{
        R2 <- Adonis$R2[1]
        p <- Adonis$`Pr(>F)`[1]
      } 
      
      res <- c("df" = df.name,
               "exposure" = expo,
               "group" = grp,
               "adonis.r2" = R2,
               "adonis.p" = p)
      
      Adonis_res <- bind_rows(Adonis_res, res)
      
      
    }# loop through Nesseria type
    
  }# loop through exposure variables
  
}# loop through microbial feature data frames


write.table(Adonis_res, file = "Adonis_multiDimMicro_Occupation_Nesseria.h.l.txt", quote = F, row.names = F, sep = "\t")

