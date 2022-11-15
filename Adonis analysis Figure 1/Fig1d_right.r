
# Right panel:  ========================================================
# VPA to parse overlapped contribution between factor and district
library(data.table)
library(dplyr)

sigVars <- fread("DAF.txt", data.table = F)$Variable
load("../_data/meta_phenotypes.RData")
meta<- merge(merge(meta.cat, meta.num, by="SampleID"),meta.geo_name, by="SampleID")
meta.sigVars <- meta %>% select(SampleID, all_of(sigVars), District) 
meta.sigVars <- meta.sigVars[complete.cases(meta.sigVars),]

colnames(meta.sigVars) <- gsub("\\W","_",colnames(meta.sigVars))
varType_df$Variable <- sub("\\W", "_", varType_df$Variable)


# load data  -------------------------
# bacteria
dat <- fread("../_data/bacteria_L6.txt")
bact.betaDivers <-  dat %>%
  tibble::column_to_rownames("#NAME") %>%
  t() %>% as.data.frame()
bact.shannon <- apply(bact.betaDivers, 1, diversity) %>% as.data.frame()
colnames(bact.shannon) <- "Shannon"

# fungi
dat <- fread("../_data/fungi_L6.txt")
fung.betaDivers <- dat %>% tibble::column_to_rownames("#NAME") %>%
  t() %>% as.data.frame()
fung.shannon <- apply(fung.betaDivers, 1, diversity) %>% as.data.frame()
colnames(fung.shannon) <- "Shannon"

# KO
KO.dat.full <- fread("../_data/metag_all_ko.txt", data.table = F)  %>% tibble::column_to_rownames("SampleID")
tmp <- KO.dat.full %>%
  apply(1, mean) %>% as.data.frame() 
colnames(tmp) <- "avg"
tmp <- tmp %>% arrange(desc(avg))

KO.betaDivers <- KO.dat.full %>% tibble::rownames_to_column("KO") %>%
  filter(KO %in% rownames(tmp)[1:500]) %>%  
  tibble::column_to_rownames("KO") %>%
  t() %>% as.data.frame()

KO.shannon <- apply(KO.betaDivers, 1, diversity) %>% as.data.frame()
colnames(KO.shannon) <- "Shannon"


# ARG
ARG.betaDivers <- fread("../_data/ARG.txt") %>% 
  tibble::column_to_rownames("SampleID") %>%
  t() %>% as.data.frame()
ARG.shannon <- apply(ARG.betaDivers, 1, diversity) %>% as.data.frame()
colnames(ARG.shannon) <- "Shannon"
ARG.abund <- apply(ARG.betaDivers, 1, sum)  %>% as.data.frame()

# VF
VF.betaDivers <- fread("../_data/VF.txt") %>%
  tibble::column_to_rownames("SampleID") %>%
  t() %>% as.data.frame()
VF.shannon <- apply(VF.betaDivers, 1, diversity) %>% as.data.frame()
colnames(VF.shannon) <- "Shannon"
VF.abund <- apply(VF.betaDivers, 1, sum)  %>% as.data.frame()


# auto optimization for 1-dimension response variables ------------------
VPA.res <- NULL
for(df.name in c( "bact.shannon","fung.shannon", "KO.shannon","ARG.shannon",
                  "ARG.abund","VF.shannon","VF.abund")){
  # df.name = "bact.shannon"
  dat.t <- eval(parse(text = df.name))
  common.sp <- intersect(rownames(dat.t),meta.sigVars$SampleID)
  dat.t <- dat.t[match(common.sp, rownames(dat.t)),]
  explainVars.w.orig <- meta.sigVars[match(common.sp, meta.sigVars$SampleID), ]
  rownames(explainVars.w.orig) <- NULL
  explainVars.w.orig <- explainVars.w.orig %>% tibble::column_to_rownames("SampleID")
  
  # optimization ----
  var.rm=NULL
  otu_cca_adj.0 = 0
  otu_cca_noadj.0 = 0
  while(T){
    
    if(!is.null(var.rm)) explainVars.w <- explainVars.w.orig %>% select(-all_of(var.rm)) else   explainVars.w <- explainVars.w.orig
    m.rda <- rda(dat.t~., explainVars.w)
    r2 <- RsquareAdj(m.rda)
    otu_cca_noadj <- r2$r.squared; otu_cca_noadj	
    otu_cca_adj <- r2$adj.r.squared; otu_cca_adj	
    
    if(otu_cca_adj < otu_cca_adj.0) {
      writeLines(paste("The optimized original r2 is: ", round(otu_cca_noadj.0,3) , "; r2adj is: ", round(otu_cca_adj.0,3), sep = ""))
      break
    } else {
      otu_cca_adj.0 = otu_cca_adj
      otu_cca_noadj.0 = otu_cca_noadj
    }
    
    otu_cca_test_term <- anova.cca(m.rda, by = 'terms', permutations = 999)
    otu_cca_test_term$padj <- p.adjust(otu_cca_test_term$`Pr(>F)`, method = 'bonferroni') #p 值校正（Bonferroni 为例）
    otu_cca_test_term <- otu_cca_test_term %>% arrange(`Pr(>F)`)
    
    var.rm <- append(var.rm, rownames(otu_cca_test_term)[nrow(otu_cca_test_term)-1])
  }
  
  # DAF Vs District contribution to effect -----
  var.rm.0 <- var.rm[1:(length(var.rm) - 1)]
  vpa <- varpart(dat.t, 
                 explainVars.w.orig %>% select(District),
                 explainVars.w.orig %>% select(-all_of(unique(c(var.rm.0,"District")))))
  
  # contribution of the overall DAFs
  vpa_part<-vpa$part
  
  t1<-vpa_part$fract
  t2 <- vpa_part$indfract
  t3<-vpa_part$contr1
  t4<-vpa_part$contr2
  
  res_c <- c("y" = df.name, "district" = t2$Adj.R.squared[1], 
             "overlapped"= t2$Adj.R.squared[2], "DAF"= t2$Adj.R.squared[3], 
             "DAF.category" = "all")

  VPA.res <- bind_rows(VPA.res, res_c)
  
  # contribution of DAFs in each variable category -----
  var.remain <- colnames(explainVars.w.orig)[!(colnames(explainVars.w.orig) %in% c(var.rm.0, "District"))]
  vartypes <- sapply(var.remain, function(x) varType_df$VarType[which(varType_df$Variable == x)])
  for(vt in unique(vartypes)){
    vpa <- varpart(dat.t, 
                   explainVars.w.orig %>% select(District),
                   explainVars.w.orig %>% select(names(which(vartypes == vt))))
    vpa_part<-vpa$part
    t2 <- vpa_part$indfract
    res_c <- c("y" = df.name, "district" = t2$Adj.R.squared[1], 
               "overlapped"= t2$Adj.R.squared[2], "DAF"= t2$Adj.R.squared[3], 
               "DAF.category" = vt)
    VPA.res <- bind_rows(VPA.res, res_c)
  }
}


VPA.res$overlapped <- as.numeric(VPA.res$overlapped)
VPA.res$district <- as.numeric(VPA.res$district)
VPA.res$DAF <- as.numeric(VPA.res$DAF)
VPA.res$DAF.proportion <- pmax(VPA.res$overlapped, 0)/(pmax(VPA.res$overlapped, 0) + VPA.res$district)

VPA.res.1 <- VPA.res


# optimization for beta diversities ------------------
# first identify which method to use: rda or cca
dat.t <- bact.betaDivers
common.sp <- intersect(rownames(dat.t),meta.sigVars$SampleID)
dat.t <- dat.t[match(common.sp, rownames(dat.t)),]
explainVars.w.orig <- meta.sigVars[match(common.sp, meta.sigVars$SampleID), ]

DCA=decorana(dat.t)
DCA # manual check length of first axis

# after manual check, below are DCA axis1 length for all microbial data ----
DCA_axis <- bind_rows(c("df" = "bact.betaDivers", "axis1" = 1.582),
                      c("df" = "fung.betaDivers", "axis1" = 4.1567),
                      c("df" = "ARG.betaDivers", "axis1" = 4.4213),
                      c("df" = "VF.betaDivers", "axis1" = 4.7823),
                      c("df" = "KO.betaDivers", "axis1" = 2.4251))

# auto optimization ---------
optimizedRes <- NULL
for(df.name in c( "bact.betaDivers","ARG.betaDivers","fung.betaDivers","KO.betaDivers","VF.betaDivers")){
  # df.name = "bact.betaDivers"
  dat.t <- eval(parse(text = df.name))
  common.sp <- intersect(rownames(dat.t),meta.sigVars$SampleID)
  dat.t <- dat.t[match(common.sp, rownames(dat.t)),]
  explainVars.w.orig <- meta.sigVars[match(common.sp, meta.sigVars$SampleID), ]
  rownames(explainVars.w.orig) <- NULL
  explainVars.w.orig <- explainVars.w.orig %>% tibble::column_to_rownames("SampleID")
  
  axis1 = DCA_axis$axis1[which(DCA_axis$df == df.name)] 
  if(axis1 > 4) md = "cca" else md = "rda"
  
  if(md == "rda"){
    var.rm=NULL
    otu_cca_adj.0 = 0
    otu_cca_noadj.0 = 0
    while(T){
      
      if(!is.null(var.rm)) explainVars.w <- explainVars.w.orig %>% select(-all_of(var.rm)) else   explainVars.w <- explainVars.w.orig
      m.rda <- rda(dat.t~., explainVars.w)
      r2 <- RsquareAdj(m.rda)
      otu_cca_noadj <- r2$r.squared; otu_cca_noadj	
      otu_cca_adj <- r2$adj.r.squared; otu_cca_adj	
      
      if(otu_cca_adj < otu_cca_adj.0) {
        writeLines(paste("The optimized original r2 is: ", round(otu_cca_noadj.0,3) , "; r2adj is: ", round(otu_cca_adj.0,3), sep = ""))
        break
      } else {
        otu_cca_adj.0 = otu_cca_adj
        otu_cca_noadj.0 = otu_cca_noadj
      }
      
      otu_cca_test_term <- anova.cca(m.rda, by = 'terms', permutations = 999)
      otu_cca_test_term$padj <- p.adjust(otu_cca_test_term$`Pr(>F)`, method = 'bonferroni') #p 值校正（Bonferroni 为例）
      otu_cca_test_term <- otu_cca_test_term %>% arrange(`Pr(>F)`)
      
      var.rm <- append(var.rm, rownames(otu_cca_test_term)[nrow(otu_cca_test_term)-1])
    }
  }else{
    var.rm = NULL
    otu_cca_adj.0 = 0
    otu_cca_noadj.0 = 0
    while(T){
      
      if(!is.null(var.rm)) explainVars.w <- explainVars.w.orig %>% select(-all_of(var.rm)) else   explainVars.w <- explainVars.w.orig
      m.cca<- cca(dat.t~., explainVars.w)
      r2 <- RsquareAdj(m.cca)
      otu_cca_noadj <- r2$r.squared; otu_cca_noadj	
      otu_cca_adj <- r2$adj.r.squared; otu_cca_adj	
      
      if(otu_cca_adj <= otu_cca_adj.0) {
        writeLines(paste("The optimized original r2 is: ", round(otu_cca_noadj.0,3) , "; r2adj is: ",  round(otu_cca_adj.0,3), sep = ""))
        break
      } else {
        otu_cca_adj.0 = otu_cca_adj
        otu_cca_noadj.0 = otu_cca_noadj
      }
      
      otu_cca_test_term <- anova.cca(m.cca, by = 'terms', permutations = 999)
      otu_cca_test_term$padj <- p.adjust(otu_cca_test_term$`Pr(>F)`, method = 'bonferroni') #p 值校正（Bonferroni 为例）
      otu_cca_test_term <- otu_cca_test_term %>% arrange(`Pr(>F)`)
      
      var.rm <- append(var.rm, rownames(otu_cca_test_term)[nrow(otu_cca_test_term)-1])
    }
  }
  
  var.rm
  vpa <- varpart(dat.t, 
                 explainVars.w.orig %>% select(District),
                 explainVars.w.orig %>% select(-all_of(unique(c(var.rm,"District")))))
  
  vpa_part<-vpa$part
  
  t1<-vpa_part$fract
  t2 <- vpa_part$indfract
  t3<-vpa_part$contr1
  t4<-vpa_part$contr2
  
  res_c <- c("y" = df.name, "model" = md, "vars.rm" = paste(var.rm, collapse = ";"),
             "district" = t2$Adj.R.squared[1],  "overlapped"= t2$Adj.R.squared[2], "DAF"= t2$Adj.R.squared[3])
  optimizedRes <- bind_rows(optimizedRes, res_c)
  
}

# calculate contribution based on the optimized models ---------

VPA.res<- NULL
for(df.name in c("bact.betaDivers", "ARG.betaDivers",  "fung.betaDivers", "KO.betaDivers")){ # "VF.betaDivers"
  #print(df.name)
  
  dat.t = eval(parse(text = df.name))
  common.sp <- intersect(rownames(dat.t),meta.sigVars$SampleID)
  dat.t <- dat.t[match(common.sp, rownames(dat.t)),]
  explainVars.w.orig <- meta.sigVars[match(common.sp, meta.sigVars$SampleID), ]
  rownames(explainVars.w.orig) <- NULL
  explainVars.w.orig <- explainVars.w.orig %>% tibble::column_to_rownames("SampleID")
  
  
  md = optimizedRes$model[which(optimizedRes$y == df.name)]
  var.rm = strsplit(optimizedRes$vars.rm[which(optimizedRes$y == df.name)], ";", fixed = T)[[1]] 
  var.rm.0 <- var.rm[1:(length(var.rm) - 1)]
  # 计算总的：
  if(md == "rda"){
    vpa <- varpart(dat.t, 
                   explainVars.w.orig %>% select(District),
                   explainVars.w.orig %>% select(-all_of(unique(c(var.rm.0,"District")))))
  }else{
    vpa <- varpart(dat.t, 
                   explainVars.w.orig %>% select(District),
                   explainVars.w.orig %>% select(-all_of(unique(c(var.rm.0,"District")))),
                   chisquare=T)
  }
  vpa_part<-vpa$part
  t2 <- vpa_part$indfract
  res_c <- c("y" = df.name, "district" = t2$Adj.R.squared[1],  "overlapped"= t2$Adj.R.squared[2], "DAF"= t2$Adj.R.squared[3], 
             "DAF.category" = "all")
  VPA.res <- bind_rows(VPA.res, res_c)
  
  
  # 分开计算各个categories
  var.remain <- colnames(explainVars.w.orig)[!(colnames(explainVars.w.orig) %in% c(var.rm.0, "District"))]
  vartypes <- sapply(var.remain, function(x) varType_df$VarType[which(varType_df$Variable == x)])
  for(vt in unique(vartypes)){
    if(md == "rda"){
      vpa <- varpart(dat.t, 
                     explainVars.w.orig %>% select(District),
                     explainVars.w.orig %>% select(names(which(vartypes == vt))))
    }else{
      vpa <- varpart(dat.t, 
                     explainVars.w.orig %>% select(District),
                     explainVars.w.orig %>% select(names(which(vartypes == vt))),
                     chisquare=T)
    }
    vpa_part<-vpa$part
    t2 <- vpa_part$indfract
    res_c <- c("y" = df.name, "district" = t2$Adj.R.squared[1],
               "overlapped"= t2$Adj.R.squared[2], "DAF"= t2$Adj.R.squared[3], 
               "DAF.category" = vt)
    VPA.res <- bind_rows(VPA.res, res_c)
  }
}
VPA.res$overlapped <- as.numeric(VPA.res$overlapped)
VPA.res$district <- as.numeric(VPA.res$district)
VPA.res$DAF <- as.numeric(VPA.res$DAF)
VPA.res$DAF.proportion <- pmax(VPA.res$overlapped, 0)/(pmax(VPA.res$overlapped, 0) + VPA.res$district)
VPA.res.2 <- VPA.res


VPA.res <- bind_rows(VPA.res.1,VPA.res.2)
write.table(VPA.res, file = "VPA_DAF.vs.district.txt",  sep = "\t", quote = F, row.names = F)
