library(data.table)
library(dplyr)
library(ggplot2)
library(vegan)

# bacteria ---------------
bact.dat <- fread("../_data/bacteria_hierarchial.txt")
bact.dat$Group <- factor(bact.dat$Group,
                    levels = c("Same_neighbor","Same_street","Same_district","Diff_district"))


Fig1c.bact <- ggplot() +
    geom_boxplot(data=bact.dat, aes(x= Group, y=Distance, fill=Group), 
                 width=0.8, outlier.shape = NA) +
    scale_fill_manual(values = c()) +  # the colors
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(label = 'Bacteria')

# fungi ---------------
fung.dat <- fread("../_data/fungi_hierarchial.txt")
fung.dat$Group <- factor(fung.dat$Group,
                         levels = c("Same_neighbor","Same_street","Same_district","Diff_district"))


Fig1c.fung <- ggplot() +
    geom_boxplot(data=fung.dat, aes(x= Group, y=Distance, fill=Group), 
                 width=0.8, outlier.shape = NA) +
   # scale_fill_manual(values = c()) +  # the colors
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(label = 'Fungi')



# kegg ko ------------------
load("../_data/microbialFeatures_relAbund_rpkm_notStandardized.RData")
KO_rel <- KO_rel %>% select(-KO.shannon)
BC.KO <- vegdist(KO_rel) %>%  as.matrix() 
#BC.KO[lower]
BC.KO.l <- BC.KO %>% reshape2::melt(value.name = "Distance")

KO.dat <- merge(bact.dat %>% select(Sample1,Sample2,Group),
                BC.KO.l, 
                by.x = c("Sample1","Sample2"), 
                by.y = c("Var1","Var2"))

Fig1c.KO <- ggplot() +
    geom_boxplot(data=KO.dat, aes(x= Group, y=Distance, fill=Group), 
                 width=0.8, outlier.shape = NA) +
    # scale_fill_manual(values = c()) +  # the colors
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(label = 'KO')


# ARG ----------
ARG_rel <- ARG_rel %>% select(-ARG.shannon, -ARG.totalRpkm)
BC.ARG <- vegdist(ARG_rel) %>%  as.matrix() 
#BC.ARG[lower]
BC.ARG.l <- BC.ARG %>% reshape2::melt(value.name = "Distance")

ARG.dat <- merge(bact.dat %>% select(Sample1,Sample2,Group),
                BC.ARG.l, 
                by.x = c("Sample1","Sample2"), 
                by.y = c("Var1","Var2"))

Fig1c.ARG <- ggplot() +
    geom_boxplot(data=ARG.dat, aes(x= Group, y=Distance, fill=Group), 
                 width=0.8, outlier.shape = NA) +
    # scale_fill_manual(values = c()) +  # the colors
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(label = 'ARG')


# VF ----------------
VF_rel <- VF_rel %>% select(-VF.shannon, -VF.totalRpkm)

BC.VF <- vegdist(VF_rel) %>%  as.matrix() 
#BC.VF[lower]
BC.VF.l <- BC.VF %>% reshape2::melt(value.name = "Distance")

VF.dat <- merge(bact.dat %>% select(Sample1,Sample2,Group),
                 BC.VF.l, 
                 by.x = c("Sample1","Sample2"), 
                 by.y = c("Var1","Var2"))

Fig1c.VF <- ggplot() +
    geom_boxplot(data=VF.dat, aes(x= Group, y=Distance, fill=Group), 
                 width=0.8, outlier.shape = NA) +
    # scale_fill_manual(values = c()) +  # the colors
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(label = 'VF')


# calculate wilcox p value 
wilcox.p_resDf <- NULL

for(df.n  in c("bact.dat", "fung.dat", "KO.dat", "ARG.dat", "VF.dat")){
    
    dat <- eval(parse(text = df.n))
    for(cmb in list(c("Same_neighbor", "Same_street" ), 
                    c("Same_neighbor", "Same_district"),
                    c("Same_neighbor" ,"Diff_district"),
                    c("Same_street",   "Same_district"),
                    c("Same_street" ,  "Diff_district"),
                    c( "Same_district", "Diff_district"))){
        res <- dat %>% 
            filter(Group %in% cmb) %>%
            mutate(nonsense =  "1") %>%
            group_by(nonsense) %>%
            do(w = wilcox.test(Distance~Group, data = .)) %>%
            summarise(wilcox.p =  w$p.value) %>%
            mutate(compare = paste(cmb, collapse = " vs. ")) %>%
            mutate(data = df.n)
        
        wilcox.p_resDf <- bind_rows(wilcox.p_resDf, res)
    }
    
}

# manually add stars to the figures
