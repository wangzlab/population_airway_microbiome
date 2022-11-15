library(data.table)


plotDat <- fread("fungi_PCA_Aspergillus.txt", data.table = F)
group_df <- fread("occupation_PC1_meta.txt", data.table = F)
plotDat.ellipse <- merge(plotDat, group_df, by.x = "V1", by.y = "#NAME")


library(ggplot2)
library(RColorBrewer )
# pick a nice color palette
palette(brewer.pal(n = 8, name = "Set2"))
myPalette <- c("#F1DE33", "#B5D236","#90C549","#56B664","#2BA67B","#229683","#268388","#2A7189",
               "#325D87","#3A4281","#442A75","#3B2955")

p0<- ggplot() +
  geom_point(data = plotDat,aes(x=PC1, y=PC2, color=Aspergillus), size=2.8, alpha=0.7) +
  stat_ellipse(data = plotDat.ellipse %>% filter(Occupation == "Y"),
               aes(x=PC1, y=PC2, group = Group3, linetype = Group),
               color="#A81F24", level = 0.8, lwd = 1) +
  stat_ellipse(data = plotDat.ellipse %>% filter(Occupation == "N"),
               aes(x=PC1, y=PC2, group = Group3, linetype = Group),
               color="darkgray",level = 0.8, lwd = 1) +
  scale_color_gradientn(colours = myPalette) +
  theme_bw() + theme(panel.grid = element_blank(),
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank())


p1 <- ggplot(data = plotDat.ellipse) +
  geom_density(aes(x = PC1, group = Group3, linetype = Group, color=Occupation),lwd = 1) +
  theme_bw() + theme(panel.grid = element_blank(), 
                     legend.position = "none",
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank()) +
  scale_color_manual(values = c("darkgray","#A81F24")) 
 
p2 <- ggplot(data = plotDat.ellipse) +
  geom_density(aes(x = PC2, group = Group3, linetype = Group, color=Occupation),lwd = 1) +
  theme_bw() + theme(panel.grid = element_blank(), 
                     legend.position = "none",
                     axis.title = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank()) +
  scale_color_manual(values = c("darkgray","#A81F24")) +
  coord_flip()


library(ggpubr)

ggarrange(
  ggarrange(p1, ggplot(), widths = c(0.8,0.2)),
  ggarrange(p0 + theme(legend.position =  c(0.85, 0.5)), p2, widths = c(0.8,0.2)),
  ncol = 1,  heights = c(0.2,0.8)
)
