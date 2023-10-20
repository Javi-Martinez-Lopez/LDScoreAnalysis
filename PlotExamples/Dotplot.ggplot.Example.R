library(ggplot2)
library(patchwork)
library(ggh4x)
library(dplyr)

#We set a filter to include only significant annotations
P1 <- subset(SScResultsDHS, FDR<0.05) %>% mutate(Description = fct_reorder(Description, desc(FDR))) %>%
  ggplot(aes(x=Description, y=-log10(FDR), size=RelEnrichment, fill=Organ))+
  geom_point(shape=21, stroke=1.1, colour="#000000")+
  coord_flip()+
  geom_hline(yintercept = -log10(0.05), linewidth=0.8, linetype="dashed", color="#676767")+
  facet_grid(factor(Organ, levels = c("Blood","Skin", "Lung", "Muscle", "Kidney", "Heart", "Gum", "Brain", "Mammary", "Mesoderm"))~., scales = "free_y", space = "free")+ #Scales free allows to set correctly the y axis. 
  force_panelsizes(rows = c(24.5,2.5,3.5,3,3,3,2,2.5,4,4.3))+
  #scale_color_manual(values = c('#c1121f','#f28482','#8ecae6','#6a040f','#f77f00','#a7c957','#cdb4db','#e9c46a','#3e5c76','#ffb703'))+
  scale_fill_manual(values = c('#c1121f','#f28482','#8ecae6','#6a040f','#f77f00','#a7c957','#cdb4db','#e9c46a','#3e5c76','#ffb703'))+
  theme_classic()+
  guides(fill = guide_legend(override.aes = list(size = 10)), 
         size=guide_legend(title = "Enrichment"),
         alpha="none", 
         color="none")+
  theme(axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=14),
        axis.title.y=element_blank(),
        axis.title.x =element_text(size = 16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=16),
        strip.text = element_blank(),
        strip.background = element_blank())
