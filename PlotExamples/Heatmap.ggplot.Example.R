DIS_RESULTS %>% mutate(Description = fct_reorder(Description, SSc))%>% 
  ggplot(aes(x=1, y=Description, fill=-log10(FDR)))+
  geom_tile(color='black')+
  facet_grid(factor(Organ, levels = c("Blood","Skin", "Lung", "Muscle", "Kidney", "Heart", "Gum", "Brain", "Mammary", "Mesoderm"))~., scales = 'free_y', space='free')+ #geom_tile does not work without scales=free
  force_panelsizes(rows = c(24.5,2.5,3.5,3,3,3,2,2.5,4,4.3))+ #This is to make sure that the title of each facet is embebbed within the panel
  xlab('SS')+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x =element_text(size = 16, face = 'bold'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=16),
        strip.text = element_text(size = 12, face = 'bold', color="white"),
        strip.background = element_rect(fill = "#333333"),
        axis.line = element_blank(),
        axis.ticks = element_blank())+
  guides(fill=guide_legend(title = "", title.position = "top"))
