TOLABEL <- SSc_DHSResults$Category[c(58,101,131,98,73,86,85,113,116)]

design <- "
  AAAABCDEEFFO
  GGHIJJJKLMNN
"

#PLOTS
P1 <- SSc_DHSResults%>% mutate(DescriptionAbr = fct_reorder(DescriptionAbr, Organ))%>%ggplot(aes(x=DescriptionAbr, 
                                                                                                 y=-log10(FDR), 
                                                                                                 size=RelEnrichment, 
                                                                                                 color=Organ,
                                                                                                 alpha=factor(FDR<0.05, levels=c(TRUE, FALSE), labels=c('p-value<0.05', 'NS'))))+
  geom_point()+
  scale_alpha_discrete(range=c(0.90, 0.30), name='')+
  scale_colour_manual('Organ',
                      values = c('#d62828', #Blood
                                 '#c8b6ff', #Brain
                                 '#fbff12', #Germ
                                 '#a8dadc', #Gum
                                 '#621708', #Heart
                                 '#a7c957', #Kidney
                                 '#4cc9f0', #Lung
                                 '#edafb8', #Mammary
                                 '#b5e48c', #Mesoderm
                                 '#e36414', #Muscle
                                 '#8338ec', #Ovary
                                 '#386641', #Pancreas
                                 '#1d3557', #Prostate
                                 '#ff9f1c', #Skin
                                 '#83c5be'))+ #Stroma
  geom_hline(yintercept = -log10(0.05), linewidth=0.5, linetype='dashed', color="#CF3245")+
  facet_manual(vars(Organ), design = design)+
  #force_panelsizes(rows = c(2,2,2), cols = c(10.5,5,1.5,1.5,2))+ 
  geom_label_repel(data = SSc_DHSResults[SSc_DHSResults$Category%in%TOLABEL,], aes(x=DescriptionAbr,
                                                                                   y=-log10(FDR),
                                                                                   label=str_wrap(AbrDESCRIPTOR, 10)),
                   size=4, force = 100, color='black')+
  ylab('-log10(FDR)')+
  theme_linedraw()+
  guides(color=guide_legend(override.aes= list(size=5)), alpha='none')+
  labs(size='Enrichment')+
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        axis.title.y=element_text(size=12, face = 'bold'),
        axis.title.x =element_blank(),
        legend.title = element_text(size=13, face = 'bold'),
        legend.text = element_text(size=12),
        strip.text = element_text(size = 11, face = 'bold', color="white"),
        strip.background = element_rect(fill = "#333333"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'right')
P1