P4 <- DIS_RESULTS2 %>% mutate(NAME = fct_reorder(NAME, ORGAN))%>% 
  ggplot(aes(x=NAME, y=DIS,  fill=RELENRICHMENT))+
  geom_tile(color='black', show.legend = T)+
  facet_grid(.~ORGAN, scales = 'free', space = 'free')+
  force_panelsizes(cols = c(60,15,2,3.5,5.5,10,8,6,2,3.5,2,2,2.5,13,5))+  
  scale_fill_gradient2(low = 'white', high = '#AF2314', mid = '#FFFFFF', midpoint = 0.1)+
  theme_classic()+
  labs(fill='Mean TAU')+
  theme(axis.text.y = element_text(size=14, hjust = 0.5, face = 'bold', color = 'black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=11, angle = 90, vjust = 0.5),
        axis.title.y =element_blank(),
        legend.title = element_text(size=13, face = 'bold'),
        legend.text = element_text(size=12),
        strip.text = element_text(size = 0, face = 'bold', color="#FFFFFF"),
        strip.background = element_rect(fill = '#333333'),
        axis.line = element_blank(),
        axis.ticks = element_blank())


P34 <- ggplot_gtable(ggplot_build(P4))
strip_both <- which(grepl('strip-', P34$layout$name))
fills <- c('#d62828', #Blood
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
           '#83c5be')

cols <- c('#a62828', #Blood
          '#98a6ff', #Brain
          '#dbdf12', #Germ
          '#88aaac', #Gum
          '#321708', #Heart
          '#97a957', #Kidney
          '#4ac9c0', #Lung
          '#cd9fa8', #Mammary
          '#a5c48c', #Mesoderm
          '#c36414', #Muscle
          '#8338ac', #Ovary
          '#385641', #Pancreas
          '#1d3557', #Prostate
          '#cf9f1c', #Skin
          '#83a59e')
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', P34$grobs[[i]]$grobs[[1]]$childrenOrder))
  P34$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  P34$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- cols[k]
  k <- k+1
}

grid.draw(P34)