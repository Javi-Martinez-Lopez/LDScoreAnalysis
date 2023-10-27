library(tidyverse)
library(geomtextpath)

# Create dataset
data <- data.frame(
  individual=SScResultsDHS$Description,
  group=SScResultsDHS$Organ,
  value=SScResultsDHS$RelEnrichment,
  alpha=-log10(SScResultsDHS$FDR)
)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(data$group), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(group)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_hline(yintercept = 0, linewidth=0.3, linetype='dashed', color='#444444')+
  geom_hline(yintercept = 0.25, linewidth=0.3, linetype='dashed', color='#444444')+
  geom_hline(yintercept = 0.5, linewidth=0.3, linetype='dashed', color='#444444')+
  geom_hline(yintercept = 0.75, linewidth=0.3, linetype='dashed', color='#444444')+
  geom_hline(yintercept = 1, linewidth=0.3, linetype='dashed', color='#444444')+
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity") +
  geom_vline(xintercept = 0, linewidth=0.5)+
  geom_texthline(label = "Enrichment", yintercept=1.01, size = 5, vjust = 1.3,
                linewidth = 1, linecolor = "white", linetype = 2, 
                color = "#000000")+
  annotate('text', x = c(nrow(SScResultsDHS)-1.8, 
                         nrow(SScResultsDHS)-1.5, 
                         nrow(SScResultsDHS)-1.2, 
                         nrow(SScResultsDHS)-0.5), 
           y = c(0.3, 0.55, 0.8, 1.05), 
           label = c('0.25', '0.50', '0.75', '1'), 
           size=3.5, 
           angle=10)+
  scale_fill_manual('Organ',
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
  ylim(-0.3,1.1) +
  #scale_alpha_continuous('-log10(FDR)', range = c(0.3,1))+
  theme_void() +
  theme(
    legend.position="right", 
    legend.box = "vertical",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = 'bold'),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.box.background = element_rect(fill = '#FFFFFF', colour = 'white', linewidth = 0.1)
  ) +
  coord_curvedpolar()
  
p
  
ggsave("Graph6.jpeg", p, width = 30, height = 20, units = 'cm', dpi = 300)


