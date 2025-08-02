###############
### Fig. 1b ###
###############

library(ggplot2)
library(ggpubr)

### read file ###
df <- read.csv('~/TcExh_source_data.csv')
df <- df[3:11,1:11]
colnames(df) <- df[1,]
df <- df[2:9,]
rownames(df) <- df[,1]
df <- df[,2:11]

df_stack <- stack(df)
df_stack$values <- as.numeric(df_stack$values)
df_stack$group <- rep(c('Isotype Ab','anti-PD1'), each = 40)
df_stack$day <- rep(c(7,14,21,28,35,42,49,56), 10)
### scatter plot ###

png('~/Desktop/fig1b.png',
    width = 15,
    height = 10,
    units = 'cm',
    res = 300)
ggplot(df_stack, aes(x = day, y = values)) +
  geom_line(aes(col = group, group = interaction(group, ind)), size = 1, linetype = 5) +
  geom_point(aes(col = group), size = 3) + 
  xlab('Days after implantation') + ylab(expression("Tumor volume"~"(mm"^"3"*")")) +
  scale_x_continuous(breaks = c(7,14,21,28,35,42,49,56)) +
  scale_color_manual(values = c("#E58C83", "#7F7F7F")) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 14, colour = 'black'),
        axis.text.y = element_text(size = 14, colour = 'black'),
        axis.title.x = element_text(size = 18, colour = 'black', vjust = 0),
        axis.title.y = element_text(size = 18, colour = 'black'),
        legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 12, colour = 'black')) +
  scale_y_continuous(limits = c(0,2500))
dev.off()
