library(ggplot2)
library(ggpubr)



### Boxplot
work_path <- '/home/rstudio/R/TcExh/figure/Figure6_new/'
### IFN-g
df <- data.frame(data1 = c(10.4, 10.3, 11.5, 10.3, 13.4, 50.7, 53.2, 52.1, 54.4, 55.7),
                 data2 = c('Vehicle','Vehicle','Vehicle','Vehicle','Vehicle','CADD522','CADD522','CADD522','CADD522','CADD522'))

### TNF-a
df <- data.frame(data1 = c(14.7,19.4,15.5,15.3,14.2,67.5,67.8,67.8,65.8,66.4),
                 data2 = c('Vehicle','Vehicle','Vehicle','Vehicle','Vehicle','CADD522','CADD522','CADD522','CADD522','CADD522'))

### GrB
df <- data.frame(data1 = c(21.4,23.5,18.5,12.5,11,81.9,79.4,79.3,78.2,79.1),
                 data2 = c('Vehicle','Vehicle','Vehicle','Vehicle','Vehicle','CADD522','CADD522','CADD522','CADD522','CADD522'))

### IL2
df <- data.frame(data1 = c(11.3,9.95,10.3,10.8,13.1,34.6,38.9,38,42.3,37.9),
                 data2 = c('Vehicle','Vehicle','Vehicle','Vehicle','Vehicle','CADD522','CADD522','CADD522','CADD522','CADD522'))

### PD1(%)
df <- data.frame(data1 = c(96.6,96.3,99.6,91,88.9,81.1),
                 data2 = c('Vehicle','Vehicle','Vehicle','CADD522','CADD522','CADD522'))

### PD1 MFI
df <- data.frame(data1 = c(75031,77385,76331,46579,47768,45499),
                 data2 = c('Vehicle','Vehicle','Vehicle','CADD522','CADD522','CADD522'))

### Cytotoxicity
df <- data.frame(data1 = c(10.55,10.65,10.85,9.82,10.85,10.85,41.95,42.85,39.85,40,42.35,40.55),
                 data2 = c('Vehicle','Vehicle','Vehicle','Vehicle','Vehicle','Vehicle','CADD522','CADD522','CADD522','CADD522','CADD522','CADD522'))

### Visualization
colnames(df) <- c('data','group')

p <- ggplot(df, aes(x=factor(group, levels = c('Vehicle','CADD522')), y=data)) + 
  geom_bar(aes(fill = group), fill = c('#B0E0D8', '#F7ABA5'), stat='summary', fun='mean', colour="black", size = 1.3, width = .7) +
  stat_summary(fun.data = mean_se,  geom = "errorbar", width = .15, colour = 'black', size = 1) + 
  geom_point(aes(fill = group), pch = 21, stroke = 1.5,
             position=position_dodge2(width=0.6), size = 3.5) + 
  scale_fill_manual(values = c('#F3877F', '#8ED3C7')) +
  geom_signif(comparisons = list(c('Vehicle','CADD522')), 
              annotations = "*", textsize = 15, size = 1, vjust = 0.5, 
              tip_length = c(0.95, 0.05), y_position = 79500) +
  theme_classic() + ylab(expression(atop("MFI of PD-1"^"+", "(% of the CD8"^"+"~"T cells)"))) + xlab("") +
  theme(text = element_text(family = 'Arial'),
        axis.line = element_line(size = 1.3),
        axis.ticks.length = unit(0.3,'cm'),
        axis.ticks = element_line(size = 1.3),
        axis.text.x = element_text(size = 22, color = 'black', vjust = 0), 
        axis.text.y = element_text(size = 22, color = 'black'),
        axis.title.y = element_text(size = 24, vjust = 1.3),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 16), legend.text = element_text(size = 15),
        legend.position = "None") +
  scale_y_continuous(breaks = c(40000,50000,60000,70000,80000), expand = c(0,0)) +
  coord_cartesian(ylim=c(40000,84000))
p
### save
png(
  filename = paste0(work_path, "Fig_6c_2.png"),
  width = 5,
  height = 4.5,
  units = 'in',
  res = 300
)

p

dev.off()
