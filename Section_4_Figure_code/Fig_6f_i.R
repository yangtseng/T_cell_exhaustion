#################
### Fig. 6F-I ###
#################

### libraries
library(ggplot2)
library(ggpubr)
library(ggsignif)


### Fig. 6g
### two-group scatter plot with standard error
### load data
df <- data.frame("group" = rep(c('IgG control','anti-PD1'), 7),
                 "avg" = c(0, 0, 424.5278,413.982,
                           531.8384,515.3424,
                           636.5026,657.5537,
                           727.3257,808.3375,
                           1281.814,1204.6422,
                           1578.4719,1682.6798),
                 "stdev" = c(0,0,203.3134262,124.2567226,
                             180.0611985,117.6163392,
                             235.0799646,177.8063884,
                             351.1884238,288.9909875,
                             740.777847,475.9369418,
                             930.662574,841.8439754),
                 'day' = rep(c('0','12','18','24','30','36','42'), each = 2))

df$day <- as.numeric(df$day)
df$group <- as.factor(df$group)
levels(df$group) <- c("IgG control", 'anti-PD1')

png(filename = "Fig_6f.png",
    width = 5,
    height = 3,
    units = 'in',
    res = 300)
  
ggplot(df, aes(x=day, y=avg, group=group, color=group)) + 
geom_line(size = .6) +
geom_point(size = 2.5)+
geom_errorbar(aes(ymin=avg-stdev, ymax=avg+stdev), width=1, size = .8) +
theme_classic() + xlab('Days post implantation') + 
ylab(expression("Tumor volume (mm"^3~")")) + 
theme(axis.title.x = element_text(size = 18, color = 'black',vjust = 0),
      axis.title.y = element_text(size = 18, color = 'black',vjust = 2),
      axis.text = element_text(size = 16, color = 'black'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = 'black'),
      axis.ticks = element_line(color = 'black'),
      legend.title = element_text(size = 15), 
      legend.text = element_text(size = 12)) + 
scale_color_manual("Group", values = c('#A8D0C7','#2F5597')) + 
scale_x_continuous(limits = c(0,43), expand = c(0,0)) +
scale_y_continuous(limits = c(0,2600), expand = c(0,0), breaks = c(0,500,1000,1500,2000,2500))

dev.off()

### Fig. 6h
### two-group scatter plot with standard error
### load data
df <- data.frame("group" = rep(c('Vehicle','CADD522'), 7),
                 "avg" = c(0, 0, 264.9567,233.4502,
                           409.43555,305.726593,
                           572.1346,373.9791198,
                           738.8693,253.3689,
                           1563.7617,517.6469,
                           1882.350389,753.4451667),
                 "stdev" = c(0,0,102.1530225,81.93117693,
                             129.3324983,106.1825265,
                             196.7047355,205.7452978,
                             249.7365549,186.1247632,
                             476.237222,506.3902983,
                             284.6975272,670.6084131),
                 'day' = rep(c('0','8','16','23','30','37','47'), each = 2))

df$day <- as.numeric(df$day)
df$group <- as.factor(df$group)
levels(df$group) <- c("Vehicle", "CADD522")

png(filename = "Fig_6h.png",
    width = 5,
    height = 3,
    units = 'in',
    res = 300)

ggplot(df, aes(x=day, y=avg, group=group, color=group)) + 
  geom_line(size = .6) +
  geom_point(size = 2.5)+
  geom_errorbar(aes(ymin=avg-stdev, ymax=avg+stdev), width=1, size = .8) +
  theme_classic() + xlab('Days post implantation') + 
  ylab(expression("Tumor volume (mm"^3~")")) + 
  theme(axis.title.x = element_text(size = 18, color = 'black',vjust = 0),
        axis.title.y = element_text(size = 18, color = 'black',vjust = 2),
        axis.text = element_text(size = 16, color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 12)) + 
  scale_color_manual("Group", values = c('#D78D82', '#A8D0C7'),
                     guide = guide_legend(reverse = TRUE)) + 
  annotate("text",label = "*", x = 23, y = 830, size = 5) +
  annotate("text",label = "***", x = 30, y = 1050, size = 5) +
  annotate("text",label = "***", x = 37, y = 2100, size = 5) +
  annotate("text",label = "***", x = 47, y = 2210, size = 5) +
  scale_x_continuous(limits = c(0,49), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,2300), expand = c(0,0), breaks = c(0,500,1000,1500,2000))

dev.off()

### Fig. 6h
### two-group scatter plot with standard error
### load data
df <- data.frame("values" = c(100,90,80,70,50,40,20,10,0,100,90,80,70,60,60),
                 "ind" = c(rep('Vehicle', 9), rep('CADD522',6)),
                 "day" = c(0,35,36,42,44,49,50,54,56,0,46,48,53,65,70))
df$values <- df$values/100
df$group <- df$ind
levels(df$group) <- c("Vehicle", 'CADD522')
png(filename = "Fig_6i.png",
    width = 5,
    height = 3,
    units = 'in',
    res = 300)

ggplot(df, aes(x=day, y=values, group=group, color=group)) + 
  geom_step(size=1) +
  theme_classic() + xlab('Days post implantation') + 
  ylab(expression("Survival rate")) + 
  theme(axis.title.x = element_text(size = 17, color = 'black', vjust = 0),
        axis.title.y = element_text(size = 17, color = 'black', vjust = 2),
        axis.text = element_text(size = 16, color = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 12)) + 
  scale_color_manual("Group", values = c('#D78D82','#A8D0C7'),
                     guide = guide_legend(reverse = TRUE)) + 
  scale_x_continuous(limits = c(0,70), breaks = c(0,10,20,30,40,50,60,70), expand = expansion(add = c(0,1))) +
  scale_y_continuous(limits = c(0,1), expand = expansion(add = c(0,.05))) +
  annotate("text",label = "p-value = 0.0009", x = 25, y = .25, size = 5)

dev.off()

