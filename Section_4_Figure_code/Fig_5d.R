########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 5d ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### The ratio of RUNX2+ exhausted T cell / CD8+ T cell in specific carcinoma region were then calculated
### 46.34 (P1T), 35 (P3T), 25.16 (P5T), 67.36 (P8T), 77.04 (P11T) [Non-responders
### 23.46 (P7T), 15 (P9T), 23.1 (P10T) [Responders]

### Boxplot
ST.df <- data.frame(data1 = c(46.34, 35, 25.16, 67.36, 77.04, 23.46, 15, 23.1),
                    data2 = c('Non-Response','Non-Response','Non-Response','Non-Response','Non-Response','Response','Response','Response'))

colnames(ST.df) <- c('data','group')

p <- ggplot(ST.df, aes(x=factor(group, levels = c('Response','Non-Response')), y=data)) + 
  geom_boxplot(fill = c('#A8D0C7','#D78D82'), outlier.shape = NA) + geom_point(position=position_jitter(width=0.03)) + 
  geom_signif(comparisons = list(c("Response", "Non-Response")), test = t.test, test.args = c('less'), textsize = 6) + 
  theme_bw() + ylab('Proportion of RUNX2+\nExhausted T cells (%)') + xlab('HCC Patient') + ylim(10, 85) +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 16), legend.text = element_text(size = 15))

### save
png(
  filename = paste0(work_path, "Fig_5d.png"),
  width = 5,
  height = 4.5,
  units = 'in',
  res = 300
)
  
p
  
dev.off()
