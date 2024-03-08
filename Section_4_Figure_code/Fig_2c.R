########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 2c ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
tcell <- readRDS("./Cellline/Cellline_motif4.rds"))

### This is the code for marker gene visualization
### We created the FeaturePlot of several well-known T cell marker genes
tcell_feature <- c('Nkg7','Ctla4','Il7r','Prf1','Pdcd1','Lef1')

### Annotation
seg <- data.frame(x = 1, xend = 1, y = 0.75, yend = 1.25)
row1 <- ggplot(seg) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_void() + geom_richtext(aes(x = 1 ,y = 1, label = "Effector markers", angle = 90), fill = "white", label.size = NA, size = 8)
row2 <- ggplot(seg) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_void() + geom_richtext(aes(x = 1 ,y = 1, label = "Exhausted markers", angle = 90), fill = "white", label.size = NA, size = 8)
row3 <- ggplot(seg) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_void() + geom_richtext(aes(x = 1 ,y = 1, label = "Naïve markers", angle = 90), fill = "white", label.size = NA, size = 8)
### Combination
row <- ggarrange(row1, row2, row3, ncol = 1, nrow = 3, widths = c(3,3,3), heights = c(1,1,1))

### FeaturePlot
fig_2c <- FeaturePlot(tcell, features = tcell_feature, reduction = 'umap',ncol = 2, order = T, combine = T, 
                      cols = c('grey80',"#D4524E")) & NoAxes() & NoLegend()
fig_2c <- fig_2c  & theme(plot.title = element_text(size = 24, face = "italic"))

### Lastly, we created the shared legend
fig2c_leg <- FeaturePlot(tcell, features = tcell_feature, order = T) + 
  scale_color_gradient(low = 'grey80',high = '#D4524E',breaks= c(0,3), labels = c('min', 'Max')) + 
  theme(legend.position = 'bottom', legend.text = element_text(size = 16))

fig2c_leg <- get_only_legend(fig2c_leg)
### Slightly modify the arragement to fit the main figure
fig2c_leg <- ggarrange(blankPlot, fig2c_leg, ncol = 2, nrow = 1, widths = c(6.5,1), heights = c(0.1,0.1))

### Plot arrangement
fig_2c <- ggarrange(fig_2c, fig2c_leg, ncol=1, nrow = 2, widths = c(10, 10), heights = c(9, 0.5))

### Save
png(
  filename  = paste(work_path, "2C.png", sep=""),
  width     = 10,
  height    = 11,
  unit = 'in',
  res = 300
)

fig_2c

dev.off()
