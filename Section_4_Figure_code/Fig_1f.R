########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 1f ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
HCC.tcell <- readRDS("./murine_tcell_Runx28.rds")

### We manually selected six features to identify the Tex progenitor population
Tex_prog_feature <- c('Ccr7','Il7r','Slamf6','Tcf7', 'Sell','Havcr2')

### FeaturePlot for each gene
Tex_prog_fig <- list()
for(i in Tex_prog_feature){
  fig <- FeaturePlot(HCC.tcell, features = i, ncol = 1, order = F, cols = c('grey80',"#D4524E")) & NoAxes() & NoLegend()
  fig <- fig  & theme(plot.title = element_text(size = 32, face = "italic"))
  fig <- fig + ggtitle(i)
  Tex_prog_fig[[i]] <-fig
}

### We add two segmentation to separate the features into two group
### (1) Markers of TEX progenitor (+) and (2) Markers of TEX progenitor (-)
seg <- data.frame(x = 0.95, xend = 1.05, y = 0.3, yend = 0.3)
col1 <- ggplot(seg) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 2) + theme_void() + 
  geom_richtext(aes(x = 1 ,y = 0.3, label = "Markers of TEX progenitor (+)"), fill = "white", label.size = NA, size = 15, angle = 0)
col2 <- ggplot(seg) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 2) + theme_void() +
  geom_richtext(aes(x = 1 ,y = 0.3, label = "Markers of TEX progenitor (-)"), fill = "white", label.size = NA, size = 15, angle = 0)

### Combine the plots and save
png(
  filename  = paste0(work_path, "Fig_1f.png"),
  width     = 9.5,
  height    = 15,
  unit = 'in',
  res = 300
)

col1/(fig[['Ccr7']] | fig[['Il7r']])/(fig[['Slamf6']] | fig[['Tcf7']]) /col2/(fig[['Sell']] | fig[["Havcr2"]]) + plot_layout(heights = c(1,3,3,1,3),guides = "collect")

dev.off()
