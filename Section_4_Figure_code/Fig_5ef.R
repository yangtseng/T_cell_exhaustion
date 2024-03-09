########################################
### Section 4, Code for main figures ###
########################################

########################
### Figure 5e and 5f ###
########################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
load("./external_validation/Human_ST_subT.RData")
### Including P1T.subT, P3T.subT, P5T.subT, P7T.subT, P8T.subT, P9T.subT, P10T.subT and P11T.subT

### Sub figures (right)
top <- SpatialFeaturePlot(P11T.sub, feature = 'RUNX2', slot = 'scale.data', image.alpha = 0,
                   crop = T, pt.size.factor = 2.5) + 
  scale_fill_gradientn(colours = c('#A7C7E0','white','#F8BAA3', '#D4524E'), 
                       breaks = c(-2,0,2,4), limits = c(-2,4), name = "RUNX2\nz-score") + 
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 20))

bottom <- SpatialFeaturePlot(P11T.sub, feature = 'Exhaust_enhanced_module', image.alpha = 0,
                   crop = T, pt.size.factor = 2.5) + 
  scale_fill_gradientn(colours = c('#387EB8','#76A3CC','white', '#EB9280','#D4524E'), 
                       breaks = c(-2, -1,0,1,2), limits = c(-2,2), name = "Exhaustion module\nz-score") + 
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 20))

### save
png(
  filename = paste('/home/rstudio/R/invivo/external_JofH_230126/runx2_pics/P11T.png',sep = ''),
  width = 15,
  height = 10,
  units = 'in',
  res = 300
)
a/b
dev.off()
