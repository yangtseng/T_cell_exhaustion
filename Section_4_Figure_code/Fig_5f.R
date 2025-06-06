########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 5f ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
load("./external_validation/Human_ST_annotated.RData"))
### Including P1T, P3T, P5T, P7T, P8T, P9T, P10T, P11T
load("./external_validation/Human_ST_sub.RData")
### Including P1T.sub, P3T.sub, P5T.sub, P7T.sub, P8T.sub, P9T.sub, P10T.sub and P11T.sub

### figure (left)
png(
  filename = paste0(work_path, "Fig_5f-1.png"),
  width = 12,
  height = 7,
  units = 'in',
  res = 300
)

SpatialDimPlot(P9T, label = F) +
  scale_fill_manual(values = c('#C7D2BD','#CB7D73','#C2C0D4','#98A8BC','#E3E2C7')) +
  theme(legend.key.size = unit(.8, 'cm'), legend.text = element_text(size = 12),
        legend.title = element_text(size = 15)) + 
  guides(fill = guide_legend(override.aes = list(size = 5), title = 'Cell type'))

dev.off()

### Sub figures (right)
top <- SpatialFeaturePlot(P9T.sub, feature = 'RUNX2', slot = 'scale.data', image.alpha = 0,
                   crop = T, pt.size.factor = 2.5) + 
  scale_fill_gradientn(colours = c('#A7C7E0','white','#F8BAA3', '#D4524E'), 
                       breaks = c(-2,0,2,4), limits = c(-2,4), name = "RUNX2\nz-score") + 
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 20))

bottom <- SpatialFeaturePlot(P9T.sub, feature = 'Exhaust_enhanced_module', image.alpha = 0,
                   crop = T, pt.size.factor = 2.5) + 
  scale_fill_gradientn(colours = c('#387EB8','#76A3CC','white', '#EB9280','#D4524E'), 
                       breaks = c(-2, -1,0,1,2), limits = c(-2,2), name = "Exhaustion module\nz-score") + 
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 20))

### save
png(
  filename = filename = paste0(work_path, "Fig_5f-2.png")
  width = 15,
  height = 10,
  units = 'in',
  res = 300
)

top/bottom

dev.off()

### This code can be modified and fit the supplementary figure [Supp. Fig. 9]
