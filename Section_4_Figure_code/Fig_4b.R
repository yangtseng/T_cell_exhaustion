########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 4b ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
tcell <- readRDS("./Cellline/Cellline_motif4.rds")

### Set default assay to chromVAR
DefaultAssay(tcell) <- 'chromvar'

### Save
png(
  filename  = "/home/rstudio/R/TcExh_fig/4B.png",
  width     = 7.5,
  height    = 4.5,
  unit = 'in',
  res = 300
)

FeaturePlot(tcell, reduction = 'umap', features = "MA0511.2", order = F,
            cols = c('grey80',"#D4524E")) + #NoAxes() +
  geom_segment(aes(x = -14, y = -10, xend = -12, yend = -10), arrow = arrow(length = unit(0.3, "cm"), type="closed"), lineend = 'round', size = 1) + #X axis arrow
  geom_segment(aes(x = -14, y = -10, xend = -14, yend = -8), arrow = arrow(length = unit(0.3, "cm"), type="closed"), lineend = 'round', size = 1) + #Y axis arrow
  theme_bw() + ggtitle("") + 
  theme(legend.position = 'right', panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0)) +
  xlab("UMAP-1") + ylab("UMAP-2") + labs(colour = 'ChromVAR') +
  theme(axis.title=element_text(size=18), legend.text = element_text(size = 18), legend.title = element_text(size = 20)) +
  theme(axis.title.x = element_text(hjust=0.05), axis.title.y=element_text(hjust= 0.05)) #+ labs(caption="Created by TTY", size = 2)

dev.off()
