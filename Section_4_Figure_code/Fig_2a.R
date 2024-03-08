########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 2a ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
tcell <- readRDS("./Cellline/Cellline_motif4.rds"))

### Modify the group names
group <- c('Active T cell','Exhausted T cell (72 hr)','Exhausted T cell (96 hr)', 'Naïve T cell')
names(group) <- levels(tcell)
tcell$group <- tcell$orig.ident
tcell <- RenameIdents(tcell, group)
Idents(tcell) <- factor(x = Idents(tcell), levels = c('Naïve T cell', 'Active T cell','Exhausted T cell (72 hr)','Exhausted T cell (96 hr)'))

### Set the default assay to SCT
DefaultAssay(tcell) <- 'SCT'

### Plot the UMAP
png(
  filename  = paste(work_path, "Fig_2a.png", sep=""),
  width     = 8,
  height    = 3.5,
  unit = 'in',
  res = 300
)

DimPlot(tcell, reduction = 'umap', group.by = 'ident', 
        cols = c('#8ed3c7','#98A8BE','#CBA5B0','#f3877f')) + #NoAxes() +
  geom_segment(aes(x = -17.5, y = -10, xend = -14.8, yend = -10), arrow = arrow(length = unit(0.3, "cm"), type="closed"), lineend = 'round', size = 1) + #X axis arrow
  geom_segment(aes(x = -17.5, y = -10, xend = -17.5, yend = -7.4), arrow = arrow(length = unit(0.3, "cm"), type="closed"), lineend = 'round', size = 1) + #Y axis arrow
  theme_bw() + ggtitle("") + 
  theme(legend.position = 'right', panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0)) +
  guides(color = guide_legend(nrow=4, byrow=F, override.aes = list(size=3), title = 'Groups')) +
  xlab("UMAP-1") + ylab("UMAP-2") +
  theme(axis.title=element_text(size=18), legend.text = element_text(size = 18), legend.title = element_text(size = 24)) +
  theme(axis.title.x = element_text(hjust=0.05), axis.title.y=element_text(hjust= 0.05)) #+ labs(caption="Created by TTY", size = 2)

dev.off()
