###########################################
### Section 1, Anlaysis of murine model ###
###########################################

#######################################################
### Further analysis of Tex progenitor in cluster 4 ###
#######################################################

set.seed(1234)
work_path = "./"
source("requirements.R")

load(paste0(work_path, "murine_tcell_modulescore4.rds"))
### It will load a pre-processed seurat object of T cells from section 1, step 3

########################################
### Step 1, Subset of Tex progenitor ###
########################################

### subclustering
c4sub <- FindSubCluster(HCC_tcell, 4, "RNA_snn", subcluster.name = "sub.cluster", resolution = 0.5, algorithm = 1)
### Cluster 4 was splited into 5 distinct sub-clusters

#########################################################
### Step 2, Visualization of cluster 4 sub-clustering ###
#########################################################

### UMAP showing only the sub-clustering result of cluster 4 [Supp. Fig. 5a]
DimPlot(c4sub, group.by = 'sub.cluster', cols = c('grey90','grey90','grey90','#f3877f','#c2e8bc','#8ed3c7','#d5a3b1','#e1b27a','grey90','grey90','grey90','grey90','grey90')) + #NoAxes() +
  geom_segment(aes(x = -5.5, y = -8, xend = -4, yend = -8), arrow = arrow(length = unit(0.3, "cm"), type="closed"), size = .8) + #X axis arrow
  geom_segment(aes(x = -5.5, y = -8, xend = -5.5, yend = -6.5), arrow = arrow(length = unit(0.3, "cm"), type="closed"), size = .8) + #Y axis arrow
  theme_bw() + ggtitle("Sub-clustering result of cluster 4") + 
  theme(legend.position = 'right', panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0), axis.title=element_text(size=15), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 15),
        axis.title.x = element_text(hjust=0.05), axis.title.y=element_text(hjust= 0.05)) + 
  xlab("UMAP-1") + ylab("UMAP-2") + 
  guides(color = guide_legend(override.aes = list(size=1.5), title = 'sub-clusters')) 

### Tex progenitor markers: Ccr7, Cxcl10, Il7r, Slamf6, Tcf7 and Xcl1
### Dot plot revealed the expression profile of Tex progenitor markers in each cluster and each sub-cluster of cluster 4 [Supp. Fig. 5b]
DotPlot(c4sub, group.by = 'sub.cluster', cols = c('grey80',"#D4524E"),features = c('Ccr7','Cxcl10','Il7r','Slamf6','Tcf7','Xcl1')) + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 18)) +
  xlab('Gene symbol') + ylab('Sub-clusters')

### We extracted the Tex progenitor population (c4_4) which expressed all the Tex progenitor markers
saveRDS(c4sub, paste0(work_path, "HCC_tex_progenitor5.rds"))
