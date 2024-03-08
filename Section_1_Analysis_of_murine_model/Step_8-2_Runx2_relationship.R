###########################################
### Section 1, Anlaysis of murine model ###
###########################################

############################################################
### Runx2 relationship with markers of T cell exhaustion ###
############################################################

set.seed(1234)
work_path = "./"
source("requirements.R")

### Load seurat object
HCC.tcell <- readRDS(paste0(work_path, "murine_tcell_Pathway7.rds"))

##########################
### Runx2 module score ###
##########################

### Runx2 module was defined as the Runx2 and its downstream target genes
### Runx2, Ctla4, Nrp1, Lgals3, Rbpj, Stat3, Klrk1 and Il18rap
Runx2_module <- c('Runx2','Ctla4','Nrp1','Lgals3','Rbpj','Stat3','Klrk1','Il18rap')

### Module score calcuation
HCC.tcell <- AddModuleScore(HCC.tcell, features = list(Runx2_sig), name = "Runx2_module")

##########################################################################
### Explore the relationship between Runx2 and exhaustion module score ###
##########################################################################

### Extract the module score from seurat object
Tex_runx2 <- as.data.frame(HCC.tcell@meta.data[["exhaustion1"]])
Tex_runx2$runx2 <- HCC.tcell@meta.data[["Runx2_module1"]]
colnames(Tex_runx2) <- c('Exhaustion_module','Runx2_module')

### Join with the clsuter information
Tex_runx2$cluster <- HCC.tcell@meta.data[["seurat_clusters"]]

### Join with the T cell subtype information
Tex_runx2$subtype <- HCC.tcell@meta.data[["cell_type"]]

### Visualize the result using scatter plot [Supp. Fig. 8a]
ggplot(Tex_runx2, aes(x = Exhaustion_module, y = Runx2_module)) + 
  geom_point(size = 1, alpha = 0.5, aes(color = subtype)) + theme_bw() + 
  geom_smooth(method=lm, color = "grey40", size=0.5) + 
  stat_cor(method = "pearson", label.x = 0.1, label.y = 1.2) + 
  scale_color_manual(values = c('#f3877f','#c2e8bc','#8ed3c7','#e2e1c3', '#bd89bd')) +
  ylab('Runx2 signature') + xlab('T exhaustion signature')+ facet_wrap(~subtype, ncol=3) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), 
        panel.grid = element_blank(), legend.position = c(0.85, 0.25), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 15)) +
  guides(color = guide_legend(title = 'Cell types'))

### Visualize the result using UMAP [Supp. Fig. 8b]
FeaturePlot(HCC.tcell, features = "Runx2_module1", cols = c('grey80',"#D4524E"), order = F) + ggtitle('Runx2 signature')

### Visualize the result with violin plot [Supp. Fig. 8c]
VlnPlot(HCC.tcell, features = "Runx2_module1", group.by = 'seurat_clusters', pt.size = 0, 
        cols = c('#f3877f','#c2e8bc','#8ed3c7','#d5a3b1','#e1b27a', '#c2c0d6', '#94a8c2','#e2e1c3', '#bd89bd')) + 
  NoLegend() + xlab('Clusters') + ylab('Module score') + ggtitle('Runx2 signature') + 
  guides(fill=guide_legend(title="Clusters")) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

### save files
saveRDS(HCC.tcell, paste0(work_path, "murine_tcell_Runx28.rds"))
  
