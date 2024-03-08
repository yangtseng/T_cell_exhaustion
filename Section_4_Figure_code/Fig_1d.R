########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 1d ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
HCC.tcell <- readRDS("./murine_tcell_Runx28.rds")

### Heatmap for well-known T cell subtype markers
### We selected the marker genes manually
genelist <- c('Ctla4','Pdcd1','Tigit','Tox','Havcr2','Lag3', 'Gzmb','Gzmk','Cd7','Il7r','Tcf7','Ccna2','Ccnb2','Ifit1','Ifit3')

### Average the expression of each gene in each cluster
avg <- as.data.frame(AverageExpression(HCC.tcell, features = genelist, assay = 'RNA', group.by = 'seurat_clusters')[["RNA"]])
### Change the column names
colnames(avg) <- paste0('C',colnames(avg))
### Generate the type columns of marker
avg$markertype <- c(rep('Exhaustion',6), rep('Effector',3), rep('Memory',2), rep('Proliferation',2), rep('IFN response',2))
### Reorder the data frame
avg <- avg[,c(1,4,5,6,2,7,3,8,9,10)]

### Basic heatmap settings
### Column annotation
annotation_col = data.frame(
  GeneType = factor(c(rep('Exhaustion',6), rep('Effector',3), rep('Memory',2), rep('Proliferation',2), rep('IFN_response',2))),
  row.names = rownames(avg)
)

### Annotation color 
ann_colors = list(
  GeneType = c(Exhaustion = "#f3877f", Effector = "#c2e8bc", Memory = "#8ed3c7", Proliferation = '#e2e1c3', IFN_response = '#bd89bd')
)

### Heatmap annotation
column_ha <- HeatmapAnnotation(
  GeneType = anno_block(gp = gpar(fill = c('#f3877f','#c2e8bc','#8ed3c7','#e2e1c3','#bd89bd')),
                        labels = c("Exh", "Eff", "Mem", "Prol","IFN"), 
                        labels_gp = gpar(col = "black", fontsize = 18)),
  gp = gpar(col = 'white')
)

### Heatmap plotting
p <- Heatmap(scale(t(avg[,1:9])), rect_gp = gpar(col = 'white',lwd = 2), 
             column_split = factor(avg$markertype, levels=unique(avg$markertype)), 
             top_annotation = column_ha, column_title = NULL, row_title = "Clusters",
             cluster_rows = F, cluster_columns = F, column_names_rot = 45, 
             row_title_gp = gpar(fontsize = 18),
             row_names_gp = gpar(fontsize = 15),
             column_names_gp = gpar(fontsize = 15),
             heatmap_height = unit(1.2, "cm")*7, 
             heatmap_width = unit(1.5, "cm")*10, 
             name = "Z-score",    
             col = colorRamp2(c(-3, 0, 3), c("#919FB7", "white", "#E58C82")),
             heatmap_legend_param = list(legend_direction = "horizontal", labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 15), legend_width = unit(6, "cm")))

### Saving
png(
  filename  = paste(work_path, "Fig_1d.png", sep=""),
  width     = 7,
  height    = 4.5,
  unit = 'in',
  res = 300
)

draw(p,  heatmap_legend_side = 'bottom')

dev.off()
