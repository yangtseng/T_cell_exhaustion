########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 2b ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
tcell <- readRDS("./Cellline/Cellline_motif4.rds")
Load("./Cellline/Cellline_DEGs.RData")

### We manually selected several marker genes
marker_genes <- c('Ccnd3','Lef1','Tcf7','Itga4','Klf3','Jak1','Sell','Il7r','Ccl4','Ccnd2','Xcl1','Igf1r','Ccl3','Irf4','Cdk6','Gzmb','Tox','Skap1')
GA_marker_genes <- c('Klf2','Klf3','Lef1','Tcf7','Klf13','Ddx5','Cxcr4','Actb','Cdk8','Lars2','Ccr6','Galnt3','Cdk6','Bcl2','Akat2','Dtx1','Zeb1','Bach1','Sppl3','Wnt5b')

### We extracted prefix of cell names
cell_cluster <- sub("\\_.*", "", colnames(counts))
cell_cluster <- as.factor(cell_cluster)

### Heatmap settings
### Color
col_fun <- colorRamp2(c(-3,-1.5, 0, 1.5, 3), c('#0E1638',"#BAD6E8", 'white','#FDCEB5','#D4524E'))

### Top annotation
top_ano <- HeatmapAnnotation(Group = cluster_info, col = list(Group = c('Naive_T' = '#8ed3c7', 'Active_T' = '#98A8BE', 'Exhausted_T_72hr' = '#CBA5B0', 'Exhausted_T_96hr' = '#f3877f')),
                             annotation_legend_param = list(nrow = 1, title = 'Group', show_annotation_name = F, show_legend = F,
                                                            title_gp = gpar(fontsize = 10, fontface = 'plain'), 
                                                            labels_gp = gpar(fontsize= 7),
                                                            labels = c('Naive T cell','Active T cell','Exhausted T cell (72hr)','Exhausted T cell (96hr)')))

### Row annotation for gene expression
row_ano_GE <- rowAnnotation(link = anno_mark(at = which(rownames(counts_scaled) %in% marker_genes), 
                                               labels = marker_genes, side = c('right'),
                                               labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

### Row annotation for gene activity
row_ano_GA <- rowAnnotation(link = anno_mark(at = which(rownames(GA_counts_scaled) %in% GA_marker_genes), 
                                               labels = GA_marker_genes, side = 'right',
                                               labels_gp = gpar(fontsize = 10), padding = unit(1, "mm")))

### Heatmap
RNA_heatmap <- Heatmap(counts_scaled, top_annotation = top_ano, column_split = cluster_info, column_title = NULL,
                       col = col_fun, show_column_names = FALSE, right_annotation = row_ano_GE, show_row_names = F, 
                       cluster_columns = F, cluster_rows = F, use_raster = F, show_heatmap_legend = F,
                       heatmap_legend_param = list(colorbar = 'continuous', 
                                                   title='z-score', direction = 'horizontal',
                                                   title_gp = gpar(fontsize = 10, fontface = 'plain'), 
                                                   labels_gp = gpar(fontsize=6), legend_width = unit(3, 'cm')))

GA_heatmap <- Heatmap(GA_counts_scaled, top_annotation = top_ano, column_split = cluster_info, column_title = NULL,
                      col = col_fun, show_column_names = FALSE, right_annotation = row_ano_GA, show_row_names = F,
                      cluster_columns = F, cluster_rows = F, use_raster = F, show_heatmap_legend = F,
                      heatmap_legend_param = list(colorbar = 'continuous', 
                                                  title='z-score', direction = 'horizontal',
                                                  title_gp = gpar(fontsize = 10, fontface = 'plain'), 
                                                  labels_gp = gpar(fontsize=6), legend_width = unit(3, 'cm')))

### save heatmap separately
png(
  filename  = paste0(work_path, "Fig_2b_GE.png"),
  width     = 3,
  height    = 5, 
  unit = 'in',
  res = 300
)

plot(RNA_heatmap, ht_gap = unit(-1,'cm'))

dev.off()

png(
  filename  = paste0(work_path, "Fig_2b_GA.png"),
  width     = 3,
  height    = 5, 
  unit = 'in',
  res = 300
)

plot(GA_heatmap, ht_gap = unit(-1,'cm'))

dev.off()
