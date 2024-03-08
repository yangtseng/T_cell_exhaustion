########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 3a ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
HCC.tcell <- readRDS("./murine_tcell_Runx28.rds")
tcell <- readRDS("./Cellline/Cellline_motif4.rds")

####################################
### Step 1, DE TF identification ###
####################################

### Set default assay as TF
DefaultAssay(HCC.tcell) <- "TF"
DefaultAssay(tcell) <- "TF"

### Extract Top 10 DE-TF
### HCC.tcell
TF.markers.HCC.tcell <- FindMarkers(HCC.tcell, only.pos = TRUE, ident.1 = c('1','4'), ident.2 = c('2','7'), min.pct = 0.1, logfc.threshold = 0.02)
### Ordered by log2FC
TF.markers.HCC.tcell <- TF.markers.HCC.tcell[order(TF.markers.HCC.tcello$avg_log2FC, decreasing = T),]
### Top 10 TF
TF.markers.HCC.tcell <- TF.markers.HCC.tcell[1:10,]

### tcell
### Change identity
Idents(tcell) <- factor(x = Idents(tcell), levels = c('Naive_T', 'Active_T','Exhausted_T_72hr','Exhausted_T_96hr'))
TF.markers.tcell <- FindMarkers(tcell, only.pos = TRUE, ident.1 = c('Exhausted_T_72hr','Exhausted_T_96hr'), ident.2 = c('Active_T'), min.pct = 0.1, logfc.threshold = 0.1)
### Ordered by log2FC
TF.markers.tcell <- TF.markers.tcell[order(TF.markers.tcell$avg_log2FC, decreasing = T),]
TF.markers.tcell <- TF.markers.tcell[1:10,]

#####################################################
### Step 2, Generate average matrix from TF assay ###
#####################################################

### Obtained average expression matrix from both HCC.tcell and tcell seurat object
### HCC.tcell
TF.exp.HCC.tcell <- AverageExpression(HCC.tcell, assays = 'TF')[["TF"]]
TF.exp.HCC.tcell.sort <- TF.exp.HCC.tcell[rownames(TF.exp.HCC.tcell)%in% rownames(TF.markers.tcell),]
### Order by alphabet
TF.exp.HCC.tcell.sort <- TF.exp.HCC.tcell.sort[order(rownames(TF.exp.HCC.tcell.sort)),]
colnames(TF.exp.HCC.tcell.sort) <- paste0('C', colnames(TF.exp.HCC.tcell.sort))

### tcell
TF.exp.tcell <- AverageExpression(tcell, assays = 'TF')[["TF"]]
TF.exp.tcell.sort <- TF.exp.tcell[rownames(TF.exp.tcell) %in% rownames(TF.markers.tcell),]
### Order by alphabet
TF.exp.tcell.sort <- TF.exp.tcell.sort[order(rownames(TF.exp.tcell.sort)),]
colnames(TF.exp.tcell.sort) <- c('\nNaïve\nT cell','\nActive\nT cell','\nExhausted\nT cell (72hr)','\nExhausted\nT cell (96hr)')

#####################################
### Step 3, Heatmap visualization ###
#####################################

### HCC.tcell
p <- Heatmap(t(scale(t(TF.exp.HCC.tcell.sort))), rect_gp = gpar(col = 'white',lwd = 2), 
             # column_split = factor(avg$markertype, levels=unique(avg$markertype)), 
             # top_annotation = column_ha, column_title = NULL,
             cluster_rows = F, cluster_columns = F, column_names_rot = 0, column_names_centered = T,
             heatmap_height = unit(0.8, "cm")*11, 
             heatmap_width = unit(1.5, "cm")*11, 
             row_names_gp = gpar(fontsize = 16),
             column_names_gp = gpar(fontsize = 16),
             name = "Z-score",    
             col = colorRamp2(c(-2, 0, 2), c("#919FB7", "white", "#E58C82")),
             heatmap_legend_param = list(legend_direction = "horizontal", 
                                         legend_width = unit(6, "cm")), 
             show_heatmap_legend = F)

hm.HCC.tcell <- grid.grabExpr(draw(p))

### tcell
p <- Heatmap(t(scale(t(TF.exp.tcell.sort))), rect_gp = gpar(col = 'white',lwd = 2), 
             # column_split = factor(avg$markertype, levels=unique(avg$markertype)), 
             # top_annotation = column_ha, column_title = NULL,
             cluster_rows = F, cluster_columns = F, column_names_rot = 0, column_names_centered = T,
            heatmap_height = unit(0.8, "cm")*11, 
            heatmap_width = unit(1.5, "cm")*11, 
            row_names_gp = gpar(fontsize = 16),
            column_names_gp = gpar(fontsize = 14),
            name = "\nZ-score",    
            col = colorRamp2(c(-2, 0, 2), c("#919FB7", "white", "#E58C82")),              
            heatmap_legend_param = list(legend_direction = "horizontal", 
                                        labels_gp = gpar(fontsize = 12),
                                        title_gp = gpar(fontsize = 15),
                                        legend_width = unit(6, "cm")))
hm.tcell <- grid.grabExpr(draw(p.ex,  heatmap_legend_side = 'bottom'))

### Combine two plots
TF.h <- wrap_plots(list(hm.HCC.tcell, hm.tcell), ncol=1)

################################
### Step 4, Add segmentation ###
################################

seg <- data.frame(x = 0.98, xend = 0.98, y = 0.6, yend = 1.55)
col1 <- ggplot(seg) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_void() + 
  geom_segment(aes(x = 0.96, y = 1.55, xend = 0.98, yend = 1.55)) + 
  geom_segment(aes(x = 0.96, y = 0.55, xend = 0.97, yend = 0.55)) + 
  geom_segment(aes(x = 0.96, y = 0.65, xend = 0.97, yend = 0.65)) +
  geom_segment(aes(x = 0.97, y = 0.55, xend = 0.97, yend = 0.65)) + 
  geom_segment(aes(x = 0.97, y = 0.6, xend = 0.98, yend = 0.6))

### Combine the segmentations
lines <- plot_spacer() / col1 / plot_spacer() + plot_layout(heights = c(4,1,3), guides = "collect")

### Save
png(
  filename  = paste0(work_path, "Fig_3a.png"),
  width     = 9.5,
  height    = 9, 
  unit = 'in',
  res = 300
)

(TF.h | plot_spacer() | (plot_spacer()/ col1 / plot_spacer() + plot_layout(widths = c(1,1,1), heights = c(2.85,4.5,2)))) + 
  plot_layout(ncol = 3, widths = c(5,-1.1, 0.3))

dev.off()
