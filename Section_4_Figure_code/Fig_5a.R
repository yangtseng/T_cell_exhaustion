########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 5a ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
load("./external_validation/murine_microarray.RData")
### Including runx2_exprs_sum.scaled, gsva.runx2_sum, Runx2_module and strain_runx2

### Color
col_fun = colorRamp2(c(-2,-1, 0,1, 2), c('#377EB8',"#BAD6E8", 'white','#FDCEB5','#D4524E'))

### Plot heatmap and save
png(
  filename = paste0(work_path, "Fig_5a.png",
  width = 5,
  height = 4,
  units = 'in',
  res = 300
)

Heatmap(runx2_exprs_mix.scaled, cluster_rows = F, cluster_columns = T, show_column_names = T, show_row_names = T, 
        rect_gp = gpar(col = "white", lwd = 2), row_names_gp = gpar(fontsize = 15, fontface = 'italic'), col = col_fun,
        column_names_gp = gpar(fontsize = 15),  column_title_gp = gpar(fontsize = 18),
        heatmap_legend_param = list(title = 'z-score', title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 12)), 
        column_title = "Runx2 expression level \n in exhausted CD8+ T cells from mice")

dev.off()
