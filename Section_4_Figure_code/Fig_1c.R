########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 1c ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load data
HCC.tcell <- readRDS(paste0(work_path, "murine_tcell_Runx28.rds"))

###################################################
### Step 1, Create doughnuts for each condition ###
###################################################

### Cluster result: HCC.tcell@meta.data[["seurat_clusters"]]
### Group result: HCC.tcell@meta.data[["group"]]
df.doughnut <- as.data.frame(prop.table(table(Idents(HCC.tcell), HCC.tcell$group), margin = 2))
df.doughnut.label <- as.character(round(df.doughnut$Freq*100, digits = 1))

### We left only top5 value in df.doughnut.label
df.doughnut.label[c(5, 7:9, 12,15, 17:19, 25:27, 29, 33, 35,36)] <- ""

### We self defined doughnut plot function
plot_doughnut <- function(df.doughnut, df.doughnut.label, section){
  ggplot(df.doughnut[section,], aes(x = 1, y = Freq, fill = Var1)) + facet_wrap(~Var2) + 
    geom_col() +
    coord_polar(theta = "y") +
    geom_text(aes(label = df.doughnut.label[section]), size = 3, position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c('#f3877f','#c2e8bc','#8ed3c7','#d5a3b1','#e1b27a', '#c2c0d6', '#94a8c2','#e2e1c3', '#bd89bd')) +
    xlim(c(0, 4)) + 
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.position = 'none',
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(), 
          strip.background = element_rect(fill="white", color = 'white'), 
          strip.text.x = element_blank())
}

doughnut.1 <- plot_doughnut(df.doughnut, df.doughnut.label, 1:9)
doughnut.2 <- plot_doughnut(df.doughnut, df.doughnut.label, 10:18)
doughnut.3 <- plot_doughnut(df.doughnut, df.doughnut.label, 19:27)
doughnut.4 <- plot_doughnut(df.doughnut, df.doughnut.label, 28:36)

##################################
### Step 2, UMAP visualization ###
##################################

### Dimplot split into 4 groups w/ 2 cols
HCC.tcell <- SetIdent(HCC.tcell, value = 'seurat_clusters') 
levels(HCC_tcell$group) <- c('short term isotype', 'short term anti-PD1', 'long term isotype', 'long term anti-PD1')

### Plotting
png(
  filename  = paste0(work_path, "Fig_1c.png"),
  width     = 9.5,
  height    = 8,
  unit = 'in',
  res = 300
)

DimPlot(HCC.tcell, reduction = 'umap', group.by = 'ident', split.by = 'group', ncol = 2,
        cols =c('#f3877f','#c2e8bc','#8ed3c7','#d5a3b1','#e1b27a', '#c2c0d6', '#94a8c2','#e2e1c3', '#bd89bd')) + #NoAxes() +
  ### Add x, y-axis arrows manually
  geom_segment(aes(x = -5.5, y = -8, xend = -4, yend = -8), arrow = arrow(length = unit(0.3, "cm"), type="closed"), size = .8) + #X axis arrow
  geom_segment(aes(x = -5.5, y = -8, xend = -5.5, yend = -6.3), arrow = arrow(length = unit(0.3, "cm"), type="closed"), size = .8) + #Y axis arrow
  xlab("UMAP-1") + ylab("UMAP-2")  + labs(c('short term isotype', 'short term anti-PD1', 'long term isotype', 'long term anti-PD1')) +
  ### Theme settings
  theme_bw() + ggtitle("")  +
  theme(legend.position = 'right', panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  
  theme(legend.text = element_text(size = 24), legend.title = element_text(size = 28)) +
  theme(axis.title.x = element_text(hjust=0.02, size = 18), axis.title.y=element_text(hjust= 0.02, size = 18)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0)) +
  guides(color = guide_legend(override.aes = list(size=3), title = 'Clusters')) + 
  theme(strip.background = element_rect(fill="white", color = 'white'), strip.text.x = element_text(size = 28, colour = "black")) + # modify the facet_grid
  ### Insert four doughnut plots for each condition
  inset_element(doughnut.1, left = 0.15, bottom = 0.38, right = 0.6, top = 0.93) +
  inset_element(doughnut.2, left = 0.67, bottom = 0.38, right = 1.12, top = 0.93) +
  inset_element(doughnut.3, left = 0.15, bottom = -0.22, right = 0.6, top = 0.42) +
  inset_element(doughnut.4, left = 0.67, bottom = -0.22, right = 1.12, top = 0.42)

dev.off()
