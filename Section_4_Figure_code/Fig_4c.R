########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 4c ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
tcell <- readRDS("./Cellline/Cellline_motif4.rds")

###############################
### Step 1, Data imputation ###
###############################

### We use MAGIC for data imputation
Runx2_magic <- magic(t(tcell@assays[["SCT"]]@data), knn=20, genes=c("Runx2"))
Runx2_magic <- as.data.frame(Runx2_magic$result$Runx2)
rownames(Runx2_magic) <- Runx2_MAGIC$params$data@Dimnames[[1]]
colnames(Runx2_magic_scale) <- 'Runx2'

### Normalized and centralized the Runx2_gep
Runx2_magic_scale <- as.data.frame(scale(Runx2_magic))

##############################################################
### Step 2, Extraction of chromatin accessibility activity ###
##############################################################

DefaultAssay(tcell) <- 'chromvar'
Runx2_cap <- as.data.frame(tcell@assays[["chromvar"]]@data)
Runx2_cap <- as.data.frame(t(Runx2_cap[rownames(Runx2_cap) == 'MA0511.2',])) ### MA0511.2 -> motif of Runx2
colnames(Runx2_cap) <- "Runx2"

### Normalized and centralized
Runx2_cap_scale <- as.data.frame(scale(Runx2_cap$Runx2))

#########################################################################
### Step 3, Generate the combined dataset including RNA and chromatin ###
#########################################################################

Runx2 <- rbind(Runx2_gep_scale, Runx2_cap_scale)
Runx2$type <- c(rep('gene_expression', 41594),rep('chromatin_accessibility', 41594))
Runx2$group <- sub("_.*", "", rownames(Runx2))
Runx2$numtype <- as.integer(factor(Runx2$group, levels = c("naive", "active", "exh72", 'exh96')))

### Plot 
p <- ggplot(Runx2, mapping = aes(x = group, y = Runx2, color = type)) +
  geom_boxplot(data = Runx2, lwd = 1.5, width=0.8, outlier.shape = NA) + ylim(-3, 3) +
  scale_x_discrete(limits=c("naive", "active", "exh72", 'exh96')) 

### Add lines
p <- p + stat_summary(fun = median,
               geom = "line",
               aes(group = type),
               col = colour, 
               position = position_dodge(width = .8), size = 3)
p <- p + stat_summary(fun = "median", aes(group = type), col = colour, 
                      size = 4, geom = "point", position = position_dodge(width = .8))

### Theme modification
p <- p + theme_bw() + scale_x_discrete(limits=c("naive", "active", "exh72", 'exh96'), 
                                       labels = c('Naïve T cell', 'Active T cell','Exhausted T cell (72hr)', 'Exhausted T cell (96hr)')) + 
  xlab("\nGroups") + ylab("Normalized values") + 
  scale_color_manual(values = c('#CCA4AF', '#B6C1E0'), name = "Types", 
                     labels = c('Chromatin Accessibility', 'Gene Expression')) +
  theme(axis.title = element_text(size = 28), axis.text = element_text(size = 20, color = 'black'), 
        axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
        legend.position = "top", legend.title = element_text(size = 24), 
        legend.justification = 'right',
        legend.text = element_text(size = 20), panel.grid = element_blank()) +
  guides(color=guide_legend(nrow=2))

### Save
png(
  filename  = paste0(work_path, 'Fig_4c.png'),
  width     = 4.5,
  height    = 7, 
  unit = 'in',
  res = 300
)

p

dev.off()


