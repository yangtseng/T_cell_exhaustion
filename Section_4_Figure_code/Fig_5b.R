########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 5b ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
load("./external_validation/murine_microarray.RData")
### Including runx2_exprs_sum.scaled, gsva.runx2_sum, Runx2_module and strain_runx2

### Boxplot
### Draw boxplot for visualization
boxplot <- function(data, group, range, map_level){
  data <- as.data.frame(t(data)[range])
  colnames(data) <- c('data')
  data$group <- factor(group[range], levels = c('Placebo','Anti-PD1'),ordered = TRUE)
  p <- ggplot(data, aes(x=group, y=data)) + 
    geom_boxplot(fill = c('#A8D0C7','#D78D82'), outlier.shape = NA, size = .8) + 
    geom_point(aes(fill = group), pch = 21, stroke = 1.5,
               position=position_dodge2(width=0.6), size = 2.5) +
    geom_signif(comparisons = list(c("Placebo", "Anti-PD1")), annotations = "*", 
                textsize = 10, size = 0.8, vjust = 0.5, y_position = 0.65) + 
    scale_fill_manual(values = c('#A8D0C7','#D78D82')) +
    scale_y_continuous(breaks = c(-0.6,-0.3,0,0.3,0.6)) +
    theme_bw() + ylab('GSVA enrichment score (Runx2)') + xlab('Treatment') + 
    theme(text = element_text(family = 'Arial'),
          axis.text.x = element_text(size = 20, color = 'black'), 
          axis.text.y = element_text(size = 18, color = 'black'),
          axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20,hjust = 1),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          legend.title = element_text(size = 14), legend.text = element_text(size = 12),
          legend.position = "None") +
    coord_cartesian(ylim=c(-0.7,0.75))
  return(p)
}

### Plot and save
png(
  filename = paste0(work_path, "Fig_5b.R",
  width = 4.5,
  height = 4.5,
  units = 'in',
  res = 300
)

boxplot(gsva.runx2_sum, strain_runx2, c(1:7), F)

dev.off()
