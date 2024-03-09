########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 5c ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
load("./external_validation/Human_bulk.RData"))
### Including gsva.runx2 and colAnn_df

### Boxplot function
### Draw boxplot for visualization
boxplot <- function(data, group, range, map_level){
  data <- as.data.frame(t(data))
  colnames(data) <- c('data')
  group <- group[match(rownames(data), rownames(group)),]
  data$group <- group[,range]
  data$group4 <- group[,range-1]
  p <- ggplot(data, aes(x=factor(group, level=c('Response', 'Non-Response')), y=data)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = group4), position=position_jitter(width=0.05)) + 
    geom_signif(comparisons = list(c('Response','Non-Response')), test = wilcox.test, test.args = c('less'),
                map_signif_level = map_level, textsize = 6) + ylim(-0.6, 0.68) +
    theme_bw() + ylab('GSVA enrichment score (RUNX2)') + xlab('HCC Patient') +
    scale_color_manual(values = c('CR' = '#90DB81','PR' = '#C5E6C0','SD' = '#F4B9B5','PD' = '#ED7268'), name = 'Patient Type') +
    theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20, hjust = 1),
          legend.title = element_text(size = 16), legend.text = element_text(size = 15))
  return(p)
}

### Plot and save
png(
  filename = paste0(work_path, "Fig_5c.png",
  width = 6.2,
  height = 4.5,
  units = 'in',
  res = 300
)

boxplot(gsva.runx2, colAnn_df, 2, F)

dev.off()
