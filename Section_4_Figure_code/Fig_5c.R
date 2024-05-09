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
  p <- ggplot(data, aes(x=factor(group, level=c('Responder','Non-Responder')), y=data)) + 
    geom_boxplot(outlier.shape = NA, size = 0.8) + 
    geom_point(aes(fill = group4), pch = 21, stroke = 1.5,
               position=position_dodge2(width=0.5), size = 3) + 
    geom_signif(comparisons = list(c('Responder','Non-Responder')), annotations = "*", 
                textsize = 10, size = 0.8, vjust = 0.6, y_position = 0.62) + 
    theme_bw() + ylab('GSVA enrichment score (RUNX2)') + xlab('HCC Patient')+
    scale_y_continuous(breaks = c(-0.6,-0.3,0,0.3,0.6)) +
    scale_fill_manual(values = c('CR' = '#90DB81','PR' = '#C5E6C0','SD' = '#F4B9B5','PD' = '#ED7268'), name = 'Patient Type') +
    theme(text = element_text(family = 'Arial'),
          axis.text.x = element_text(size = 20, color = 'black', vjust = 1), 
          axis.text.y = element_text(size = 18, color = 'black'),
          axis.title.x = element_text(size = 20, vjust = -1), 
          axis.title.y = element_text(size = 20, hjust = 1),
          legend.title = element_text(size = 16), legend.text = element_text(size = 15),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(),) +
    coord_cartesian(ylim=c(-0.62,0.72)) +
    scale_y_continuous(breaks = c(-0.6,-0.3,0,0.3,0.6)) +
    scale_x_discrete(labels=c('Responder','Non-Responder'))
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
