########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 4e ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
tcell <- readRDS("./Cellline/Cellline_motif4.rds")

# Dot plot for both DA and DE in Runx2 module
Runx2_gene <- c('Runx2','Ctla4','Il18rap','Klrk1','Lgals3','Nrp1','Rbpj','Stat3')

### Extract data from seurat
DefaultAssay(tcell) <- 'GA'
GA_data <- DotPlot(tell, features = Runx2_gene)$data

DefaultAssay(tcell) <- 'SCT'
GE_data <- DotPlot(tcell, features = Runx2_gene)$data

### Combine
dot_data <- rbind(GA_data, GE_data)

### Create annotation 
dot_data$label <- c(rep('Gene Activity',32), rep("Gene Expression", 32))

### Change id order
dot_data <- dot_data %>%
  arrange(id = factor(id, levels=c("Naive_T", "Active_T", "Exhausted_T_72hr", "Exhausted_T_96hr")))
         
### Create new id name
dot_data$name <- rep(c('Naive T cell','Active T cell','Exhausted T cell 72hr','Exhausted T cell 96hr'), each = 16)

dot_data <- dot_data %>%
  mutate(name = factor(name, levels = c('Exhausted T cell 96hr','Exhausted T cell 72hr','Active T cell','Naive T cell')))

### Plot and save
png(
  filename  = paste0(work_path, "Fig_4e.png",
  width     = 6,
  height    = 8,
  unit = 'in',
  res = 300
)

ggplot(dot_data, aes(x=features.plot, y = name, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_size('%\nDetected', range = c(0,10)) +
  facet_wrap(~label, ncol = 1) + theme_bw() +
  xlab('Gene Symbol') + ylab('') + labs(colour="Scaled\nExpression\nLevel") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.y = element_text(size = 16, colour = 'black'),
        axis.text.x = element_text(size = 16, face = 'italic', colour = 'black', angle = 40, hjust = 1), 
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 24),
        axis.title.x = element_text(vjust = -1),
        legend.text = element_text(size = 14),
        legend.position = 'bottom',
        legend.justification = "right",
        legend.direction = 'horizontal',
        legend.title = element_text(vjust = 0.2, size = 16, hjust = 0.95)) +
  scale_color_gradient2(breaks = c(-1,0,1), low = '#919FB7', high = '#E58C82')
  
dev.off()
