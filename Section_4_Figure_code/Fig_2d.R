########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 2d ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
tcell <- readRDS("./Cellline/Cellline_motif4.rds")

### The violin plots were based on the module score of exhaustion and effectory
Eff <- VlnPlot(tcell, features = c('effectory1'), pt.size = 0, cols = c('#9DD2C6','#98A8BE','#CBA5B0','#E98A82')) + xlab('') + ylab('Effectory Module Score') +
  ggtitle('') + theme(legend.position = 'none', axis.line = element_line(), 
                      axis.title.y = element_text(hjust = 0.3, size = 20), axis.text.y = element_text(size = 18),
                      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18)) +
  coord_capped_cart(bottom=capped_horizontal(), left=capped_vertical(capped="both")) + 
  scale_x_discrete(labels=c('Naive\nT cell', 'Active\nT cell', 'Exhausted\nT cell (72hr)', 'Exhausted\nT cell (96hr)'))

Exh <- VlnPlot(tcell, features = c('exhaustion1'), pt.size = 0, cols = c('#9DD2C6','#98A8BE','#CBA5B0','#E98A82')) + xlab('') + ylab('Exhaustion Module Score') +
  ggtitle('') + theme(legend.position = 'none', axis.line = element_line(), 
                      axis.title.y = element_text(hjust = 0.5, size = 20), axis.text.y = element_text(size = 18),
                      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 18)) +
  coord_capped_cart(bottom=capped_horizontal(), left=capped_vertical(capped="both")) + 
  scale_x_discrete(labels=c('Naive\nT cell', 'Active\nT cell', 'Exhausted\nT cell (72hr)', 'Exhausted\nT cell (96hr)'))

### Save
png(
  filename  = paste(work_path, "Fig_2d.png", sep=""),
  width     = 13.5,
  height    = 4,
  unit = 'in',
  res = 300
)

Eff + Exh + plot_layout(ncol = 2)

dev.off()
