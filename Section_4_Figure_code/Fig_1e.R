########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 1e ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
HCC.tcell <- readRDS(paste0(work_path, "murine_tcell_Runx28.rds"))

### Please refer to Section_1 Step_3_Module_score for related analysis code

####################################################
### Step 1, Generate the dataset for ridge plots ###
####################################################

rdp_exh <- cbind(as.data.frame(as.factor(HCC.tcell$seurat_clusters)), as.data.frame(HCC.tcell$exhaustion1), as.data.frame(HCC.tcell$checkpoint1),
                 as.data.frame(HCC.tcell$effectory1), as.data.frame(HCC.tcell$cell_type))

### Modify the names in the data frame
colnames(rdp_exh) <- c('y','exh','cp','eff', 'group')
rdp_exh$group <- sub('Effector T cell','Effector\n T cell', rdp_exh$group)
rdp_exh$group <- sub('Exhausted T cell','Exhausted\n T cell', rdp_exh$group)
rdp_exh$group <- sub('Memory T cell','Memory\n T cell', rdp_exh$group)
rdp_exh$group <- sub('Proliferative T cell','Proliferative\n T cell', rdp_exh$group)
rdp_exh$group <- sub('High IFN response T cell','High IFN\n response\n T cell', rdp_exh$group)
rdp_exh$y <- paste0('C', rdp_exh$y)

################################################
### Step 2, Ridge plot for each module score ###
################################################

### Exhaustion module
rdp_exh <- rdp_exh[order(rdp_exh$group, decreasing = T),]
group.n <- c('High IFN\n response\n T cell', 'Proliferative\n T cell', 'Memory\n T cell', 'Effector\n T cell', 'Exhausted\n T cell')
names(group.n) <- c('High IFN response T cell',  'Proliferative T cell', 'Memory T cell', 'Effector T cell', 'Exhausted T cell')
rdp_exh_p <- ggplot(rdp_exh, aes(y = y, fill = y)) + geom_density_ridges(aes(x = exh), scale = 0.9) +
  scale_y_discrete(expand = expansion(add = c(0.025, 0.9))) + scale_x_continuous(expand = expansion(add = c(0.01, 0.5))) +
  facet_grid(factor(group, levels=c('High IFN\n response\n T cell','Proliferative\n T cell', 'Memory\n T cell', 
                                    'Effector\n T cell', 'Exhausted\n T cell'))~., scales = "free", space = "free", switch = "y") + 
  theme_bw() + ylab("") + xlab('Exhaustion\nModule Score') +
  scale_fill_manual(values = c('#f3877f','#c2e8bc','#8ed3c7','#d5a3b1','#e1b27a', '#c2c0d6', '#94a8c2','#e2e1c3', '#bd89bd')) +
  theme(strip.placement = "outside", strip.background = element_rect(colour=NA, fill=NA), strip.text.y = element_blank(),
        axis.line.y.left   = element_line(color = 'black'), legend.position = 'none', panel.border = element_rect(fill = NA, color = NA), 
        panel.grid = element_line(color = NA), axis.title.x = element_text(size = 18, color = 'black', hjust = 0), 
        axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 14, color = 'black'))

### Effectory module
rdp_eff_p <- ggplot(rdp_exh, aes(y = y, fill = y)) + geom_density_ridges(aes(x = eff), scale = 0.9) +
  scale_y_discrete(expand = expansion(add = c(0.025, 0.9))) + scale_x_continuous(expand = expansion(add = c(0.01, 0.5))) +
  facet_grid(factor(group, levels=c('High IFN\n response\n T cell','Proliferative\n T cell', 'Memory\n T cell', 
                                    'Effector\n T cell', 'Exhausted\n T cell'))~., scales = "free", space = "free", switch = "y") + 
  theme_bw() + ylab("") + xlab('Effectory\nModule Score') +
  scale_fill_manual(values = c('#f3877f','#c2e8bc','#8ed3c7','#d5a3b1','#e1b27a', '#c2c0d6', '#94a8c2','#e2e1c3', '#bd89bd')) +
  theme(strip.placement = "outside", strip.background = element_rect(colour=NA, fill=NA), axis.line.y.left   = element_line(color = 'black'),
        legend.position = 'none', panel.border = element_rect(fill = NA, color = NA), panel.grid = element_line(color = NA),
        axis.title.x = element_text(size = 18, color = 'black', hjust = 0), axis.ticks = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_text(size = 14, color = 'black'), strip.text.y = element_blank())

### Checkpoint module
rdp_cp_p <- ggplot(rdp_exh, aes(y = y, fill = y)) + geom_density_ridges(aes(x = cp), scale = 0.9) +
  scale_y_discrete(expand = expansion(add = c(0.025, 0.9))) + scale_x_continuous(expand = expansion(add = c(0.01, 0.5))) +
  facet_grid(factor(group, levels=c('High IFN\n response\n T cell','Proliferative\n T cell', 'Memory\n T cell', 
                                    'Effector\n T cell', 'Exhausted\n T cell'))~., scales = "free", space = "free", switch = "y") + 
  theme_bw() + ylab("") + xlab('Checkpoint\nModule Score') +
  scale_fill_manual(values = c('#f3877f','#c2e8bc','#8ed3c7','#d5a3b1','#e1b27a', '#c2c0d6', '#94a8c2','#e2e1c3', '#bd89bd')) +
  theme(strip.placement = "outside", strip.background = element_rect(colour=NA, fill=NA), legend.position = 'none',
        panel.border = element_blank(), panel.grid = element_line(color = NA), axis.line.y.left   = element_line(color = 'black'),
        axis.text.y = element_text(size = 18, color = 'black', vjust = -1.5), axis.title.x = element_text(size = 18, color = 'black', hjust = 0),
        axis.text.x = element_text(size = 14, color = 'black'), axis.ticks = element_blank(), strip.text.y.left = element_text(hjust = 1, size = 18, angle=0))

###############################################
### Step 4, Bar plot for cell cycle scoring ###
###############################################

rdp_ccy <- as.data.frame(table(HCC.tcell$seurat_clusters, HCC.tcell$Phase))

### Modify the lable
rdp_ccy_ct <- c('Exhausted\n T cell','Effector\n T cell','Memory\n T cell', 'Exhausted\n T cell', 'Exhausted\n T cell', 'Exhausted\n T cell', 'Effector\n T cell',
                'Proliferative\n T cell','High IFN\n response\n T cell')
rdp_ccy$cell_type <- rep(rdp_ccy_ct, 3)
colnames(rdp_ccy) <- c('y','ccy', 'freq','group')
rdp_ccy$y <- paste0('C', rdp_ccy$y)
rdp_ccy$group <-  factor(rdp_ccy$group, levels = c('High IFN\n response\n T cell','Proliferative\n T cell', 'Memory\n T cell', 'Effector\n T cell', 'Exhausted\n T cell'))

### Reordering
rdp_ccy <- rdp_ccy[order(rdp_ccy$group),]

### Barplot
rdp_ccy_p <- ggplot(rdp_ccy, aes(x = ccy, y = freq, fill = y)) + geom_bar(stat='identity', width = 0.8) + 
  scale_y_discrete(expand = expansion(add = c(0.025, 0.9))) + scale_x_discrete(expand = expansion(add = c(0.01, 0.4))) +
  facet_nested(group + factor(y, levels = c('C9','C8','C3','C7','C2','C6','C5','C4','C1'))~., 
               scales = "free_x", space = "free", switch = "y") + theme_bw() + ylab("") + xlab('Cell Cycle\nPhase') +
  scale_fill_manual(values = c('#f3877f','#c2e8bc','#8ed3c7','#d5a3b1','#e1b27a', '#c2c0d6', '#94a8c2','#e2e1c3', '#bd89bd')) +
  theme(strip.placement = "outside", strip.background = element_rect(colour=NA, fill=NA), legend.position = 'none',
        panel.border = element_blank(), panel.grid = element_line(color = NA), axis.text.y = element_blank(),  strip.text.y = element_blank(),
        axis.title.x = element_text(size = 18, color = 'black', hjust = 0), axis.text.x = element_text(size = 14, color = 'black'), axis.ticks = element_blank())

#################################
### Step 4, Combine the plots ###
#################################

png(
  filename  = paste0(work_path, "Fig_1e.png", sep=""),
  width     = 9,
  height    = 11,
  unit = 'in',
  res = 300
)

rdp_cp_p + plot_spacer() + rdp_exh_p + plot_spacer() + rdp_eff_p + plot_spacer() + rdp_ccy_p +plot_layout(widths = c(4,-1.5,4,-1.5,4,-1.5,2), guides = "collect")

dev.off()
