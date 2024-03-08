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

rdp_exh <- cbind(as.data.frame(as.factor(HCC.tcell$seurat_clusters)), as.data.frame(HCC.tcell$exhaust1), as.data.frame(HCC.tcell$cp1),
                 as.data.frame(HCC.tcell$eff1), as.data.frame(HCC.tcell$cell_type))

### Modify the names in the data frame
colnames(rdp_exh) <- c('y','exh','cp','eff', 'group')
rdp_exh$group <- sub('Effector T cell','Effector\n T cell', rdp_exh$group)
rdp_exh$group <- sub('Exhausted T cell','Exhausted\n T cell', rdp_exh$group)
rdp_exh$group <- sub('Memory T cell','Memory\n T cell', rdp_exh$group)
rdp_exh$group <- sub('Proliferative T cell','Proliferative\n T cell', rdp_exh$group)
rdp_exh$group <- sub('High IFN response T cell','High IFN\n response\n T cell', rdp_exh$group)
rdp_exh$y <- paste0('C', rdp_exh$y)
