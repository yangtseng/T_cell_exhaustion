########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 3b ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
PE <- read.csv(paste0(work_path, "Pathway_enrichment_Tex.csv"), header = 1)
gs.go.kegg <- escape::getGeneSets(species = 'Mus musculus', library = 'C2', subcategory = 'KEGG')

### Plot the result of pathway enrichment
plot_GSEA <- function(PE, gs, title){
  PE$name = rownames(PE)
  
  ### Add gene number    
  for(i in 1:nrow(PE)){
    a <- PE$name[i]
    j <- gsub("-","_", a)
    PE$genenum[i] <- length(gs[[j]]@geneIds)
  }  
  PE$genenum <- as.numeric(PE$genenum)
  
  ### Color with adjust_pval
  PE$adj_p <- -log10(PE$p_val_adj)
  
  ### Adjust name of each pathway
  PE$name <- gsub("KEGG-", "", PE$name)
  PE$name <- tolower(PE$name)
  PE$name <- Hmisc::capitalize(PE$name)
  PE$name <- gsub("-", " ", PE$name)
  PE$name <- gsub("cams", "CAMs", PE$name)
  PE$name <- gsub("Jak stat", "JAK-STAT", PE$name)
  PE$name <- gsub("Type i",'Type I', PE$name)
  PE$name <- gsub("And", "and", PE$name)
  PE %>%
    arrange(FE) %>%    ### First sort by val. This sort the dataframe but NOT the factor levels
    mutate(name=factor(name, levels=name)) %>% 
    ggplot(aes(x = FE, y = name, color = adj_p, fill = adj_p)) +
    geom_bar(stat = 'identity', width = .2) +
    geom_point(aes(x = FE, size = genenum, color = adj_p)) + 
    theme_bw() + xlim(c(0,1.7)) +
    scale_color_gradientn(colors = c("#919FB7", "#E58C82"),
                          limits = c(150, 300), 
                          na.value = "#E58C82") + 
    scale_fill_gradientn(colors = c("#919FB7", "#E58C82"),
                         limits = c(150, 300), 
                         na.value = "#E58C82") +
    scale_size_continuous(range = c(10,15)) + labs(size = 'Number of gene', color = '-log(adj_p_val)\n', fill = '-log(adj_p_val)\n') + 
    theme(panel.grid = element_blank(), axis.title = element_text(size = 32), axis.text = element_text(size = 28, color = 'black'), 
          plot.title = element_text(size = 28)) + 
    ylab('KEGG Gene set\n') + xlab('Fold Enrichment') + ggtitle(title) + 
    theme(legend.text = element_text(size = 24), legend.title = element_text(size = 28)) +
    guides(
      color = guide_colorbar(order = 0),
      size = guide_legend(order = 1)
    )
}

### Save
png(
  filename  = paste0(work_path, "Fig_3b.png"),
  width     = 14,
  height    = 12,
  unit = 'in',
  res = 300
)

plot_GSEA(PE, gs.go.kegg, 'C1, 4 vs. C2, 7')

dev.off()
