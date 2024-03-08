########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 3c ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
PE <- read.csv(paste0(work_path, "Pathway_enrichment_Tex.csv"), header = 1)
gs.go.kegg <- escape::getGeneSets(species = 'Mus musculus', library = 'C2', subcategory = 'KEGG')

### Filter the KEGG pathway based on the result of pathway enrichment 
KEGG <- gs.go.kegg[names(gs.go.kegg) %in% gsub("-", "_", rownames(PE))]

### DEGenes that were regulated by Runx2
DEGs.runx2.both <- c('Ctla4','Il18rap','Nrp1','Lgals3','Klrk1','Rbpj','Stat3','Runx2')
DEGs.runx2.HCC.tcell <- c('AA467197','Capg','Havcr2','Icos','Id2','Lag3','Mdfic','Prf1','Rora','Sdcbp2','Sdf4','Uhrf2','Zeb2')
DEGs.runx2.tcell <- c('Fkbp5','Itpr1','Anxa6','Runx3','Ezh2','Ahnak','Il12rb2','Gzmb','Mctp2','Lclat1')

##################################
### Step 1, Stack form of KEGG ###
##################################
### Create stack KEGG 
KEGG_stack <- data.frame()
### For loop
for(i in names(KEGG)){
  a <- as.data.frame(KEGG[[i]]@geneIds)
  a$names <- rep(KEGG[[i]]@setName,nrow(a))
  KEGG_stack <- rbind(KEGG_stack, a)
}
colnames(KEGG_stack) <- c('genes','names')

KEGG_stack <- KEGG_stack[order(KEGG_stack$genes),]

################################
### Step 2, Pathways network ###
################################

### KEGG network
g <- igraph::graph_from_data_frame(KEGG_stack, directed=T)

V(g)$type <- bipartite_mapping(g)$type  ## Add the "type" attribute
g$layout = igraph::layout_with_fr(g) #create a circular layout

### To the network
V(g)$color <- ifelse(names(page_rank(g)$vector) %in% names(KEGG), "#F1DBE3", 
                     ifelse(names(page_rank(g)$vector) %in% DEGs.runx2.both, '#FDF7F1', 
                            ifelse(names(page_rank(g)$vector) %in% DEGs.runx2.HCC.tcell,'#E58C82',
                                   ifelse(names(page_rank(g)$vector) %in% DEGs.runx2.tcell, "#919FB7", "grey80"))))
V(g)$shape <- ifelse(V(g)$type, "square", "circle")
V(g)$size <- ifelse(V(g)$type, 1.5, 0.5)
V(g)$label <- ifelse(names(page_rank(g)$vector) %in% KEGG_AVG27$name, gsub("KEGG ", "", gsub("_", " ", names(page_rank(g)$vector))), NA)
V(g)$label.cex <- 1
E(g)$color <- "lightgray"

### Generate the KEGG network
pkegg <- MetamapsDB::ig2ggplot(g, dfOnly = F, labels = T)

### Modify the color
color_name <- names(page_rank(g)$vector)
color_name <- color_name[order(color_name)]
color_manual <- ifelse(color_name %in% names(KEGG), "#F1DBE3", 
                       ifelse(color_name %in% DEGs.runx2.both, '#FDF7F1', 
                              ifelse(color_name%in% DEGs.runx2.HCC.tcell,'#E58C82',
                                     ifelse(color_name%in% DEGs.runx2.tcell, "#919FB7", "grey80"))))
### Add shape argument
pkegg[["labels"]][['shape']] <- 'type'

### Change types of the network
pkegg_new_size <- names(page_rank(g)$vector)
pkegg_new_size <- pkegg_new_size[order(pkegg_new_size)]
pkegg_new_type <- names(page_rank(g)$vector)
pkegg_new_type <- pkegg_new_type[order(pkegg_new_type)]

pkegg_new_size <- ifelse(grepl("KEGG_", pkegg_new_size), 8, 6)
pkegg_new_type <- ifelse(grepl("KEGG_", pkegg_new_type), "0", "1")

### Plot

png(
  filename  = paste0(work_path, "Fig_3c.png"),
  width     = 20,
  height    = 12,
  unit = 'in',
  res = 300
)

pkegg + scale_color_manual(values = color_manual) + theme(legend.position = 'none') + 
  scale_shape_manual(values = pkegg_new_type) +
  scale_size_manual(values = pkegg_new_size)  +
  # Lines
  geom_segment(aes(x = 0.08, y = 0.48, xend = 0.132, yend = 0.63), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.13, y = 0.3, xend = 0.255, yend = 0.505), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.28, y = 0.12, xend = 0.318, yend = 0.458), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.33, y = 0.32, xend = 0.348, yend = 0.435), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.255, y = 0.82, xend = 0.326, yend = 0.66), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.442, y = 0.622, xend = 0.45, yend = 0.82), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.372, y = 0.422, xend = 0.55, yend = 0.4), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.41, y = 0.282, xend = 0.55, yend = 0.322), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.415, y = 0.49, xend = 0.56, yend = 0.47), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.762, y = 0.58, xend = 0.762, yend = 0.662), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.708, y = 0.04, xend = 0.71, yend = 0.12), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.772, y = 0.915, xend = 0.75, yend = 0.82), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.885, y = 0.815, xend = 0.895, yend = 0.72), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  geom_segment(aes(x = 0.915, y = 0.375, xend = 0.8, yend = 0.285), 
               lineend = 'round', size = 1.2, colour = 'grey50') +
  # Labels
  geom_richtext(aes(x = 0.06, y = 0.48, label = "JAK-STAT signaling pathway"), 
                fill = "white", label.size = NA, size = 10, angle = 0) +
  geom_richtext(aes(x = 0.085, y = 0.28, label = "Natural killer cell mediated cytotoxicity"), 
                fill = "white", label.size = NA, size = 10, angle = 0) +
  geom_richtext(aes(x = 0.27, y = 0.1, label = "Graft versus host disease"), 
                fill = "white", label.size = NA, size = 10, angle = 0) +
  geom_richtext(aes(x = 0.28, y = 0.3, label = "Autoimmune thyroid disease"), 
                fill = "white", label.size = NA, size = 10, angle = 0) +
  geom_richtext(aes(x = 0.63, y = 0.32, label = "Cell adhesion molecules CAMs"), 
                fill = "white", label.size = NA, size = 10, angle = 0) + 
  geom_richtext(aes(x = 0.6, y = 0.4, label = "Type I diabetes mellitus"), 
                fill = "white", label.size = NA, size = 10, angle = 0) +
  geom_richtext(aes(x = 0.65, y = 0.47, label = "Antigen processing and presentation"), 
                fill = "white", label.size = NA, size = 10, angle = 0) +
  geom_richtext(aes(x = 0.28, y = 0.84, label = "T cell receptor signaling pathway"), 
                fill = "white", label.size = NA, size = 10, angle = 0) +
  geom_richtext(aes(x = 0.48, y = 0.84, label = "Primary immunodeficiency"), 
                fill = "white", label.size = NA, size = 10, angle = 0) +
  geom_richtext(aes(x = 0.75, y = 0.8, label = "Glycolysis glucogenesis"), 
                fill = "white", label.size = NA, size = 10, angle = 0) +
  geom_richtext(aes(x = 0.92, y = 0.7, label = "Pentose phosphate pathway"), 
                fill = "white", label.size = NA, size = 10, angle = 0) + 
  geom_richtext(aes(x = 0.78, y = 0.58, label = "Protein export"), 
                fill = "white", label.size = NA, size = 10, angle = 0) +
  geom_richtext(aes(x = 0.78, y = 0.28, label = "Ribosome"), 
                fill = "white", label.size = NA, size = 10, angle = 0) +
  geom_richtext(aes(x = 0.75, y = 0.12, label = "Circadian rhythm mammal"), 
                fill = "white", label.size = NA, size = 10, angle = 0) +
  
  # gene labels
  geom_richtext(aes(x = 0.04, y = 0.715, label = "<i>Il12rb2</i>"), 
                fill = NA, label.size = NA, size = 4, angle = 0) +
  geom_richtext(aes(x = 0.065, y = 0.745, label = "<i>Stat3</i>"), 
                fill = NA, label.size = NA, size = 4, angle = 0) +
  geom_richtext(aes(x = 0.165, y = 0.405, label = "<i>Klrk1</i>"), 
                fill = NA, label.size = NA, size = 4, angle = 0) +
  geom_richtext(aes(x = 0.288, y = 0.432, label = "<i>Gzmb</i>"), 
                fill = NA, label.size = NA, size = 4, angle = 0) +
  geom_richtext(aes(x = 0.308, y = 0.432, label = "<i>Prf1</i>"), 
                fill = NA, label.size = NA, size = 4, angle = 0) +
  geom_richtext(aes(x = 0.37, y = 0.47, label = "<i>Ctla4</i>"), 
                fill = NA, label.size = NA, size = 4, angle = 0) +
  geom_richtext(aes(x = 0.405, y = 0.495, label = "<i>Icos</i>"), 
                fill = NA, label.size = NA, size = 4, angle = 0)

dev.off()
