########################################
### Section 4, Code for main figures ###
########################################

#################
### Figure 4d ###
#################

set.seed(1234)
work_path = "./main_figures/"
source("requirements.R")

### Load latest file
tcell <- readRDS("./Cellline/Cellline_motif4.rds")

### Set default assay to peaks
DefaultAssay <- "peaks"

### Design the function for coverage plot
cover <- function(data, gene, ranges.show, ranges.label, sub, multiple){
  ranges.show$color <- "orange"
  CP <- CoveragePlot(
    object = data,
    region = gene,
    idents = idents.plot,
    extend.upstream = 2000,
    extend.downstream = 1000,
    region.highlight = ranges.show,
    heights = c(8,1,1,1),
    sep = c(":", "-")
  )
  CP <- CP & scale_fill_manual(values = c('#9DD1C6','#98A8BE','#CBA4B0','#E98A82'))
  CP <- CP & theme(strip.text.y.left = element_blank(),strip.background = element_blank(),
                   axis.title = element_text(size = 18), axis.text = element_text(size = 0),
                   legend.position = 'None') 
  if(multiple == 1){
    CP[[1]][[2]]$layers[[5]]$data$dodge <- c(0,2,3)
    CP[[1]][[2]]$layers[[5]]$aes_params$size <- 2
  }else if(multiple == 2){
    CP[[1]][[2]]$layers[[5]]$data$dodge <- c(0,1.5)
    CP[[1]][[2]]$layers[[5]]$aes_params$size <- 2
  }else{
    CP[[1]][[2]]$layers[[4]]$aes_params$size <- 0
  }
  linkplot <- CP[[1]][[4]]
  # to assign back into the assembled coverage plot
  CP[[1]][[4]] <- linkplot + scale_color_gradient2(low = '#919FB7', # whatever colors you want here
                                                       mid = "grey50",
                                                       high = '#E58C82')
  CP[[1]][[4]]$labels$x <- gsub('bp','kb',CP[[1]][[4]]$labels$x)
  CP[[1]][[4]]$labels$x <- gsub('chr','\nchr',CP[[1]][[4]]$labels$x)

  if(sub == T){
    CP <- CP & ylab('')
    CP[[1]][[1]]$labels$y <- ranges.label
    CP[[1]][[1]]$theme$axis.title <- element_text(size = 18)
  }
  CP[[1]][[1]]$labels$title <- gene
  CP[[1]][[1]]$theme$plot.title <- element_text(size = 22, face = 'italic')
  
  return(CP)
}

### We plot four genes as a group
### Runx2
ranges.show.runx2 <- StringToGRanges(c('chr17-44611802-44612806','chr17-44605855-44606630','chr17-44608540-44609596','chr17-44674555-44675444',
                                       'chr17-44670830-44671630','chr17-44727490-44729293','chr17-44556857-44557728'))
### Ctla4
ranges.show.ctla4 <- StringToGRanges(c("chr1-60898220-60900102", 'chr1-60908480-60909374','chr1-60910813-60911105', 
                                       'chr1-60911280-60911968','chr1-60913133-60913702', 'chr1-60886596-60887470'))
### Il18rap
ranges.show.il18rap <- StringToGRanges("chr1-40523650-40524669")
### Klrk1
ranges.show.klrk1 <- NULL

runx2 <- cover(Alltcellcombined_new, "Runx2", ranges.show.runx2, '(range 0-5700)', sub = T, multiple = 1)
ctla4 <- cover(Alltcellcombined_new, 'Ctla4', ranges.show.ctla4, '(range 0-3300)', sub = T, multiple = 0)
il18rap <- cover(Alltcellcombined_new, "Il18rap", ranges.show.il18rap, '(range 0-2600)', sub = T, multiple = 0)
klrk1 <- cover(Alltcellcombined_new, 'Klrk1', ranges.show.klrk1, '(range 0-5700)', sub = T, multiple = 0)

### Save
png(
  filename  = paste0(work_path, "Fig_4d-1.png"),
  width     = 10.5,
  height    = 4.5,
  unit = 'in',
  res = 300
)

runx2|ctla4|il18rap|klrk1

dev.off()

### Another four Runx2 regulated genes' coverage plots
### Lgals3
ranges.show.lgals3 <- StringToGRanges(c("chr14-47371596-47372036",'chr14-47373466-47374684'))
### Nrp1
ranges.show.nrp1 <- StringToGRanges(c("chr8-128390821-128391805",'chr8-128387061-128388131','chr8-128418387-128418873',
                                      'chr8-128420207-128421290','chr8-128480687-128481597','chr8-128483190-128485109'))
### Rbpj
ranges.show.rbpj <- StringToGRanges(c("chr5-53572172-53573256",'chr5-53565607-53567350'))
### Stat3
ranges.show.stat3 <- StringToGRanges(c("chr11-100933781-100934472"))

lgals3 <- cover(Alltcellcombined_new, "Lgals3", ranges.show.lgals3, '(range 0-1600)', sub = T, multiple = 0)
nrp1 <- cover(Alltcellcombined_new, 'Nrp1', ranges.show.nrp1, '(range 0-1900)', sub = T, multiple = 0)
rbpj <- cover(Alltcellcombined_new, "Rbpj", ranges.show.rbpj, '(range 0-6900)', sub = T, multiple = 0)
stat3 <- cover(Alltcellcombined_new, 'Stat3', ranges.show.stat3, '(range 0-7600)', sub = T, multiple = 2)

png(
  filename  = paste0(work_path, "Fig_4d-2.png"),
  width     = 10.5,
  height    = 4.5,
  unit = 'in',
  res = 300
)

lgals3|nrp1|rbpj|stat3

dev.off()

