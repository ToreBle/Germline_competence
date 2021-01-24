library("ggplot2")
library("tidyverse")
library("rGREAT")
library("GenomicRanges")
library("ggpubr")

dir <- "PATH/Germline_competence/"

#Load gene sets and H3K27ac peaks
sets <- c("EpiLC", "EpiSC")
for(set in sets){
  gene_set <- read.table(paste0(paste0(dir,"/scRNAseq_PGCLC_differentiation/",set,"_genes.txt"), sep=""), h=F)
  colnames(gene_set) <- c("id", "gene")
  peaks <- read.table(paste0(paste0(dir,"Enhancer_definition/H3K27ac_peaks_", set, ".bed"), sep=""), h=F)
  
  #Using GREAT (McLean et al. 2010) to link the enhancers to the proximal genes (all genes)
  great_genes_Enhancer <- submitGreatJob(peaks, species = "mm10", rule="twoClosest", adv_twoDistance = 500, version = "4")
  
  pdf(paste0(paste0(dir,"GREAT_",set, "_enhancer_distance.pdf"), sep=""), height = 4, width = 5)
  assignments <- plotRegionGeneAssociationGraphs(great_genes_Enhancer, type=2)
  dev.off()
  
  #Extract the assignments and keep enhancers that are linked to stage specific genes with a minimum distance of 3.5 kb to the TSS
  assignments <- as.data.frame(na.omit(assignments))
  linked_enhancer <- inner_join(assignments, gene_set, by="gene")
  linked_enhancer <- subset(linked_enhancer, linked_enhancer$distTSS>3500 | linked_enhancer$distTSS< (-3500))
  print(paste(set, "enhancers:", nrow(linked_enhancer), sep=" "))
  write.table(data.frame(linked_enhancer), row.names = F, col.names = F, quote = FALSE, sep = "\t", file=paste0(paste0(dir,"Enhancer_definition/",set,"_linked_enhancers.txt"), sep=""))
  write.table(data.frame(linked_enhancer[,1:3]), row.names = F, col.names = F, quote = FALSE, sep = "\t", file=paste0(paste0(dir,"Enhancer_definition/",set,"_linked_enhancers.bed"), sep=""))
  
}
