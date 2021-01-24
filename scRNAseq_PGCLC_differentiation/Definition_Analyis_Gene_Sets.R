library("ggplot2")
library("Seurat")
library("ggpubr")
library("tidyverse")

dir <- "PATH/Germline_competence/"

#Load complete UMI count matrix and define cell stages
Expression <- read.csv(paste0(dir,"scRNAseq_PGCLC_differentiation/UMI_counts.csv"), row.names = 1)

ESC <- Expression[,c(grep(".1", names(Expression)))]
SL <- Expression[,c(grep(".7", names(Expression)))]
d1_EpiLC <- Expression[,c(grep(".3", names(Expression)))]
d2_EpiLC <- Expression[,c(grep(".2", names(Expression)))]
d3_EpiLC <- Expression[,c(grep(".8", names(Expression)))]
EpiSC <- Expression[,c(grep(".5", names(Expression)))]


#Define d2/d4 PGCLC cluster and the EB cluster (non-PGCLC cluster of d2/d4 EB)
d2_EB_all <- Expression[,c(grep(".4", names(Expression)))]
d4_EB_all <- Expression[,c(grep(".6", names(Expression)))]

#Define d2 PGCLC
d2_PGCLC_cluster <- read.csv(paste0(dir,"/scRNAseq_PGCLC_differentiation/d2_PGCLC_cluster.csv")) %>% dplyr::select(-row.mean, -X)
d2_PGCLC <- Expression %>% dplyr::select(colnames(d2_PGCLC_cluster))
d2_EB <- d2_EB_all %>% dplyr::select(-one_of(names(d2_PGCLC)))

#Define d4 PGCLC
d4_PGCLC_cluster <- read.csv(paste0(dir,"/scRNAseq_PGCLC_differentiation/d4_PGCLC_cluster.csv")) %>% dplyr::select(-row.mean, -X)
d4_PGCLC <- Expression %>% dplyr::select(colnames(d4_PGCLC_cluster))
d4_EB <- d4_EB_all %>% dplyr::select(-one_of(names(d4_PGCLC)))


#Prepare the metadata for Seurat
ES <- data.frame(id = colnames(ESC), Stage = "ESC", Sample= "ESC")
SL <- data.frame(id = colnames(SL), Stage = "others", Sample= "SL")
d1 <- data.frame(id = colnames(d1_EpiLC), Stage = "others", Sample= "d1 EpiLC")
d2 <- data.frame(id = colnames(d2_EpiLC), Stage = "EpiLC", Sample= "d2 EpiLC")
d3 <- data.frame(id = colnames(d3_EpiLC), Stage = "others", Sample= "d3 EpiLC")
SC <- data.frame(id = colnames(EpiSC), Stage = "EpiSC", Sample= "EpiSC")

PGC_negative1 <- data.frame(id = colnames(d2_EB), Stage = "negative", Sample= "d2 EB")
PGCLC1 <- data.frame(id = colnames(d2_PGCLC), Stage = "PGCLC", Sample= "d2 EB")
PGC_negative2 <- data.frame(id = colnames(d4_EB), Stage = "negative", Sample= "d4 EB")
PGCLC2 <- data.frame(id = colnames(d4_PGCLC), Stage = "PGCLC", Sample= "d4 EB")

scRNAseq_data_differentiation <- rbind(PGC_negative1, PGC_negative2, PGCLC1, PGCLC2, ES, SL, d1, d2, d3, SC)
rownames(scRNAseq_data_differentiation) <- scRNAseq_data_differentiation$id
save(scRNAseq_data_differentiation, file=paste0(dir,"/scRNAseq_PGCLC_differentiation/scRNAseq_data_differentiation.Rdata"))


#Seurat object to select upregulated genes for different stages
Seurat <- CreateSeuratObject(Expression, project = "PGCLC", meta.data = scRNAseq_data_differentiation)
Seurat <- NormalizeData(Seurat)
#Seurat <- ScaleData(Seurat)


#####################################
#Definition of the different gene sets
#####################################

for(set in c("ESC", "EpiLC", "EpiSC", "PGCLC")){
  if(set=="ESC"){
    up <- FindMarkers(Seurat, group.by = Seurat$Stage, ident.1 = "ESC", ident.2 = c("EpiLC","EpiSC"), test.use = "negbinom", only.pos = TRUE)
  }
  if(set=="EpiLC"){
    up <- FindMarkers(Seurat, group.by = Seurat$Stage, ident.1 = "EpiLC", ident.2 = c("ESC","EpiSC"), test.use = "negbinom", only.pos = TRUE)
  }
  if(set=="EpiSC"){
    up <- FindMarkers(Seurat, group.by = Seurat$Stage, ident.1 = "EpiSC", ident.2 = c("ESC","EpiLC"), test.use = "negbinom", only.pos = TRUE)
  }
  if(set=="PGCLC"){
    up <- FindMarkers(Seurat, group.by = Seurat$Stage, ident.1 = "PGCLC", ident.2 = "negative", test.use = "negbinom", only.pos = TRUE)
  }
  
  set_up <- subset(up, up$p_val_adj<0.005)
  set_up <- subset(set_up, set_up$pct.1>0.2)
  set_up <- subset(set_up, set_up$pct.2<0.4)
  
  
  library("biomaRt")
  ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")#, host = "http://mar2016.archive.ensembl.org")
  genes <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"), filters="ensembl_gene_id", values=rownames(set_up), mart=ensembl, useCache = FALSE)
  detach("package:biomaRt", unload=TRUE)
  
  genes <- dplyr::rename(genes, ens_gene = ensembl_gene_id, gene = external_gene_name)
  write.table(genes, row.names = F, col.names = F, quote = FALSE, sep = "\t", file=paste0(paste0(dir,"/scRNAseq_PGCLC_differentiation/",set,"_genes.txt"), sep=""))
  }



##################################
#Analysis of the defined gene sets
#################################


  ESC_genes <- read.table(paste0(paste0(dir,"/scRNAseq_PGCLC_differentiation/ESC_genes.txt"), sep=""), h=F)
  colnames(ESC_genes) <- c("id", "gene")
  EpiLC_genes <- read.table(paste0(paste0(dir,"/scRNAseq_PGCLC_differentiation/EpiLC_genes.txt"), sep=""), h=F)
  colnames(EpiLC_genes) <- c("id", "gene")
  EpiSC_genes <- read.table(paste0(paste0(dir,"/scRNAseq_PGCLC_differentiation/EpiSC_genes.txt"), sep=""), h=F)
  colnames(EpiSC_genes) <- c("id", "gene")
  PGCLC_genes <- read.table(paste0(paste0(dir,"/scRNAseq_PGCLC_differentiation/PGCLC_genes.txt"), sep=""), h=F)
  colnames(PGCLC_genes) <- c("id", "gene")

  pdf("Expression_dynamics_gene_sets.pdf", height=5, width = 8)
  gene_sets <- list(Naive=ESC_genes, Formative=EpiLC_genes, Primed=EpiSC_genes, PGCLC=PGCLC_genes)
  all <- data.frame()
  
  for(N in 1:length(gene_sets)){
    List <- list(ESC=ESC,d1_EpiLC=d1_EpiLC, d2_EpiLC=d2_EpiLC, d3_EpiLC=d3_EpiLC, EpiSC=EpiSC, d2_EB=d2_EB, d2_PGCLC=d2_PGCLC, d4_EB=d4_EB, d4_PGCLC=d4_PGCLC)
    genes <- gene_sets[[N]]
    rownames(genes) <- genes$id
  
    genes_expression <- data.frame()
    cell_expression <- data.frame()
    #noise <- data.frame()
  
  for(i in 1:length(List)){
    subset <- List[[i]][rownames(genes),]
    subset_genes <- data.frame(mean = rowMeans(subset), State = names(List[i]))
    subset_cell <- data.frame(mean = rowMeans(t(subset)), State = names(List[i]))
    subset_cell$set <- names(gene_sets[N])
    subset_genes$set <- names(gene_sets[N])
    subset_genes$id <- rownames(subset_genes)
    
    if(names(List[i])=="ESC"){ESC_ex <- subset}
    if(names(List[i])=="d1_EpiLC"){d1EpiLC_ex <- subset}
    if(names(List[i])=="d2_EpiLC"){d2EpiLC_ex <- subset}
    if(names(List[i])=="d2_PGCLC"){d2PGCLC_ex <- subset}
    if(names(List[i])=="d4_PGCLC"){d4PGCLC_ex <- subset}
    rownames(subset_genes) <- NULL
    genes_expression <- rbind(genes_expression, subset_genes)
    cell_expression <- rbind(cell_expression, subset_cell)
  }
  
  color_dynamics <- c("#0101DF", "#8904B1", "#DF0101", "#B40404", "#FF8000", "darkgray", "#01DF01", "gray", "#088A08")
  labels <- c("ESC" = "ESC", "d1_EpiLC" = "d1 EpiLC", "d2_EpiLC" = "d2 EpiLC", "d3_EpiLC" = "d3 EpiLC", "EpiSC" = "EpiSC", "d2_EB" = "d2 EB", "d2_PGCLC" = "d2 PGCLC", "d4_EB" = "d4 EB", "d4_PGCLC" = "d4 PGCLC")
  genes_expression$log <- log10(genes_expression$mean+1) 
  
  p<- ggboxplot(genes_expression, x = "State", y = "log", color = "State", add="jitter", outlier.shape=NA)+labs(x =" ", y=expression(paste("log"[2],"(Expression per gene + 1)", sep="")), title=names(gene_sets[N]))+scale_color_manual(values=color_dynamics)+theme_classic()+
    theme(legend.position="", axis.text.x= element_text(color = "black", size=18, angle=45, hjust=1), axis.text.y= element_text(color = "black", size=14), axis.line = element_line(colour = "white"), axis.title=element_text(size=18), plot.title = element_text(size=16))+
    scale_x_discrete(labels=labels)
  q <- ggboxplot(cell_expression, x = "State", y = "mean", color = "State", add="jitter", outlier.shape=NA)+labs(x =" ", y="Expression per cell", title = names(gene_sets[N]))+scale_color_manual(values=color_dynamics)+theme_classic()+
    theme(legend.position="", axis.text.x= element_text(color = "black", size=21, angle=45, hjust=1), axis.text.y= element_text(color = "black", size=18), axis.line = element_line(colour = "white"), axis.title=element_text(size=21), plot.title = element_text(size=16))+
    scale_x_discrete(labels=labels)
  print(p)
  print(q)
  all <- rbind(all, cell_expression)
}
dev.off()

pdf("Overview_gene_sets.pdf", height=4, width = 8)
overview <- all %>% group_by(set,State) %>% summarise(m=mean(mean), sd=var(mean))
overview$set <- factor(overview$set, levels=c("Naive", "Formative", "Primed", "PGCLC"))
overview$State <- factor(overview$State, levels = c("d4_PGCLC", "d2_PGCLC", "d4_EB", "d2_EB", "EpiSC", "d3_EpiLC", "d2_EpiLC", "d1_EpiLC", "ESC"))
plot_overview <- ggplot(overview, aes(set, State, size=log2(m+1), color=set)) + geom_point()+theme_classic()+scale_color_manual(values=c("#A4A4A4", "#848484", "#6E6E6E", "#088A08"))+
  theme(legend.position="right", axis.text.x= element_text(color = "black", size=12), axis.text.y= element_text(color = "black", size=12), axis.line = element_line(colour = "white"), axis.title=element_text(size=16), plot.title = element_text(size=16))+
  scale_y_discrete(labels=labels)+labs(size="log2 (UMI +1)", x="Gene set", y="Stage")+scale_x_discrete(labels=c("Naive"=paste("2i ESC\nn = ",nrow(ESC_genes),sep=""), "Formative"=paste("EpiLC \nn = ",nrow(EpiLC_genes),sep=""), "Primed"=paste("EpiSC \nn = ",nrow(EpiSC_genes),sep=""), "PGCLC"=paste("PGCLC \nn = ",nrow(PGCLC_genes),sep="")))+guides(color = FALSE)
print(plot_overview)
dev.off()
