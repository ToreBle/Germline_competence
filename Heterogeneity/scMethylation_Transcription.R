library(biomaRt)
library(edgeR)
library(dplyr)
library(ggpubr)
library(Rtsne)

dir <- "PATH/Germline_competence"

#scRNAseq data from https://doi.org/10.1038/s41586-019-1825-8
scRNA_counts <- read.table(paste0(dir,"Heterogeneity/GSE121650_rna_counts.tsv"), header=T, row.names = "ensembl_id")


#RPKM normalization
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl", host = "http://dec2016.archive.ensembl.org")
s=getBM(attributes=c("ensembl_gene_id", "start_position","end_position"), mart=ensembl, useCache = F)
s$length=s$end_position - s$start_position
detach("package:biomaRt", unload=TRUE)

genes <- s %>% dplyr::filter(ensembl_gene_id %in% rownames(scRNA_counts))
genes <- genes[order(rownames(scRNA_counts), genes$ensembl_gene_id),]
sc_rpkm <- rpkm(scRNA_counts, genes$length)


#Select only assigned epiblast cells from meta data: https://github.com/rargelaguet/scnmt_gastrulation
meta <- read.table(paste0(dir,"/Heterogeneity/sample_metadata.txt"), header=T)
meta_epi <- subset(meta, meta$lineage10x=="Epiblast" & meta$stage!="E7.5")
meta_epi$sample <- as.character(meta_epi$sample)


#tSNE plot of the scRNAseq data
set.seed(1989)
epi <- meta_epi %>% dplyr::filter(meta_epi$sample %in% colnames(sc_rpkm))
tsne_prep <- data.frame(sc_rpkm) %>% dplyr::select(epi$sample)
d <- dist(t(tsne_prep))
r <- Rtsne(d, is.matrix=T)
tnse_plot <- data.frame(x = r$Y[,1], y = r$Y[,2], State = epi$stage)
pdf("tSNE_Epiblast.pdf", height = 5, width = 5.5)
ggplot(tnse_plot, aes(y, x, colour = State))+geom_point(size=2) +  scale_color_manual(values = c("#0101DF", "#DF0101", "#FF8000"), "Stage") + 
  labs(x = "Component 2: t-SNE", y = "Component 1: t-SNE")+
  theme_classic()+
  theme(axis.title=element_text(size=21), axis.text= element_text(size=16, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="top", legend.key.size = unit(1.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))
dev.off()

#Expression of PGCLC genes within single-cell epiblast
PGCLC_genes <- read.table(paste0(dir,"scRNAseq_PGCLC_differentiation/PGCLC_genes.txt"))
colnames(PGCLC_genes) <- c("ensembl_id", "gene")

scPGC_genes <- as.data.frame(sc_rpkm) %>% dplyr::filter(rownames(sc_rpkm) %in% PGCLC_genes$ensembl_id)
epi <- meta_epi %>% dplyr::filter(meta_epi$sample %in% colnames(scPGC_genes))
scPGC <- scPGC_genes %>% dplyr::select(epi$sample)

#tSNE plot for PGCLC genes only
set.seed(1989)
epi <- meta_epi %>% dplyr::filter(meta_epi$sample %in% colnames(scPGC_genes))
tsne_prep <- data.frame(scPGC_genes) %>% dplyr::select(epi$sample)
d <- dist(t(tsne_prep))
r <- Rtsne(d, is.matrix=T)
tnse_plot <- data.frame(x = r$Y[,1], y = r$Y[,2], State = epi$stage)
pdf("tSNE_PGCLC_genes_Epiblast.pdf", height = 5, width = 5.5)
ggplot(tnse_plot, aes(y, x, colour = State))+geom_point(size=2) +  scale_color_manual(values = c("#0101DF", "#DF0101", "#FF8000"), "Stage") + 
  labs(x = "Component 2: t-SNE", y = "Component 1: t-SNE")+
  theme_classic()+
  theme(axis.title=element_text(size=21), axis.text= element_text(size=16, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="top", legend.key.size = unit(1.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=6)))
dev.off()


#single-cell DNA methylation of the PGCLC enhancers
load(paste0(dir, "Heterogeneity/meta_data_epiblast_scMT.Rdata"))

#Analysis of the PGCLC gene expression and single-cell DNA methylation
Stages <- c("E4.5", "E5.5", "E6.5")
meta_data_enhancer <- data.frame()
for(Stage in Stages){
  enhancer <- meta_data %>%
    filter(stage==Stage) %>%
    arrange(mCpG_local)
  rownames(enhancer) <- enhancer$id
  enhancer$id <- NULL
  meta_data_enhancer <- rbind(meta_data_enhancer, enhancer)
}

scM <- meta_data_enhancer


#single-cell expression of the PGCLC genes (linked to an enhancer)
Enhancer_linked_genes <- read.table(paste0(dir, "Enhancer_definition/PGCLC_enhancer_genes.txt"), header = T)
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
PGCLC_genes_enhancer <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), filters = 'external_gene_name', values = Enhancer_linked_genes$gene, mart = ensembl, useCache = F)
colnames(PGCLC_genes_enhancer) <- c("ensembl_id", "gene")
detach("package:biomaRt", unload=TRUE)


scPGC_enhancer <- as.data.frame(sc_rpkm) %>% dplyr::filter(rownames(sc_rpkm) %in% PGCLC_genes_enhancer$ensembl_id)
epi <- meta_epi %>% dplyr::filter(meta_epi$sample %in% colnames(scPGC_enhancer))
scPGC_genes_enhancer <- scPGC_enhancer %>% dplyr::select(epi$sample)


Stages <- c("E4.5", "E5.5", "E6.5")
sc_Expression <- data.frame()
for(i in Stages){
  tmp <- data.frame()
  sub <- subset(meta_data, meta_data$stage==i)
  tmp <- scPGC_genes_enhancer[,match(sub$sample, colnames(scPGC_genes_enhancer))]
  tmp <- data.frame(Expression=colMeans(tmp, na.rm = T), Stage=sub$stage)
  sc_Expression<- rbind(sc_Expression, tmp)}

#remove outliers that are most likely artfacts
sc_Expression <- subset(sc_Expression, sc_Expression$Expression>0 & sc_Expression$Expression<40)

#combine the single-cell expression and DNA methylation data
scM_sorted <- scM[match(rownames(sc_Expression), scM$sample),]
scMT <- cbind(scM_sorted, sc_Expression)
scMT$Stage <- NULL

pdf("single-cell expression and DNA methylation and PGCLC genes and enhancers")
scMT_plot <- ggplot(scMT, aes(mCpG_local, log2(Expression+1), color=stage))+geom_point()+
  theme_classic()+theme(legend.position="top", axis.text.x= element_text(color = "black", size=12), axis.text.y= element_text(color = "black", size=12), axis.line = element_line(colour = "white"), axis.title=element_text(size=13))+
  labs(x="CpG methylation", y="Expression - log2 (fpkm+1)")+scale_color_manual(values = c("#0101DF", "#DF0101", "#FF8000"))+
  scale_x_continuous(labels = scales::percent_format(accuracy = 1))
print(scMT_plot)
dev.off()