library(Seurat)
library(ggplot2)
library(Hmisc)
library(ggpubr)
library(tidyverse)

dir <- "PATH/Germline_competence/"

#Load cellranger files and remove cells with low transcript abundancy
R1_WT <- Read10X(data.dir = paste0(dir,"scRNAseq_dCD/WT/filtered_feature_bc_matrix/"))
R1_WT <- CreateSeuratObject(counts = R1_WT, project = "PGCLC")
R1_WT <- AddMetaData(R1_WT, "WT", col.name = "Cell")
Mll34_dCD <- Read10X(data.dir = paste0(dir, "scRNAseq_dCD/dCD/filtered_feature_bc_matrix/"))
Mll34_dCD <- CreateSeuratObject(counts = Mll34_dCD, project = "PGCLC")
Mll34_dCD <- AddMetaData(Mll34_dCD, "dCD", col.name = "Cell")
scRNA_d4EB <- merge(R1_WT, y = Mll34_dCD, add.cell.ids = c("WT", "Mll3/4 dCD"), project = "PGCLC")

mt.genes <- rownames(scRNA_d4EB)[grep("^mt-",rownames(scRNA_d4EB))]
C<-GetAssayData(object = scRNA_d4EB, slot = "counts")
percent.mito <- apply(C[mt.genes,], 2, sum)/Matrix::colSums(C)*100
scRNA_d4EB <- AddMetaData(scRNA_d4EB, percent.mito, col.name = "percent.mito")

rb.genes <- rownames(scRNA_d4EB)[grep("^Rp[sl]",rownames(scRNA_d4EB))]
percent.ribo <- apply(C[rb.genes,], 2, sum)/Matrix::colSums(C)*100
scRNA_d4EB <- AddMetaData(scRNA_d4EB, percent.ribo, col.name = "percent.ribo")

pdf("/Users/patrick/Desktop/CMMC/Mll34_CD/scRNAseq/pdf/QC_filter_scRNAseq_dCD.pdf")
print(FeatureScatter(scRNA_d4EB, feature1="percent.ribo", feature2="percent.mito", group.by = "Cell"), )
print(FeatureScatter(scRNA_d4EB, feature1="nCount_RNA", feature2="percent.mito", group.by = "Cell"))
print(FeatureScatter(scRNA_d4EB, feature1="nFeature_RNA", feature2="percent.mito", group.by = "Cell"))
dev.off()

selected <- WhichCells(scRNA_d4EB, expression = nFeature_RNA > 2000)
data.filt <- subset(scRNA_d4EB, cells = selected)
scRNA_d4EB <- data.filt

#Apply standard tSNE, UMAP and clustering
set.seed(1989)
scRNA_d4EB <- ScaleData(scRNA_d4EB)
scRNA_d4EB <- FindVariableFeatures(object = scRNA_d4EB)
scRNA_d4EB <- RunPCA(object = scRNA_d4EB)
scRNA_d4EB <- RunTSNE(object = scRNA_d4EB, dims = 1:10)
scRNA_d4EB <- RunUMAP(object = scRNA_d4EB, dims = 1:10)

#keep the Seurat object without clustering for the filtering of PGCLC later
no_cluster_scRNA_d4EB <- scRNA_d4EB

#Find clusters in the data set
scRNA_d4EB <- FindNeighbors(object = scRNA_d4EB)
scRNA_d4EB <- FindClusters(object = scRNA_d4EB, resolution = 0.3)

#Cluster identification
load(paste0(dir,"scRNAseq_PGCLC_differentiation/tissues_E8_25.Rdata"))
all <- data.frame()
for(cluster in levels(scRNA_d4EB@meta.data$seurat_clusters)){
  cluster_counts <- as.data.frame(GetAssayData(subset(scRNA_d4EB, idents = cluster, slot="scale.data")))
  print(ncol(cluster_counts))
  cluster_counts$mean <- rowMeans(cluster_counts)
  
  Cluster_subset <- data.frame()
  for(i in 1:length(tissues)){
    summary <- cluster_counts %>% mutate(Tissue = names(tissues[i]), Cluster=as.character(cluster), X=rownames(cluster_counts))
    Tissue_subset <- inner_join(summary, tissues[[i]], by="X") %>% dplyr::select(Tissue, Cluster, mean, X)
    Cluster_subset <- rbind(Cluster_subset, Tissue_subset)}
  all <- rbind(all, Cluster_subset)
}

#remove the tissue with low expression in all clusters and visualize the remiaing ones
all <- subset(all, all$Tissue!="SomaticMesoderm")
all <- subset(all, all$Tissue!="PreSomaticMesoderm")
all <- subset(all, all$Tissue!="PharyngealMesoderm")
all <- subset(all, all$Tissue!="mixedMesoderm")
all <- subset(all, all$Tissue!="midHindBrain")
all <- subset(all, all$Tissue!="midHindGut")
all <- subset(all, all$Tissue!="Notochord")
all <- subset(all, all$Tissue!="NeuralTube")
all <- subset(all, all$Tissue!="MesodermProgenitors")
all <- subset(all, all$Tissue!="Cardiac")
ggplot(all, aes(Tissue, log2(mean+1), fill=Cluster))+geom_boxplot(outlier.shape = NA)+
  theme_classic()

#Additional cluster, like the 2-cell-LC cluster can be assessed by marker gene expression
markers <- FindAllMarkers(scRNA_d4EB, only.pos = TRUE)
subset(markers, markers$cluster=="7")

#Transcriptional correlation of the most 1000 variable genes in the main cluster
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

transcriptional_correlation <- data.frame()
for(cluster in 0:5){
  cluster_counts <- as.data.frame(GetAssayData(subset(scRNA_d4EB, idents = cluster)))
  print(ncol(cluster_counts))
  cluster_counts$var <- apply(cluster_counts,1,var)
  cluster_counts <- top_n(cluster_counts, 1000, var) %>% dplyr::select(-var)
  corr <- rcorr(as.matrix(cluster_counts), type = c("spearman"))
  c <- flattenCorrMatrix(corr$r) %>% dplyr::mutate(Cluster=as.character(cluster))
  transcriptional_correlation  <- rbind(transcriptional_correlation , c)}


pdf("Transcriptional_correlation.pdf", height = 4, width = 4)
transcriptional_correlation $Cluster <- factor(transcriptional_correlation $Cluster, levels= c(2,0,4,1,3,5))
plot <- ggplot(transcriptional_correlation, aes(cor, Cluster, fill=Cluster))+geom_violin()+geom_boxplot(width=0.3, outlier.shape = NA, fill="white")+
  labs(x = "Transcriptional correlation", y="")+
  theme_classic()+
  theme(axis.title=element_text(size=12), axis.text= element_text(size=12, color = "black"), axis.text.x = element_text(size = 12), axis.line = element_blank(), legend.position="")+
  scale_y_discrete(labels=c("0" = "undefined\n dCD 1", "1" = "ExM-LC 2", "2" = "undefined\n dCD 2", "3" = "ExE-LC", "4"="ExM-LC 1", "5"="PGCLC"))+scale_fill_manual(values=c("darkgray", "darkgray", "#04B404", "#04B404", "#04B404", "#04B404"))
print(plot)
dev.off()


#Graphical output of the analysis, like UMAP plot of the cell types
pdf("scRNAseq_d4EB_UMAP_clusters.pdf", height = 6, width = 6)
UMAP_plot <- DimPlot(object = scRNA_d4EB, reduction = "umap", group.by ="Cell", pt.size = 0.1, cols = c("#0101DF", "#DF0101"))+theme_classic()+
  theme(axis.title=element_text(size=21), axis.text= element_text(size=18, color = "black"), axis.line = element_blank(), legend.title = element_blank(), legend.text=element_text(size=21), legend.position="top")+labs(x="UMAP 1", y="UMAP 2", title = "")
print(UMAP_plot)

#Assigning the cluster identity
new.cluster.ids <- c("undefined\n dCD 1", "ExM-LC 2", "undefined\n dCD 2", "ExE-LC","ExM-LC 1", "PGCLC", "ExEndoderm/\n Gut-like", "2-cell like", "Endothelial")
names(new.cluster.ids) <- levels(scRNA_d4EB)
scRNA_d4EB <- RenameIdents(scRNA_d4EB, new.cluster.ids)
plot_cluster <- DimPlot(object = scRNA_d4EB, reduction = "umap", label = TRUE, pt.size = 0.3)+theme_classic()+
  theme(axis.title=element_text(size=16), axis.text= element_text(size=12, color = "black"), axis.line = element_blank(), legend.title = element_blank(), legend.text=element_text(size=21), legend.position="")+labs(x="UMAP 1", y="UMAP 2", title = "")
print(plot_cluster)

#Percentage of the cells in each cluster per cell line
stats <- data.frame(table(scRNA_d4EB$seurat_clusters, scRNA_d4EB$Cell))
stats$norm_Freq <- c(stats$Freq[1:9]/1699*100, stats$Freq[10:18]/1416*100)
stats$Var2 <- factor(stats$Var2, levels = c("WT", "dCD"))
plot <- ggplot(stats, aes(x=Var2, y=Var1, size=norm_Freq, color=Var2))+geom_point()+theme_classic()+
  theme(axis.title=element_text(size=13), axis.text= element_text(size=12, color = "black"), axis.line = element_blank(), legend.title = element_blank(), legend.text=element_text(size=12), legend.position="top")+
  scale_color_manual(values=c("gray", "#DF0101"))+scale_y_discrete(labels=c("0" = "undefined\n dCD 1", "1" = "ExM-LC 2", "2" = "undefined\n dCD 2", "3" = "ExE-LC", "4"="ExM-LC 1", "5"="PGCLC", "6"="ExEndoderm/\n Gut-like", "7"="2-cell like", "8"="Endothelial"))+
  labs(x="", y="")+guides(color=FALSE)
print(plot)
dev.off()


#Filter for putative PGCLC cells by Prmd1+ Dppa3+ and Klf4- (pluripotency marker)
pdf("Dppa3_Prdm1_Klf4_plots.pdf", height = 4, width = 4)
FeaturePlot(object = scRNA_d4EB, features = c("Dppa3", "Tfap2c", "Prdm1", "Klf4"), pt.size = 0.1, cols = c("orange", "blue"), ncol = 1)

Dpp_Klf4 <- subset(no_cluster_scRNA_d4EB, subset = (Dppa3 > 0.1|Prdm1>1)&Klf4<1)
subset_Dpp_Klf4 <- DimPlot(object = Dpp_Klf4, pt.size = 0.2, group.by = "Cell", cols = c("#0101DF", "#DF0101"))+theme_classic()+
  theme(axis.title=element_text(size=21), axis.text= element_text(size=18, color = "black"), axis.line = element_blank(), legend.title = element_blank(), legend.text=element_text(size=21), legend.position="top")+labs(x="UMAP 1", y="UMAP 2", title = "")
print(subset_Dpp_Klf4)
dev.off()

#Expression of PGCLC genes within the Dppa3/Prdm1+ and Klf4- cells
Enhancer_linked_genes <- read.table(paste0(dir, "Enhancer_definition/PGCLC_enhancer_genes.txt"), header = T)
Expression_Enhancer_linked_genes <- data.frame(AverageExpression(Dpp_Klf4, features = Enhancer_linked_genes$gene, slot = "scale.data", group.by = "Cell"))
colnames(Expression_Enhancer_linked_genes) <- c("dCD", "WT")

PGCLC_genes <- read.table(paste0(dir, "scRNAseq_PGCLC_differentiation/PGCLC_genes.txt"))
colnames(PGCLC_genes) <- c("Ensembl_id", "gene")
PGCLC_genes_TSS_only <- anti_join(data.frame(PGCLC_genes), Enhancer_linked_genes, by="gene")
PGCLC_genes_TSS_only <- subset(PGCLC_genes_TSS_only, PGCLC_genes_TSS_only$gene!="Rps23rg1"&PGCLC_genes_TSS_only$gene!="Hrob"&PGCLC_genes_TSS_only$gene!="Septin1")
Expression_PGCLC_genes_TSS <- data.frame(AverageExpression(Dpp_Klf4, features = PGCLC_genes_TSS_only$gene, slot = "scale.data", group.by = "Cell"))
colnames(Expression_PGCLC_genes_TSS) <- c("dCD", "WT")

linked_genes <- stack(Expression_Enhancer_linked_genes)
TSS_genes <- stack(Expression_PGCLC_genes_TSS)
linked_genes$Group <- "PGCLC genes \n linked to enhancers"
TSS_genes$Group <- "PGCLC genes \n w/o enhancers"

comparision <- rbind(linked_genes, TSS_genes)

pdf("Comparision of the PGCLC gene expression of Dppa3/Prdm1 cells linked to enhancers")
comparision$ind <- factor(comparision$ind, levels = c("dCD", "WT"))
print(ggplot(comparision, aes(y=values, x=ind, fill=ind))+geom_boxplot()+labs(x="", y="scaled expression"))+facet_wrap(~Group)+stat_compare_means(comparisons =  list(c("WT","dCD" )), method = "t.test", label = "p.format", paired=T, size=5, label.y = 2)+theme_classic()+theme(strip.text.x = element_text(size = 18))+theme(strip.background = element_blank())+
  theme(axis.title=element_text(size=16), axis.text= element_text(size=14, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="", legend.key.size = unit(1, "cm"))+annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+scale_fill_manual(values=c("#0101DF", "#DF0101"))
dev.off()