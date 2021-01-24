library("ggplot2")
library("dplyr")
library("ggpubr")
library("Hmisc")

dir <- "PATH/Germline_competence/"

#Load complete UMI count matrix and define cell stages
Expression <- read.csv(paste0(dir,"/scRNAseq/UMI_counts.csv"), row.names = 1)

#Transfer a correlation matrix into a data.frame
flattenCorrelationMatrix <- function(cmatrix, pvalue) {
  ut <- upper.tri(cmatrix)
  data.frame(
    row = rownames(cmatrix)[row(cmatrix)[ut]],
    column = rownames(cmatrix)[col(cmatrix)[ut]],
    cor  =(cmatrix)[ut],
    p = pvalue[ut]
  )
}

#Load the meta data for the scRNAseq data set
load(paste0(dir,"/scRNAseq/meta_scRNAseq_PGCLC.Rdata"))

#Calculation of the transcriptional noise for each stage by:
#Selection of the stage
#Calculation the variance for all genes of a stage and select the top 500 (top_n)
#Spearman correlation of all highly varibale genes
#Transformation of the Spearman correlation into the transcriptional noise by sqrt((1-cor)/2)
States <- c("ESC", "d1 EpiLC", "d2 EpiLC", "d3 EpiLC", "EpiSC", "d2 PGCLC", "d4 PGCLC")
noise <- data.frame()
for(i in States){
  print(i)
  all <- data.frame(Expression)  %>% dplyr::select((subset(meta_seurat_PGCLC, meta_seurat_PGCLC$Sample==i)$id))
  all$var <- apply(all,1,var)
  all <- top_n(all, 500, var) %>% dplyr::select(-var)
  corr <- rcorr(as.matrix(all), type = c("spearman"))
  c <- flattenCorrelationMatrix(corr$r, corr$P) %>% dplyr::mutate(State=i, cor=sqrt((1-cor)/2))
  noise <- rbind(noise, c)
}

noise$State <- factor(noise$State, levels=c("ESC", "d1 EpiLC", "d2 EpiLC", "d3 EpiLC", "EpiSC", "d2 PGCLC", "d4 PGCLC"))

pdf("Transcriptional Noise.pdf", height=5, width = 6)
ggplot(noise, aes(State, cor, fill=State))+geom_violin()+geom_boxplot(width=0.3, outlier.shape = NA, fill="white")+
  labs(x = "", y = "Transcriptional noise")+
  theme_classic()+scale_fill_manual(values = c("#0101DF", "#8904B1", "#DF0101", "#B40404", "#FF8000", "#01DF01", "#088A08"))+
  theme(axis.title=element_text(size=21), axis.text= element_text(size=18, color = "black"), axis.text.x = element_text(angle = 45, hjust = 1, size = 21), axis.line = element_blank(), legend.position="")+
  stat_compare_means(ref.group = "d2 EpiLC", label = "p.signif", label.y = 0.55, size=14, symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("*", "*C","*B","*A", "ns")))+
  ylim(0.17,0.6)
dev.off()