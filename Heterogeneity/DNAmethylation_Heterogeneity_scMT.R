library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(scales)

dir <- "PATH/Germline_competence"

#Preparation of the meta data from https://github.com/rargelaguet/scnmt_gastrulation and public GEO ids:
meta <- read.table(paste0(dir,"/Heterogeneity/sample_metadata.txt"), header=TRUE)
meta_epi <- subset(meta, meta$lineage10x=="Epiblast" & meta$stage!="E7.5")

GEO_id <- read.table(paste0(dir,"/Heterogeneity/GEO_id_met.txt"), header=TRUE)
meta_epi <- inner_join(meta_epi, GEO_id, by="sample")
rownames(meta_epi) <- meta_epi$id
meta_epi <- meta_epi %>% dplyr::select(sample, stage, id)


#load the Rdata from the Dissimilarity_matrices.R script containing
#Average CpG methylation for the single-cells
#Coverage of the CpGs
load(paste0(dir, "Heterogeneity/DNAmethylation_genome.Rdata"))
meta_global <- meta_data

load(paste0(dir, "Heterogeneity/DNAmethylation_PGCLC_enhancer.Rdata"))
mCpG_comparision <- merge(meta_data, meta_global, by=0)#, all=T)
rownames(mCpG_comparision) <- mCpG_comparision$Row.names
mCpG_comparision$Row.names <- NULL

meta_data <- merge(mCpG_comparision, meta_epi, by=0)
rownames(meta_data) <- meta_data$Row.names
meta_data$Row.names <- NULL
meta_data <- subset(meta_data, meta_data$stage!="E7.5")


#Correlation of CpG methylation at enhancers vs genome-wide and 
#QC for the coverage at PGCLC enhancers
pdf("mCpG_genome_wide_and_cov_correlation.pdf", height = 4, width = 5)
p1 <- ggscatter(meta_data, x = "mCpG_local", y = "mCpG_global", xlab = "PGCLC Enhancer - mCpG/CpG", ylab = "Genome wide - mCpG/CpG", na.rm=T, color = "stage")+
  scale_color_manual(values=c("#0101DF", "#DF0101", "#FF8000"))+theme_classic()+
  theme(axis.title=element_text(size=16), axis.text= element_text(size=14, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="top", legend.key.size = unit(1, "cm"))+
  stat_cor(method = "spearman", size=5)+guides(colour = guide_legend(override.aes = list(size=6)))

p2 <- ggscatter(meta_data, x = "mCpG_local", y = "coverage_local",
                xlab = "PGCLC Enhancer - mCpG/CpG", ylab = "PGCLC enhancer - CpG coverage", na.rm=T, color = "stage")+
  scale_color_manual(values=c("#0101DF", "#DF0101", "#FF8000"))+theme_classic()+
  theme(axis.title=element_text(size=16), axis.text= element_text(size=14, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="top", legend.key.size = unit(1, "cm"))+
  stat_cor(method = "spearman", size=5, label.y=3.80)+guides(colour = guide_legend(override.aes = list(size=6)))
print(p1)
print(p2+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))))
dev.off()



#PGCLC enhancer Heatmap, Dissimilarity and Heterogeneity
#Preparation for the ordering of the dissimility matrix 
#by increasing CpG methylation at PGCLC enhancers in each stage
E4_naive <- meta_data %>%
  filter(stage=="E4.5") %>%
  arrange(mCpG_local)

rownames(E4_naive) <- E4_naive$id

E5_formative <- meta_data %>%
  filter(stage=="E5.5") %>%
  arrange(mCpG_local)
rownames(E5_formative) <- E5_formative$id


E6_primed <- meta_data %>%
  filter(stage=="E6.5") %>%
  arrange(mCpG_local)
rownames(E6_primed) <- E6_primed$id

meta_diss <- rbind(E4_naive, E5_formative, E6_primed)

#Dissimilarity matrix for all PGCLC enhancers
load(paste0(dir, "Heterogeneity/diss_scM_Epiblast_all_PGCLC_enhancer.Rdata"))
diss_ordered <- diss_ordered[match(rownames(meta_diss), rownames(diss_ordered)), ]
diss_ordered <- diss_ordered[,match(rownames(meta_diss), colnames(diss_ordered))]

meta_diss <- meta_diss %>% dplyr::select(mCpG_local, stage)

anno <- HeatmapAnnotation(df = meta_diss, 
                          col = list(stage = c("E4.5" = "#0101DF", "E5.5" = "#DF0101", "E6.5"="#FF8000"),
                                     mCpG_local = colorRamp2(c(0, 0.6, 0.7), c("blue", "#FFFFFF", "orange"))),
                          
                          annotation_legend_param = list(stage = list(nrow=1, title = "Stage", title_position = "topcenter", title_gp = gpar(fontsize = 14),
                                                                      labels_gp = gpar(fontsize = 12), grid_height = unit(0.4, "cm"), space = unit(3, "mm")),
                                                         mCpG_local = list(direction = "horizontal", title = "mCpG/CpG", title_position = "leftcenter", title_gp = gpar(fontsize = 14), 
                                                                           labels_gp = gpar(fontsize = 12), grid_height = unit(0.4, "cm"))),
                          show_annotation_name = FALSE, show_legend =c(TRUE, FALSE))

pdf("Dissimilarity_matrix_PGCLC_enhancer_Epiblast.pdf", height = 5.5, width = 6)
enhancer_diss <- Heatmap(diss_ordered, show_row_names = F, show_column_names = F, cluster_rows = F, cluster_columns = F, top_annotation = anno, name="Dissimilarity", col = colorRamp2(c(0.2,0.4,0.6), c("#0101DF", "white", "#DF0101")), heatmap_legend_param = list(title = "mCpG heterogeneity", title_gp = gpar(fontsize = 14), labels_gp = gpar(fontsize = 12), title_position = "leftcenter-rot", legend_height = unit(4.5, "cm")))#use_raster = T to reduce the image complexity
draw(enhancer_diss, annotation_legend_side = "top")#merge_legend = TRUE)
dev.off()

#Summarizing the dissmiliarty matrix
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

#Dissimilarity matrix for the different enhancers
Enhancer <- c("genome", "GroupI_PGCLC", "GroupII_PGCLC", "EpiLC_enhancer", "EpiSC_enhancer")
Hetero <- data.frame()
for(e in Enhancer){
  load(paste0(dir, "Heterogeneity/diss_scM_Epiblast_", e, ".Rdata"))
  print(ncol(diss_ordered))
  tmp <- flattenCorrMatrix(diss_ordered)
  E4 <- tmp %>% filter(row %in% E4_naive$id & column %in% E4_naive$id) %>% mutate(Stage = "E4.5", Enhancer= e)
  E5 <- tmp %>% filter(row %in% E5_formative$id & column %in% E5_formative$id) %>% mutate(Stage = "E5.5", Enhancer= e)
  E6 <- tmp %>% filter(row %in% E6_primed$id & column %in% E6_primed$id) %>% mutate(Stage = "E6.5", Enhancer= e)
  all <- rbind(E4, E5, E6)
  Hetero <- rbind(Hetero, all)
}

Hetero$Enhancer <- factor(Hetero$Enhancer, levels = c("genome", "EpiLC_enhancer", "EpiSC_enhancer", "GroupI_PGCLC", "GroupII_PGCLC"))

pdf("Epigenetic heterogeneity Comparision.pdf", height = 4, width = 5)
plot <- ggplot(Hetero, aes(x=Enhancer, y=cor, fill=Stage))+geom_violin()+geom_boxplot(position = position_dodge(), width=0.3, fill="white", outlier.shape = NA)+theme_classic()+facet_wrap(~Stage)+
  theme(strip.text.x = element_text(size = 12))+theme(strip.background = element_blank())+
  theme(axis.title=element_text(size=12), axis.text.x = element_text(angle=45, hjust = 1, size=10), axis.text= element_text(size=11, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 18), legend.text=element_text(size=16), legend.position="", legend.key.size = unit(1, "cm"))+
  labs(x="", y="Epigenetic heterogeneity \n (Dissimilarity)")+scale_fill_manual(values = c("#0101DF", "#DF0101", "#FF8000"))+scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,0.45))+
  scale_x_discrete(labels=c("genome" = "Genome", "EpiLC_enhancer" = "EpiLC", "EpiSC_enhancer" = "EpiSC", "GroupI_PGCLC"="PGCLC I", "GroupII_PGCLC"="PGCLC II"))
print(plot)
dev.off()
