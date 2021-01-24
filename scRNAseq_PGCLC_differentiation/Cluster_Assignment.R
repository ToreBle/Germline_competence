library("ggplot2")
library("tidyverse")
library("ggpubr")
  
  
  dir <- "PATH/Germline_competence/"
  
  #All cell types from E8.25 embryos defined by Ibarra-Soria et al. 2018 
  load(paste0(dir,"scRNAseq_PGCLC_differentiation/tissues_E8_25.Rdata"))
  
  #Get the expression of all tissue genes within each cluster of the d4EB
  pdf(file="Cluster_analysis_d4EB.pdf", height=5, width = 10)
  for(c in 1:cluster){
    #read cluster
    print(c)
    cluster_name <- paste("cluster", c, sep = "")
    cluster_i <- read.csv(paste(dir, "scRNAseq_PGCLC_differentiation/d4_fullCluster_", c, ".csv", sep=""))
    cluster_i <- cluster_i %>% mutate(Group=rep(cluster_name, nrow(ci))) %>% dplyr::select(X, Group, row.mean)
    
    tmp <- data.frame()
    for(i in 1:length(tissues)){
      cluster_i <- cluster_i %>% mutate(Tissue = paste(names(tissues[i]), " (n=", nrow(tissues[[i]]),")", sep=""), Cluster=rep(cluster_name, nrow(ci)))
      cluster_tissue_expression <- inner_join(cluster_i, tissues[[i]], by="X")
      tmp <- rbind(tmp,cluster_tissue_expression )
    }
    
    tmp<- tmp %>% mutate("Gene"=X, mean=log2(row.mean+1))%>% dplyr::select(Gene, Tissue, Cluster, Group, mean)
    #print(ggboxplot(tmp, x = "Tissue", y = "mean", color = "Group", palette = "jco", outlier.shape = NA)+labs(Title=cluster_name, x =" ", y="log2(mean Expression+1)")+coord_flip())
    assign(cluster_name, tmp)
  }
  

#combine all clusters and discard tissues with low expression in all clusters
d4 <- rbind(cluster1, cluster2, cluster3, cluster4)

d4 <- subset(d4, d4$Tissue!="SomaticMesoderm (n=12)")
d4 <- subset(d4, d4$Tissue!="PreSomaticMesoderm (n=16)")
d4 <- subset(d4, d4$Tissue!="PharyngealMesoderm (n=14)")
d4 <- subset(d4, d4$Tissue!="mixedMesoderm (n=6)")
d4 <- subset(d4, d4$Tissue!="midHindBrain (n=6)")
d4 <- subset(d4, d4$Tissue!="midHindGut (n=7)")

pdf("d4EB_clusters_summary.pdf", width = 7, height = 9)
ggplot(d4, aes(log2(mean+1), Tissue, fill=Cluster))+geom_boxplot(outlier.shape = NA)+
  theme_classic()+theme(legend.position="top", axis.text.y= element_text(color = c("black", "black", "#088A08", "black", "#088A85", "darkgray", "black", "#01DFA5", rep("black",8)), size=10), axis.text.x= element_text(color = "black", size=12), axis.line = element_line(colour = "white"), axis.title=element_text(size=14))+
  labs(x="log2 (UMI counts+1)", y="Tissues-specific gene expression \n(Ibarra-Soria et al. 2018)")+
  scale_fill_manual(name = "Cluster", labels = c("Endothelial", "PGCLC", "ExE-LC", "ExM-LC"), values = c("#088A85", "#088A08", "darkgray", "#01DFA5"))
dev.off()


#Tile plots
c1 <- cluster1 %>% group_by(Tissue,Group) %>% summarise(m=mean(mean), sd=var(mean))
c1$Cluster <- "Endothelial"
c2 <- cluster2 %>% group_by(Tissue,Group) %>% summarise(m=mean(mean), sd=var(mean))
c2$Cluster <- "PGCLC"
c3 <- cluster3 %>% group_by(Tissue,Group) %>% summarise(m=mean(mean), sd=var(mean))
c3$Cluster <- "ExE-LC"
c4 <- cluster4 %>% group_by(Tissue,Group) %>% summarise(m=mean(mean), sd=var(mean))
c4$Cluster <- "ExM-LC"
d4_clusters <- rbind(c1,c2,c3,c4)
d4_clusters <- subset(d4_clusters, d4_clusters$Tissue!="SomaticMesoderm (n=12)")
d4_clusters <- subset(d4_clusters, d4_clusters$Tissue!="PreSomaticMesoderm (n=16)")
d4_clusters <- subset(d4_clusters, d4_clusters$Tissue!="PharyngealMesoderm (n=14)")
d4_clusters <- subset(d4_clusters, d4_clusters$Tissue!="mixedMesoderm (n=6)")

d4_clusters$Cluster <- factor(d4_clusters$Cluster, levels=c("Endothelial", "PGCLC", "ExE-LC", "ExM-LC"))
pdf("/Users/patrick/Desktop/CMMC/singleCell_RNAseq/Cell Classification/d4_clusters_summary.pdf", width = 6, height = 8)
ggplot(d4_clusters, aes(Cluster, Tissue, size=(m), color=Cluster))+geom_point()+
  theme_classic()+theme(legend.position="top", axis.text.x= element_text(color = "black", size=12, angle=45, hjust=1), axis.text.y= element_text(color = "black", size=10), axis.line = element_line(colour = "white"), axis.title=element_text(size=14))+
  labs(size="UMI counts", x="Clusters", y="Tissues-specific gene expression \n(Ibarra-Soria et al. 2018)")+guides(color = FALSE)+scale_color_manual(values = c("#088A85", "#088A08", "darkgray", "#01DFA5"))+scale_size_continuous(limits =  c(0,3))
dev.off()



