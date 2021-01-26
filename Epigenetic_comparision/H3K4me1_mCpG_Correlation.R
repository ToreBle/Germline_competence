library(ggplot2)
library(dplyr)


dir <- "PATH/Germline_competence/"

H3K4me1_replicates <- read.table(paste0(dir, "Epigenetic_comparision/H3K4me1_replicates_PGCLC_enhancer.txt"), h=T)


H3K4me1_replicates$ratio_rep1 <- H3K4me1_replicates$H3K4me1_EpiLC_rep1/H3K4me1_replicates$H3K4me1_EpiSC_rep1
H3K4me1_replicates$ratio_rep2 <- H3K4me1_replicates$H3K4me1_EpiLC_rep2/H3K4me1_replicates$H3K4me1_EpiSC_rep2
H3K4me1_replicates$ratio_rep3 <- H3K4me1_replicates$H3K4me1_EpiLC_rep3/H3K4me1_replicates$H3K4me1_EpiSC_rep3

cut_off <- 1.2

broad <- H3K4me1_replicates

high <- subset(broad, broad$ratio_rep1>=cut_off)
middle <- subset(broad, broad$ratio_rep1<cut_off & broad$ratio_rep1>(1/cut_off))
low <- subset(broad, broad$ratio_rep1<(1/cut_off))

high$Group_rep1 <- "EpiLC > EpiSC"
middle$Group_rep1 <- "steady"
low$Group_rep1 <- "EpiLC < EpiSC"
broad <- rbind(low,middle,high)

high <- subset(broad, broad$ratio_rep2>=cut_off)
middle <- subset(broad, broad$ratio_rep2<cut_off & broad$ratio_rep2>(1/cut_off))
low <- subset(broad, broad$ratio_rep2<(1/cut_off))

high$Group_rep2 <- "EpiLC > EpiSC"
middle$Group_rep2 <- "steady"
low$Group_rep2 <- "EpiLC < EpiSC"
broad <- rbind(low,middle,high)

high <- subset(broad, broad$ratio_rep3>cut_off)
middle <- subset(broad, broad$ratio_rep3<cut_off & broad$ratio_rep3>(1/cut_off))
low <- subset(broad, broad$ratio_rep3<=(1/cut_off))

high$Group_rep3 <- "EpiLC > EpiSC"
middle$Group_rep3 <- "steady"
low$Group_rep3 <- "EpiLC < EpiSC"
all <- rbind(low,middle,high)

summary_rep1 <- data.frame(all %>% group_by(Group_rep1) %>% summarise(g = n()))
summary_rep1$Group_rep1 <- factor(summary_rep1$Group_rep1, levels = c("EpiLC > EpiSC", "steady", "EpiLC < EpiSC"))

summary_rep2 <- all %>% group_by(Group_rep2) %>% summarise(g = n())
summary_rep2$Group_rep2 <- factor(summary_rep2$Group_rep2, levels = c("EpiLC > EpiSC", "steady", "EpiLC < EpiSC"))

summary_rep3 <- all %>% group_by(Group_rep3) %>% summarise(g = n())
summary_rep3$Group_rep3 <- factor(summary_rep3$Group_rep3, levels = c("EpiLC > EpiSC", "steady", "EpiLC < EpiSC"))

#pdf("Pie_charts_Ratio_H3K4me1.pdf")
p1 <- ggplot(summary_rep1, aes(x="", y=g, fill=Group_rep1))+
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0)+
  theme_classic()+
  theme(axis.title =element_blank(), axis.text= element_blank(), axis.line = element_blank(), legend.title = element_text(size=21), legend.text=element_text(size=14), legend.position="right")+
  labs(x="", y="", fill="H3K4me1 rep1")+scale_fill_manual(values=c("#013ADF", "#088A08", "#40FF00"))

p2 <- ggplot(summary_rep2, aes(x="", y=g, fill=Group_rep2))+
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0)+
  theme_classic()+
  theme(axis.title =element_blank(), axis.text= element_blank(), axis.line = element_blank(), legend.title = element_text(size=21), legend.text=element_text(size=14), legend.position="right")+
  labs(x="", y="", fill="H3K4me1 rep2")+scale_fill_manual(values=c("#013ADF", "#088A08", "#40FF00"))
p3 <- ggplot(summary_rep3, aes(x="", y=g, fill=Group_rep3))+
  geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0)+
  theme_classic()+
  theme(axis.title =element_blank(), axis.text= element_blank(), axis.line = element_blank(), legend.title = element_text(size=21), legend.text=element_text(size=14), legend.position="right")+
  labs(x="", y="", fill="H3K4me1 rep3")+scale_fill_manual(values=c("#013ADF", "#088A08", "#40FF00"))
p1+p2+p3+plot_layout(ncol=1)
dev.off()

####

GroupI <- subset(all, all$Group_rep1=="EpiLC > EpiSC"|all$Group_rep2=="EpiLC > EpiSC"&all$Group_rep3=="EpiLC > EpiSC")
all_granges <- makeGRangesFromDataFrame(all, keep.extra.columns=T)
GroupII <- data.frame(all_granges[-queryHits(findOverlaps(all_granges, GRanges(GroupI), type="any")),])
colnames(GroupII)[1:3] <- c("chr", "start", "end")
GroupII$width <- NULL
GroupII$strand <- NULL
GroupI$Group <- "Group I"
GroupII$Group <- "Group II"

PGCLC_enhancer_grouped <- rbind(GroupI,GroupII) %>% dplyr::select(chr:end,Group)
write.table(PGCLC_enhancer_grouped, row.names = F, col.names = F, quote = FALSE, sep = "\t", file=paste0(dir, "Epigenetic_comparision/PGCLC_enhancer_grouped.bed"))