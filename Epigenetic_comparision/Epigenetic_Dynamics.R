library(tidyverse)
library(patchwork)
library(coin)
library(rstatix)
library(ggplot2)

dir <- "PATH/Germline_competence/"

#Load all quantifications for ATAC-, ChIP- and WGB-seq for:
#EpiLC enhancer
#EpiSC enhancer
#PGCLC enhancer
#PGCLC TSS

Epigenetic_data <- read.table(paste0(dir, "Epigenetic_comparision/Epigenetic_data_summary_1kb.txt"), h=T)

#omit the coordinates and stack for visualization
Epigenetic_data <- Epigenetic_data  %>% dplyr::select(enhancer:H3K27ac_EpiSC)
stacked_epigenetic_data <- melt(Epigenetic_data)
colnames(stacked_epigenetic_data) <- c("Enhancer", "Method_Stage", "RPGC")
overview <- stacked_epigenetic_data %>% group_by(Enhancer, Method_Stage) %>% summarise(mean_RPGC=mean(RPGC))
overview <- overview %>% separate(Method_Stage, c("Method", "Stage"))
ChIPs <- unique(overview$Method)

#split by Method (individual ChIP-seq experiment, ATAC-seq and Whole-genome-bisulfite sequencing)
for(c in 1:10){
  tmp <- subset(overview, overview$Method==ChIPs[c])
  tmp$Stage <- factor(tmp$Stage, levels=c("ESC", "EpiLC", "EpiSC"))
  tmp$Enhancer <- factor(tmp$Enhancer, levels=c("PGCLC_TSS", "PGCLC_enhancer", "EpiSC_enhancer", "EpiLC_enhancer"))
  if(c==1){
    ATAC_plot <- ggplot(tmp, aes(x=Stage, y=Enhancer, fill=log10(mean_RPGC+1)))+geom_tile(stat="identity", width=1, height=.95)+coord_equal()+
      theme_classic()+scale_fill_gradient2(low="#0404B4", midpoint = 1.05, high="#DF0101")+
      theme(plot.title = element_text(size=16, color="black", hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.text= element_text(size=11, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 9), legend.text=element_text(size=6), legend.position="top", legend.key.size = unit(0.4, "cm"))+
      labs(x="",y="", fill="log10 (RPGC+1)", title=ChIPs[c])
  }
  if(c==2){
    H3K4me1_plot <- ggplot(tmp, aes(x=Stage, y=Enhancer, fill=log10(mean_RPGC+1)))+geom_tile(stat="identity", width=1, height=.95)+coord_equal()+
      theme_classic()+scale_fill_gradient2(low="#0404B4", midpoint = 0.85, high="#DF0101")+
      theme(plot.title = element_text(size=16, color="black", hjust = 0.5), axis.text.y = element_text(face = c('bold', 'bold', 'plain', 'plain'), colour = c("#0101DF", "#0101DF", "black", "black"), size=c(12,12,11,11)), axis.text= element_text(size=11, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 9), legend.text=element_text(size=6), legend.position="top", legend.key.size = unit(0.4, "cm"))+
      labs(x="",y="", fill="log10 (RPGC+1)", title=ChIPs[c])+scale_y_discrete(labels=c("PGCLC \nTSS", "PGCLC \nenhancer", "EpiSC \nenhancer", "EpiLC \nenhancer"))
  }
  if(c==3){
    H3K4me2_plot <- ggplot(tmp, aes(x=Stage, y=Enhancer, fill=log10(mean_RPGC+1)))+geom_tile(stat="identity", width=1, height=.95)+coord_equal()+
      theme_classic()+scale_fill_gradient2(low="#0404B4", midpoint = 1.05, high="#DF0101")+
      theme(plot.title = element_text(size=16, color="black", hjust = 0.5), axis.text= element_text(size=11, color = "black"), axis.line = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.title = element_text(size = 9), legend.text=element_text(size=6), legend.position="top", legend.key.size = unit(0.4, "cm"))+
      labs(x="",y="", fill="log10 (RPGC+1)", title=ChIPs[c])
  }
  if(c==4){
    H3K4me3_plot <- ggplot(tmp, aes(x=Stage, y=Enhancer, fill=log10(mean_RPGC+1)))+geom_tile(stat="identity", width=1, height=.95)+coord_equal()+
      theme_classic()+scale_fill_gradient2(low="#0404B4", midpoint = 1.05, high="#DF0101")+
      theme(plot.title = element_text(size=16, color="black", hjust = 0.5), axis.text.y = element_text(face = c('bold', 'bold', 'plain', 'plain'), colour = c("#0101DF", "#0101DF", "black", "black"), size=c(12,12,11,11)), axis.text= element_text(size=11, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 9), legend.text=element_text(size=6), legend.position="top", legend.key.size = unit(0.4, "cm"))+
      labs(x="",y="", fill="log10 (RPGC+1)", title=ChIPs[c])+scale_y_discrete(labels=c("PGCLC \nTSS", "PGCLC \nenhancer", "EpiSC \nenhancer", "EpiLC \nenhancer"))
    
  }
  if(c==5){
    H3K9me2_plot <- ggplot(tmp, aes(x=Stage, y=Enhancer, fill=log10(mean_RPGC+1)))+geom_tile(stat="identity", width=1, height=.95)+coord_equal()+
      theme_classic()+scale_fill_gradient2(low="#0404B4", midpoint = 0.2, high="#DF0101")+
      theme(plot.title = element_text(size=18, color="black", hjust = 0.5), axis.text= element_text(size=11, color = "black"), axis.line = element_blank(),  axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.title = element_text(size = 9), legend.text=element_text(size=6), legend.position="top", legend.key.size = unit(0.4, "cm"))+
      labs(x="",y="", fill="log10 (RPGC+1)", title=ChIPs[c])
  }
  if(c==6){
    H3K9me3_plot <- ggplot(tmp, aes(x=Stage, y=Enhancer, fill=log10(mean_RPGC+1)))+geom_tile(stat="identity", width=1, height=.95)+coord_equal()+
      theme_classic()+scale_fill_gradient2(low="#0404B4", midpoint = 0.35, high="#DF0101")+
      theme(plot.title = element_text(size=18, color="black", hjust = 0.5), axis.text= element_text(size=11, color = "black"), axis.line = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.title = element_text(size = 9), legend.text=element_text(size=6), legend.position="top", legend.key.size = unit(0.4, "cm"))+
      labs(x="",y="", fill="log10 (RPGC+1)", title=ChIPs[c])
  }
  if(c==7){
    H3K27me2_plot <- ggplot(tmp, aes(x=Stage, y=Enhancer, fill=log10(mean_RPGC+1)))+geom_tile(stat="identity", width=1, height=.95)+coord_equal()+
      theme_classic()+scale_fill_gradient2(low="#0404B4", midpoint = 0.28, high="#DF0101")+
      theme(plot.title = element_text(size=18, color="black", hjust = 0.5), axis.text= element_text(size=11, color = "black"), axis.line = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.title = element_text(size = 9), legend.text=element_text(size=6), legend.position="top", legend.key.size = unit(0.4, "cm"))+
      labs(x="",y="", fill="log10 (RPGC+1)", title=ChIPs[c])
  }
  if(c==8){
    H3K27me3_plot <- ggplot(tmp, aes(x=Stage, y=Enhancer, fill=log10(mean_RPGC+1)))+geom_tile(stat="identity", width=1, height=.95)+coord_equal()+
      theme_classic()+scale_fill_gradient2(low="#0404B4", midpoint = 0.42, high="#DF0101")+
      theme(plot.title = element_text(size=16, color="black", hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.text= element_text(size=11, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 9), legend.text=element_text(size=6), legend.position="top", legend.key.size = unit(0.4, "cm"))+
      labs(x="",y="", fill="log10 (RPGC+1)", title=ChIPs[c])
  }
  if(c==9){
    mCpG_plot <- ggplot(tmp, aes(x=Stage, y=Enhancer, fill=mean_RPGC))+geom_tile(stat="identity", width=1, height=.95)+coord_equal()+
      theme_classic()+scale_fill_gradient2(low="#0404B4", midpoint = 48, high="#DF0101")+
      theme(plot.title = element_text(size=18, color="black", hjust = 0.5), axis.text= element_text(size=11, color = "black"), axis.line = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.title = element_text(size = 9), legend.text=element_text(size=6), legend.position="top", legend.key.size = unit(0.4, "cm"))+
      labs(x="",y="", fill="% mCpG", title=ChIPs[c])
  }
  if(c==10){
    H3K27ac_plot <- ggplot(tmp, aes(x=Stage, y=Enhancer, fill=log10(mean_RPGC+1)))+geom_tile(stat="identity", width=1, height=.95)+coord_equal()+
      theme_classic()+scale_fill_gradient2(low="#0404B4", midpoint = 1.1, high="#DF0101")+
      theme(plot.title = element_text(size=16, color="black", hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.text= element_text(size=11, color = "black"), axis.line = element_blank(), legend.title = element_text(size = 9), legend.text=element_text(size=6), legend.position="top", legend.key.size = unit(0.4, "cm"))+
      labs(x="",y="", fill="log10 (RPGC+1)", title=ChIPs[c])  }
  
}

#pdf("tile_plots_EpiLC_EpiSC_CHIP.pdf", height = 4, width = 16)
print(H3K4me1_plot+
      H3K4me2_plot+
      ATAC_plot+
      plot_layout(nrow = 1, byrow = FALSE))

print(H3K27ac_plot+
      H3K27me3_plot+
      H3K9me3_plot+
      plot_layout(nrow = 1, byrow = FALSE))

print(H3K4me3_plot+
      H3K9me2_plot+
      H3K27me2_plot+
      plot_layout(nrow = 1, byrow = FALSE))
dev.off()


#Subset for the PGCLC enhancers and all data from EpiLC and EpiSC only to determine the differences between
#germline competent EpiLC and non-germline competent EpiSC
Epigenetics_PGCLC_enhancer <- subset(Epigenetic_data, Epigenetic_data$enhancer=="PGCLC_enhancer")

d <- data.frame()
r <- data.frame()
EpiLC <- list(ATAC=Epigenetics_PGCLC_enhancer$ATAC_EpiLC, H3K27ac=Epigenetics_PGCLC_enhancer$H3K27ac_EpiLC, H3K4me1=Epigenetics_PGCLC_enhancer$H3K4me1_EpiLC, H3K4me2=Epigenetics_PGCLC_enhancer$H3K4me2_EpiLC, H3K4me3=Epigenetics_PGCLC_enhancer$H3K4me3_EpiLC, H3K9me2=Epigenetics_PGCLC_enhancer$H3K9me2_EpiLC, H3K9me3=Epigenetics_PGCLC_enhancer$H3K9me3_EpiLC, H3K27me2=Epigenetics_PGCLC_enhancer$H3K27me2_EpiLC, H3K27me3=Epigenetics_PGCLC_enhancer$H3K27me3_EpiLC, mCpG=Epigenetics_PGCLC_enhancer$mCpG_EpiLC)
EpiSC <- list(ATAC=Epigenetics_PGCLC_enhancer$ATAC_EpiSC, H3K27ac=Epigenetics_PGCLC_enhancer$H3K27ac_EpiSC, H3K4me1=Epigenetics_PGCLC_enhancer$H3K4me1_EpiSC, H3K4me2=Epigenetics_PGCLC_enhancer$H3K4me2_EpiSC, H3K4me3=Epigenetics_PGCLC_enhancer$H3K4me3_EpiSC, H3K9me2=Epigenetics_PGCLC_enhancer$H3K9me2_EpiSC, H3K9me3=Epigenetics_PGCLC_enhancer$H3K9me3_EpiSC, H3K27me2=Epigenetics_PGCLC_enhancer$H3K27me2_EpiSC, H3K27me3=Epigenetics_PGCLC_enhancer$H3K27me3_EpiSC, mCpG=Epigenetics_PGCLC_enhancer$mCpG_EpiSC)

#Effect size for each ChIP-seq, ATAC-seq and genome-wide CpG methylation
for(i in 1:length(EpiLC)){
  print(names(EpiLC[i]))
  LC <- data.frame(RPGC=EpiLC[[i]], State="EpiLC")
  SC <- data.frame(RPGC=EpiSC[[i]], State="EpiSC")
  eff <- rbind(LC, SC)
  weff <- wilcox_effsize(eff, RPGC~State, paired = T, ci=T, ci.type = "basic")
  weff$ChIP <- names(EpiLC[i])
  qnorm(wilcox.test(eff$RPGC~eff$State, paired = T)$p.value/2)
  print(mean(LC$RPGC/SC$RPGC, na.rm=T))
  if(mean(LC$RPGC/SC$RPGC, na.rm=T)<1){
    weff <- weff %>% mutate(effsize=effsize*(-1), conf.low=conf.low*(-1), conf.high=conf.high*(-1))}
  if(names(EpiLC[i])=="ATAC"|names(EpiLC[i])=="H3K4me1"|names(EpiLC[i])=="H3K4me2"|names(EpiLC[i])=="H3K4me3"|names(EpiLC[i])=="H3K27ac"){
    weff$Group <- "Active"}
  if(names(EpiLC[i])=="H3K27me2"|names(EpiLC[i])=="H3K27me3"|names(EpiLC[i])=="H3K9me2"|names(EpiLC[i])=="H3K9me3"|names(EpiLC[i])=="mCpG"){
    weff$Group <- "Repressive"}
  r <- rbind(r, weff)
}
r$ChIP <- factor(r$ChIP, levels=(r %>% arrange(effsize))$ChIP)
pdf("Effect size PGCLC enhancer.pdf", height = 4, width = 5)
ggplot(r, aes(y=effsize, x=ChIP, fill=Group))+geom_bar(stat = "identity")+coord_flip()+scale_y_reverse()+theme_classic()+
  theme(axis.title=element_text(size=14), axis.text= element_text(size=12, color = "black"), axis.line = element_blank(), legend.title = element_blank(), legend.text=element_text(size=12), legend.position="top", legend.key.size = unit(0.5, "cm"))+
  scale_fill_manual(values = c("#0101DF", "#DF0101"))+labs(x="", y="Effect size (EpiLC vs. EpiSC)")+geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=.2,position=position_dodge(.9), size=1)+ylim(0.6,-0.85)
dev.off()

