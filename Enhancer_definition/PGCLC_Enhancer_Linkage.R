library(dplyr)
library(InteractionSet)
library(IRanges)


dir <- "PATH/Germline_competence/"

#Preparing Genomic Interaction from Atlasi et al, 2019.
source <- read.table(paste0(dir, "Enhancer_definition/hglft_genome_mm10_Atlasi_2019_source.bed"))
colnames(source) <- c("chr", "start", "end")
target <- read.table(paste0(dir, "Enhancer_definition/hglft_genome_mm10_Atlasi_2019_target.bed"))
colnames(target) <- c("chr", "start", "end")

source <- makeGRangesFromDataFrame(source, ignore.strand = T)
target <- makeGRangesFromDataFrame(target, ignore.strand = T)

#Final GInteractions object with all orignial interactions from Atlasi et al. 2019 but with mm10 coordinates
gi <- GInteractions(source, target)

#Define all TSS and filtering for PGCLC TSS
genes <- read.table(paste0(dir, "Enhancer_definition/genes_v86.txt"))
TSS <- genes[c("V3", "V5", "V6", "V4", "V13")]

plus <- subset(TSS, TSS$V4=="+")
plus <- plus %>% rowwise() %>% mutate(start = V5-3000, end =V5+3000)
minus <- subset(gene, gene$V4=="-")
minus <- minus %>% rowwise() %>% mutate(start = V6-3000, end =V6+3000)

tmp <- rbind(plus,minus)
TSS_extended <- tmp[c("V3",  "start",  "end",  "V4",  "V13")]
colnames(TSS_extended) <- c("chr", "start", "end", "strand", "name")
TSS_extended <- unique(makeGRangesFromDataFrame(TSS_extended, keep.extra.columns = T))

#Filtering of PGCLC TSS from all TSS
PGCLC_genes <- read.table(paste0(dir, "scRNAseq_PGCLC_differentiation/PGCLC_genes.txt"), h=F)
colnames(PGCLC_genes) <- c("ID", "gene")
PGCLC_genes <- subset(TSS_extended, name %in% PGCLC_genes$gene)
PGCLC_genes <- makeGRangesFromDataFrame(PGCLC_genes, keep.extra.columns=T)


#Find interactions to PGCLC genes
olap <- findOverlaps(gi, PGCLC_genes)


#All interactions to PGCLC-genes, then exclude interactions with other TSS
first <- cbind(data.frame(anchors(gi[queryHits(olap)])$first), data.frame(PGCLC_genes[subjectHits(olap),]$name))
second <- cbind(data.frame(anchors(gi[queryHits(olap)])$second), data.frame(PGCLC_genes[subjectHits(olap),]$name))
PGCLC_genes_interactions <- rbind(first, second)
PGCLC_genes_interactions <- makeGRangesFromDataFrame(PGCLC_genes_interactions, keep.extra.columns=T)
enhancer_only <- PGCLC_genes_interactions[!PGCLC_genes_interactions %over% TSS_extended]

#Subset the PGCLC-genes-interactions to d6PGCLC H3K27ac peaks, resulting in PGCLC enhancers
d6PGCLC_H3K27ac_peaks <- read.table(paste0(dir, "Enhancer_definition/H3K27ac_peaks_d6PGCLC.bed"))
colnames(d6PGCLC_H3K27ac_peaks) <- c("chr", "start", "end")
PGCLC_enhancer <- subsetByOverlaps(GRanges(d6PGCLC_H3K27ac_peaks), enhancer_only, maxgap = 950)

hits <- findOverlaps(GRanges(d6PGCLC_H3K27ac_peaks), enhancer_only, maxgap = 950)
PGCLC_enhancer_genes <- cbind(DataFrame(GRanges(d6PGCLC_H3K27ac_peaks)[queryHits(hits)]), DataFrame(enhancer_only[subjectHits(hits)]$PGCLC_genes.subjectHits.olap.....name)) 
PGCLC_enhancer_genes <- as.data.frame(PGCLC_enhancer_genes)
colnames(PGCLC_enhancer_genes) <- c("chr", "start", "end", "width", "s", "gene")

write.table(unique(PGCLC_enhancer_genes[,1:3]), row.names = F, col.names = F, quote = FALSE, sep = "\t", file=paste0(dir, "Enhancer_definition/PGCLC_enhancer.bed"))
write.table(unique(PGCLC_enhancer_genes), row.names = F, col.names = T, quote = FALSE, sep = "\t", file=paste0(dir, "Enhancer_definition/PGCLC_enhancer_genes.txt"))
