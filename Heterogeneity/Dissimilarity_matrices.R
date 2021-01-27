library("PDclust")

#set for region to be processed (e.g. PGCLC for PGCLC_enhancer)
args = commandArgs()
set <- args[7]
print(set)

#get all single-cell DNA methylation bedGraph files
list <- list.files(pattern = paste("*_met_", set, "_enhancer.bedGraph", sep=""))
data <- lapply(setNames(list, make.names(gsub(paste("*_met_", set, "_enhancer.bedGraph$", sep=""), "", list))), read.table)
colnames_data <- c("chr","start","end", "meth") 
data <- lapply(data, setNames, colnames_data)

#Calculate the mean mCpG levels for all region in a single-cell
mCpG <- lapply(data, function(x) {mean(x$meth, na.rm = TRUE)})
mCpG <- as.data.frame(unlist(mCpG))
rownames(mCpG) <- substr(rownames(mCpG), 0, 10)
colnames(mCpG) <- "mCpG_local"

#Get the CpG coverage
cov <- lapply(data, nrow)
cov <- as.data.frame(unlist(cov))
rownames(cov) <- substr(rownames(cov), 0, 10)
colnames(cov) <- "coverage_local"

#Create meta data data frame of mean CpG methylation and coverage
meta_data <- merge(mCpG, cov, by=0, all=T)
rownames(meta_data) <- meta_data$Row.names
meta_data$Row.names <- NULL
save(meta_data, file=paste("meta_data_local_scM_Epiblast_", set, ".Rdata", sep=""))

#Check the pairwise comparisions of CpGs
data_pairwise <- create_pairwise_master(data)
pdf(paste("Number of CpGs in two samples of ", set, "for diss.pdf", sep=""), height = 4, width = 4)
hist(log10(data_pairwise$num_cpgs), xlab = "log10 (number of bi-covered CpGs)", main = "", col = "blue")
dev.off()

#Order and save the dissimilarity matrix
diss <- convert_to_dissimilarity_matrix(data_pairwise)
rownames(diss) <- substr(rownames(diss), 0, 10)
colnames(diss) <- substr(colnames(diss), 0, 10)
diss_ordered <- diss[order(row.names(diss)),order(colnames(diss))]
save(diss_ordered, file=paste("diss_scM_Epiblast_", set, ".Rdata", sep=""))
