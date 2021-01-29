Content
-------

* `sample_metadata.txt`: Meta data file for the [scNMT data](https://github.com/rargelaguet/scnmt_gastrulation)
* `GEO_id_met.txt`: Linkage of the GEO id with the meta data file
* `meta_data_epiblast_scMT.Rdata`: Processed meta data file containing the GEO entry, average CpG methylation levels and stages for all epiblast cells (E4.5-E6.5)
* `DNAmethylation_genome.Rdata`: Genome-wide average of CpG methylation and its coverage in each single cell.
* `DNAmethylation_PGCLC_enhancer.Rdata`: Average of CpG methylation and its coverage in each single cell for all PGCLC enhancers
* `Dissimilarity_matrices.R`: Calculation of the dissimilarity matrices with [PDclust](https://github.com/hui-tony-zk/PDclust)
* `diss_scM_Epiblast_ENHANCER.Rdata` : Output of the `Dissimilarity_matrices.R` containing the dissimilarity matrix for the for all CpGs of the ENHANCERS
* `DNAmethylation_Heterogeneity_scMT.R`: R-script for the mCpG methylation heterogeneity
* `GSE121650_rna_counts.tsv`: single-cell RNA-seq count matrix of the scNMT data
* `scMethylation_Transcription.R`: single-cell transcriptomes and DNA methylation analysis
