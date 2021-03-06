# Enhancer-associated H3K4 methylation safeguards in vitro germline competence
In the following folders scripts are depsoited that were used for the investigation of germline competence: [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.07.07.192427v2)

# Summary
Germline specification in mammals occurs through an inductive process whereby competent cells in the post-implantation epiblast differentiate into primordial germ cells (PGC). The intrinsic factors that endow epiblast cells with the competence to respond to germline inductive signals remain unknown. Single-cell RNA sequencing across multiple stages of an *in vitro* PGC-like cells (PGCLC) differentiation system shows that PGCLC genes initially expressed in the naïve pluripotent stage become homogeneously dismantled in germline competent epiblast like-cells (EpiLC). In contrast, the decommissioning of enhancers associated with these germline genes is incomplete. Namely, a subset of these enhancers partly retain H3K4me1, accumulate less heterochromatic marks and remain accessible and responsive to transcriptional activators. Subsequently, as *in vitro* germline competence is lost, these enhancers get further decommissioned and lose their responsiveness to transcriptional activators. Importantly, using H3K4me1 deficient cells, we show that the loss of this histone modification reduces the germline competence of EpiLC and decreases PGCLC differentiation efficiency. Our work suggests that, although H3K4me1 might not be essential for enhancer function, it can facilitate the (re)activation of enhancers and the establishment of gene expression programs during specific developmental transitions.

<p align="center">
<img src="image/Model_Germline_competence.png" width="750" align="center">
</p>

*Figure: Model illustrating the partial decommissioning of PGCLC enhancers during PGCLC differentiation.*


## Content
* `/scRNAseq_PGCLC_differentiation/`: Analysis of the single-cell RNA-seq data of the PGCLC differentiation
* `/Enhancer_definition/`: Assignment of EpiLC, EpiSC and PGCLC enhancers
* `/Epigenetic_comparision/`: Comparative ChIP-seq analysis for germline competence
* `/scRNAseq_dCD/`: Analysis of the single-cell RNA-seq data of the H3K4me1-deficient dCD d4 EB
* `/Heterogeneity/`: CpG methylation heterogeneity anaylsis for PGCLC enhancer using [scNMT data](https://doi.org/10.1038/s41586-019-1825-8)

## Data
All data from the different stages of PGCLC differentiation have been deposited in the [GEO database](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155089) under the follwoing accession numbers:
* [GSE155058 - ATAC-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155058)
* [GSE155062 - ChIP-seq (Comparision EpiLC vs. EpiSC)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155062)
* [GSE155069 - ChIP-seq (*Mll3/4 dCD* and *Otx2* KO)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155069)
* [GSE155079 - bulk RNA-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155079)
* [GSE155083 - Genome-wide DNA methylation](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155083)
* [GSE155088 - scRNA-seq (10x Genomics)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155088)

## Contact
Álvaro Rada-Iglesias (alvaro.rada@unican.es)  
Tore Bleckwehl (tbleckwe@uni-koeln.de)
