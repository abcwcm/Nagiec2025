Transcriptomics
=======================

Processing and differential expression determination
-----------------------------------------------------

Transcriptomics data was generated with GeneChip Human Transcriptome Array 2.0 (Affymetrix), which contains 70% exon probes and 30% exon-exon spanning probes.
The CEL files were read into [R](https://www.R-project.org/) with the help of the [`oligo` package](https://pubmed.ncbi.nlm.nih.gov/20688976/), the Platform Design Info for Affymetrix HTA-2_0 package ([pd.hta.2.0](http://www.bioconductor.org/packages/release/data/annotation/html/pd.hta.2.0.html)), and the annotation package [`hta20transcriptcluster.db`]((http://www.bioconductor.org/packages/release/data/annotation/html/hta20transcriptcluster.db.html)).
The raw intensity values were normalized with the [Robust Multichip Average algorithm](dx.doi.org/10.1093/biostatistics/4.2.249) as implemented in the `oligo::rma()` function.

Statistically significantly differentially expressed genes were determined with [`limma`](dx.doi.org/10.1093/nar/gkv007),
comparing the values of the DN replicates to those of the control conditions (EV, YA) for every day of the time course, using an adjusted p-value threshold of 5%.
The resulting log-transformed fold change (logFC) values were used for downstream analyses and visualizations.
In addition, we applied the Transcript Time Course Analysis algorithm as implemented in the [`TTCA package`](https://pubmed.ncbi.nlm.nih.gov/28088176/) to identify genes that change when taking the entire time course into consideration.

Clustering and dimensionality reduction
----------------------------------------

We wanted to identify groups of genes that showed similar expression changes across time as well as functional similarities.
To this end, we used the function cluster_walktrap (k = 60) of the [`igraph` package](https://static1.squarespace.com/static/5b68a4e4a2772c2a206180a1/t/5cd1e3cbb208fc26c99de080/1557259212150/c1602a3c126ba822d0bc4293371c.pdf) to identify clusters of genes with similar changes across time focussing on genes that were found to be differentially expressed by both limma and TTCA analyses (see above).

To assess which number of neighbors (*k*) returned robust results, we additionally obtained tSNE coordinates with the [`Rtsne` package](https://cran.r-project.org/web/packages/Rtsne/index.html) using the logFC values for each day of the time course.
We subsequently collated the 12 clusters identified this way into 7 clusters based on (i) the direction of the expression changes and (ii) great similarity of GO term enrichments (see below). 

Gene set overrepresentation analyses
------------------------------------

To determine overrepresented gene sets as shown in the manuscripts, we used the R package [`clusterProfiler`](https://pubmed.ncbi.nlm.nih.gov/22455463/) with the functions `compareCluster()` and `enricher()` and the gene sets as defined by the [MSigDB Hallmark gene set](https://pubmed.ncbi.nlm.nih.gov/22455463/) (MSigDB database v6).
For visualization, the clusterProfiler functions `cnetplot()`, `heatplot()`, and `dotplot()` were used. 

Proteomics
==============

Identifying significantly changing proteinsThe protein abundances were loadedinto R and re-normalized, using a method established in the Krumsiek Lab. In brief, for every protein, its median abundances across all samples is determined, followed by a dilution factor, which is based on the median-normalized values. For details, see SuppFigS3_px_plots_proteome.Rmd in the github repo.To identify proteins with statistically significantlychanging values across the different conditions, the limma package was used [Ref]. Missing values were replaced with a pseudo-value (log2(2^min(<observed values>))/2) before limma::lmFit was used to fit a linear model to the observed protein abundances over all replicates and samples of interest. LogFC values were obtained via estimating the coefficients of the model and statistical significantly changing proteins were identified via limma::eBayes() function and limma::topTable() with Benjamini-Hochberg correction for multiple testing and a defined p value cut-off of 0.1 for results presented in Supplementary Figure S3 and 0.05 for the time course analyses.ClusteringTo identify groups of proteins with similar changes over time, we applied kNN-based clustering to the logFC values for DN-vs-YA as calculated by limma (described in previous paragraph): igraph::cluster_walktrap(scran::buildSNNGraph(t(prots.lf), k =40))

Metabolites
==============

Metabolite raw values were median-normalized (per run) and quotient-normalized (across all samples; according to [Dieterle et al., 2006]( https://www.ncbi.nlm.nih.gov/pubmed/16808434), missing values were imputed using a KNN-based approach (following the approach described as the most robust [here](https://pubmed.ncbi.nlm.nih.gov/30830398/).
Limma was applied to the log2-transformed, normalized values to determine significantly changing metabolites in pairwise comparisons for each day of the time course. 

