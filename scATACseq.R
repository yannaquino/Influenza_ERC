################################################################################
################################################################################
# File name: scATACseq.R
# Author:
# Created: mar. 28 avril 2020 11:06:00 CEST
# Last modified: ven. 12 juin 2020 15:28:20 CEST
################################################################################
################################################################################

################################################################################
## Section: preamble
# Set common single-cell R library
.libPaths(
   c('/pasteur/projets/policy01/evo_immuno_pop/R_3.6_scLib',
     '/pasteur/entites/Geh/Shared_Programs/Rlib/3.6.0/')
)

# Load packages
library(dplyr)
library(magrittr)
library(forcats)
library(ggplot2)
library(data.table)
library(SummarizedExperiment)
library(GenomicRanges)
library(diffloop)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(TFBSTools)
library(motifmatchr)
library(circlize)
library(viridis)
library(grid)
library(ComplexHeatmap)
library(Matrix)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(UpSetR)
library(BiocSingular)
library(LiblineaR)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(goseq)
library(GO.db)
# Set data path
project <- "/pasteur/projets/policy01/evo_immuno_pop"
user <- "Mary"
experiment <- "scATACseq_pilot"
sample <- "ATAC2"
# Set 'not in' operator
`%nin%` <- Negate(`%in%`)
# Initialize gg_color_hue function
gg_color_hue <- function(n) {
	  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# Set ggplot theme preferences
theme_set(theme_bw())
theme_update(text=element_text(family="serif",size=12),
	     panel.grid=element_blank()
)
################################################################################

 ################################################################################
 ## Section: quality control and pre-processing
 # For the quality control steps, all data structures will be integrated in a
 # Seurat object, as their pipeline provides comfortable high-level functions.
 
 # Read in the sample's raw peak-barcode cutsite matrix. Each value in the matrix
 # represents the number of Tn5 cut sites for each single barcode that map within
 # each peak
 cutsite <- Seurat::Read10X_h5(
   filename = paste0(project,"/",user,"/",experiment,
 	     "/cellRanger/",sample,"/outs/raw_peak_bc_matrix.h5")
 )
 # Read in meta-data from the singlecell.csv file
 metadata <- read.csv(
   file = paste0(project,"/",user,"/",experiment,
 	 "/cellRanger/",sample,"/outs/singlecell.csv"),
   header = T,
   row.names = 1
 )
 # Build Seurat object
 ATAC <- Seurat::CreateSeuratObject(
   counts = cutsite,
   assay = 'peaks',
   project = sample,
   min.cells = 1,
   meta.data = metadata
 )
 
 # Associate to the Seurat object the corresponding file containing the
 # individual fragments in the library. The advantage of retaining this file is
 # that it contains all fragments associated with each single cell, as opposed to
 # only reads that map to peaks.
 ATAC <- Signac::SetFragments(
   object = ATAC,
   file = paste0(project,"/",user,"/",experiment,
 	 "/cellRanger/",sample,"/outs/fragments.tsv.gz")
 )
 # Associate cellranger's and demuxlet's pre-processing meta-data
 DMX <- readr::read_tsv(
   file = paste0(project,"/",user,"/",experiment,"/demuxlet/popscle/",
          sample,"_demuxlet_8ind.best"),
   col_names = T
 )
 # Filter the ATAC barcodes for those that were detected by demuxlet
 ATAC <- ATAC[,DMX$BARCODE]
 # Add demuxet's droplet type classification as metadata to the barcodes
 ATAC@meta.data$droplet.type <-
   DMX %>% dplyr::filter(BARCODE %in% colnames(ATAC)) %>% dplyr::pull(DROPLET.TYPE)
 # Add demuxlet's most likely identification as metadata to the barcodes
 ATAC@meta.data$individual <-
   DMX %>% dplyr::filter(BARCODE %in% colnames(ATAC)) %>% dplyr::pull(SNG.BEST.GUESS)
 
 # First QC filter: total number of fragments per barcode. Hard filter. Barcodes
 # associated to too few fragments (i.e., <1,000 are uninformative (and
 # classified as AMB by demuxlet). Barcodes associates to too many fragments
 # (i.e., <30,000) are most likely doublets.
 ATAC <- subset(
   ATAC,
   subset = peak_region_fragments >= 1000 & peak_region_fragments < 30000
 )
 
 # All following filters will be soft. Barcodes failing each filter will be
 # flagged; barcodes will be definitively filtered out according to the overlap
 # of failed filters.
 
 # Initialize a logical table. One column per filter, one row per barcode.
 # Filter-wise logical vectors will be added as columns progressively.
 passFilters <- tibble::tibble(
   "Barcodes" = colnames(ATAC)
 ) 
 
 # Second QC filter: fraction of fragments in peak regions. Represents the
 # fraction of total fragments that fall within ATAC-seq peaks. Barcodes with low
 # values (i.e., <15-20%) often represent low-quality cells or technical
 # artifacts that should be removed.
 # Compute fraction of fragments in peak regions per barcode.
 ATAC$pct_reads_in_peaks <-
   ATAC$peak_region_fragments / ATAC$passed_filters * 100
 # Set filter threshold.
 FPRFilter <- 15
 # Plot FPR distribution for SNGs and DBLs
 png(paste0("~/DATA/scATACseqMonocytes/FPR",FPRFilter,sample,".png"),width=11,height=8.5,units="in", res=480)
 ggplot(data = ATAC@meta.data) + 
   geom_point(mapping=aes(x=droplet.type,
                          y=pct_reads_in_peaks,
                          color=log10(peak_region_fragments)),
              position="jitter") +
   geom_violin(mapping=aes(x=droplet.type,
                           y=pct_reads_in_peaks),
               color="black",fill="black",alpha=0.3) +
   geom_hline(yintercept=FPRFilter, linetype=2) +
   xlab("Droplet type") + ylab("Fraction of fragments in peak regions") +
   scale_y_continuous(limits=c(0,100)) +
   scale_color_viridis(name="log10(Fragments in peak regions)",
                       guide=guide_colourbar(title.vjust=0.75)) +
   theme_bw() + theme(text=element_text(family="serif",size=14), legend.position="bottom",
                      panel.grid=element_blank())
 dev.off()
 # Add logical vector of barcodes passing second filter
 passFilters <-
   (ATAC$pct_reads_in_peaks > FPRFilter) %>%
    which(. == T) %>% names() %>%
   (function(x) {
      passFilters %<>% dplyr::mutate("FPR" = case_when(
        Barcodes %in% x ~ TRUE,
        Barcodes %nin% x ~ FALSE
     )
    )
   }
 )
 
 # Third QC filter: ratio of fragments in blacklist regions. The ENCODE project
 # has provided a list of blacklist regions, representing reads which are often
 # associated with artifactual signal.  Cells with a high proportion of reads
 # mapping to these areas (compared to reads mapping to peaks) often represent
 # technical artifacts.
 # Compute blacklist ratio per barcode
 ATAC$blacklist_ratio <-
   ATAC$blacklist_region_fragments/ATAC$peak_region_fragments
 # Set filter threshold
 RBRFilter <- 0.05
 # Add logical vector of barcodes passing third filter
 passFilters <-
   (ATAC$blacklist_ratio < RBRFilter) %>%
    which(. == T) %>% names() %>%
   (function(x) {
      passFilters %<>% dplyr::mutate("RBR" = case_when(
        Barcodes %in% x ~ TRUE,
        Barcodes %nin% x ~ FALSE
     )
    )
   }
  )
 
 # Fourth QC filter:nucleosome banding pattern. The histogram of fragment sizes
 # (determined from the paired-end sequencing reads) should exhibit a strong
 # nucleosome banding pattern. We calculate this per single cell, and quantify
 # the approximate ratio of mononucleosomal (i.e., 200 bp-long) to
 # nucleosome-free (i.e., <100 bp-long) fragments, and store it
 # as nucleosome_signal
 # Compute nucleosome signal
 ATAC <- NucleosomeSignal(object = ATAC)
 # Get nucleosome_signal value marking the 99th percentile. Outlier-defining
 # threshold.
 NSFilter <-
   ATAC$nucleosome_signal %>% quantile(x=.,probs=0.99)
 # Group cells by high or low NS strength. Barcodes which are outliers for the
 # mononucleosomal/nucleosome-free ratio have different banding patterns. The
 # remaining cells exhibit a pattern that is typical for a successful ATAC-seq
 # experiment.
 ATAC$nucleosome_group <- ifelse(
   ATAC$nucleosome_signal > NSFilter,
   "Outliers (99th percentile)",
   "Non outliers"
 )
 # By default, NS is calculated only on chromosome 1
 readsInChr1 <- Signac::GetReadsInRegion(
   object = ATAC,
   region = "chr1-1-2000000",
   group.by = "nucleosome_group"
 )
 # Plot NS distribution for Outliers and Non outliers
 png(paste0("~/DATA/scATACseqMonocytes/NS",sample,".png"),width=11,height=8.5,units="in", res=480)
 ggplot(data = readsInChr1,
        mapping = aes(x = length, fill = group)) +
   geom_histogram(bins = 200) +
   facet_wrap(~group, scales = "free_y") + 
   xlab("Fragment length (bp)") + ylab("Fragments") +
   xlim(0,800) +
   theme_bw() + theme(text=element_text(family="serif",size=14),
                      legend.position="none",
                      panel.grid=element_blank())
 dev.off()
 # Add logical vector of barcodes passing fourth filter
 passFilters <-
   (ATAC$nucleosome_signal < NSFilter) %>%
    which(. == T) %>% names() %>%
   (function(x) {
      passFilters %<>% dplyr::mutate("NS" = case_when(
        Barcodes %in% x ~ TRUE,
        Barcodes %nin% x ~ FALSE
    )
   )
  }
 )
 # Fifth qc filter: tss enrichment score. the enrichment of tn5 integration
 # events at transcription start sites (tsss) can also be an important quality
 # control metric to assess the targeting of tn5 in atac-seq experiments. the
 # tss enrichment score is the number of tn5 integration sites around the tss
 # normalized to the number of tn5 integration
 # sites in flanking regions.
 # Create GRanges object with TSS positions
 gene.ranges <- genes(EnsDb.Hsapiens.v86)
 gene.ranges <- gene.ranges[
   gene.ranges$gene_biotype == 'protein_coding',
   
 ]
 # Get ranges corresponding to TSSs
 tss.ranges <- GRanges(
   seqnames = seqnames(gene.ranges),
   ranges = IRanges(start = start(gene.ranges), width = 2),
   strand = strand(gene.ranges)
 )
 
 seqlevelsStyle(tss.ranges) <- 'UCSC'
 tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')
 # To save time, use only the 2000 first TSSs
 ATAC <- TSSEnrichment(object = ATAC, tss.positions = tss.ranges[1:2000])
 # Get 99th percentile to define outliers
 TSSFilter <-
   ATAC$TSS.enrichment %>% quantile(x=.,probs=0.05)
 ATAC$high.tss <- ifelse(
   ATAC$TSS.enrichment > TSSFilter,
   'High',
   'Low'
 )
 # Barcode-position matrix containing the number of fragments associated to each
 # barcode at each base position in a 2000-bp range centered around the TSS.
 tss.matrix <- ATAC@assays$peaks@misc$TSS.enrichment.matrix
 # Mean accessibility across "High" TSS Enrichment score barcodes at each position
 tss.high <- which(ATAC$high.tss == "High") %>% names(.) %>%
   tss.matrix[.,] %>% Matrix::colMeans()
 # Mean accessibility across "Low" TSS Enrichment score barcodes at each position
 tss.low <- which(ATAC$high.tss == "Low") %>% names(.) %>%
   tss.matrix[.,] %>% Matrix::colMeans()
 # Regroup information in a single table. Restructure table to fit ggplot usage.
 tss.means <-
   tibble(position=-1000:1000,High=tss.high,Low=tss.low) %>%
   tidyr::pivot_longer(data=.,cols=2:3,names_to="Class",values_to="Means")
 # Plot TSS Enrichment score distribution for barcodes with high and low mean TSS
 # Enrichment score
 png(paste0("~/DATA/scATACseqMonocytes/TSSEnrichment",sample,".png"),width=11,height=8.5,units="in", res=480)
 ggplot(data=tss.means) +
   geom_line(mapping=aes(x=position,y=Means,color=Class)) +
   facet_wrap(~Class,
              labeller=labeller(
                Class=setNames(
                       object=c("High TSS Enrichment score",
                                "Low TSS Enrichment score (5th percentile)"),
                       nm=c("High","Low")
                          ))) +
   xlab("Distance from TSS (bp)") +
   ylab("Mean TSS Enrichment score") +
   theme_bw() +
   theme(panel.grid=element_blank(),text=element_text(family="serif",size=14),
         legend.position="none")
 dev.off()
 # Add logical vector of barcodes passing fifth filter
 passFilters <-
   (ATAC$TSS.enrichment > TSSFilter) %>%
    which(. == T) %>% names() %>%
   (function(x) {
      passFilters %<>% dplyr::mutate("TSS" = case_when(
        Barcodes %in% x ~ TRUE,
        Barcodes %nin% x ~ FALSE
    )
   )
  }
 )
 
 # Get overlaps of barcodes not passing QC filters. Store as list (one entry per
 # filter, because that is how UpSetR and VennDiagram like it).
 listFilters <-
 lapply(X=c("FPR","RBR","NS","TSS"), FUN=function(x) {
   passFilters %>% dplyr::filter(!!as.name(x) == F) %>% dplyr::pull(Barcodes)
  }
 )
 names(listFilters) <- c("FPR","RBR","NS","TSS")
 
 # Plot upset barplot with overlaps of barcodes not passing QC filter
 pdf(paste0("~/DATA/scATACseqMonocytes/upSetQC",sample,".pdf"),width=11,height=8.5)
 UpSetR::upset(
   fromList(listFilters),
   order.by = "freq",
   sets.bar.color=gg_color_hue(4),
   sets.x.label="Barcodes filtered out",
   mainbar.y.label="Intersection size"
 )
 dev.off()
 
 # Define barcodes passing QC. Here, I chose to filter out all barcodes that
 # failed either the FPR or the TSS filters
 # Add another logical vector, whether each barcode definitively passed QC
 passFilters %<>% dplyr::mutate("passQC" = case_when(
   FPR==F|TSS==F~F,
   TRUE~T
   )
 )
 
 # Other sources (ArchR Project) seem to agree that only FPR and TSS are
 # pertinent filters.
 # Scatter plot TSS Enrichment score vs. peak region fragments
 
 # Subset ATAC for the barcodes passing QC. Store in new object ATACQC
 ATACQC <-
 passFilters %>%
   dplyr::filter(passQC == T) %>%
   dplyr::pull(Barcodes) %>%
   ATAC[,.]
  
 # Post-QC clean-up
 rm(cutsite,metadata,DMX,readsInChr1,gene.ranges,tss.ranges,tss.matrix,tss.high,
    tss.low,tss.means,listFilters)
 ################################################################################
 
 ################################################################################
 # Section: normalization and dimensionality reduction
 # Term frequency-inverse document frequency (TF-IDF) normalization and latent
 # semantic indexing (LSI) linear dimensionality reduction.  LSI is a technique
 # in natural language processing for analyzing relantionships between a set of
 # documents and the terms they contain by producing a set of concepts related
 # the documents and terms. LSI assumes that words that are close in meaning will
 # occur in similar pieces of text. A matrix containing word counts per document
 # is constructed from a large piece of text and singular value decomposition
 # (SVD) is used to reduce the number of rows while preserving the similarity
 # structure among columns.  In our application, "terms" and "documents" are to
 # be substituted by "peaks/genomic regions" and "cell barcodes" respectively.
 # Hence, LSA on a peak-barcode cutsite count matrix assumes that similar cells
 # will contain chromatin open in peaks that are close in meaning. 
 
 # Extract the peak-barcode cutsite matrix and binarize it
 binaryMatrix <- ATACQC@assays$peaks@counts
 binaryMatrix@x[binaryMatrix@x > 0] <- 1
 
 # TF-IDF normalization. Implementation used by Stuart, Butler et al. 2019, which
 # computes the log of the TF-IDF matrix, adding a pseudocount of 1 to avoid
 # computing the log of 0. 
 # TF-IDF is a two-step normalization procedure that both normalizes across cells
 # to correct for differences in cellular sequencing depth, and across peaks to
 # give higher values to more rare peaks.
 # Compute term frequency. Weigh all sites for individual cells by the total
 # number of sites accessible in that cell.
 TF <- t(t(binaryMatrix) / Matrix::colSums(binaryMatrix))
 # Compute inverse document frequency. Weigh peaks inversely to their occurence
 # in cells.
 IDF <- ncol(binaryMatrix)/Matrix::rowSums(binaryMatrix)
 # Logarithmic scaling of the inverse fraction of barcodes that contain the peak.
 # Applied when performing the standard LSI implementation used by Cusanovich,
 # Hill et al. 2018.
 # IDF <- log(1 + IDF)
 # Compute TF-IDF matrix.
 TFIDF <- Matrix::Diagonal(n=length(IDF),x=IDF) %*% TF
 # Compute the log of TF-IDF matrix, adding a pseudocount of 1 to avoid computing
 # log of 0
 scale.factor <- 10^4
 TFIDF@x <- log1p(x = TFIDF@x * scale.factor)
 
 # Find top features. Calculate the empirical cumulative distribution function of
 # the number of cells each peak (feature) is open in. Keep the top 95%.
 feature.counts <- Matrix::rowSums(TFIDF)
 empirical.distribution <- ecdf(feature.counts)
 feature.distro <- data.frame(
   peaks = names(feature.counts),
   count = feature.counts,
   percentile = empirical.distribution(feature.counts)
 )%>% dplyr::as_tibble()
 # Subset binary matrix to keep only those peaks
 TFIDF <-
 feature.distro %>% 
   dplyr::filter(percentile > 0.05) %>%
   dplyr::pull(peaks) %>%
   TFIDF[.,]
 
 # Find top cells. Same as what was done peak-wise.
 cell.counts <- Matrix::colSums(TFIDF)
 empirical.distribution <- ecdf(cell.counts)
 cell.distro <- data.frame(
   cells = names(cell.counts),
   count = cell.counts,
   percentile = empirical.distribution(cell.counts)
 )%>% dplyr::as_tibble()
 # Subset binary matrix to keep only those cells
 TFIDF <-
 cell.distro %>%
   dplyr::filter(percentile > 0.05) %>%
   dplyr::pull(cells) %>%
   TFIDF[,.]
 
 # Subset ATAQC and the binary matrix to keep the same peaks and cells as in the
 # TFIDF  matrix
 ATACQC <- ATAC[rownames(TFIDF),colnames(TFIDF)]
 binaryMatrix <- binaryMatrix[rownames(TFIDF),colnames(TFIDF)]
 
 # LSI dimensionality reduction. Dimensionality reduction by random IRLBA
 
 # approximation of SVD to 50 dimensions. Because the normalization method is
 # different (i.e., TF-IDF vs. center-scaling) we speak of LSI and not principal
 # components analysis.
 LSI <- BiocSingular::runSVD(
   x=t(TFIDF),k=50,
   BSPARAM = BiocSingular::IrlbaParam(extra.work=50),
   center = F, scale = F)
 
 feature.loadings <- LSI$v
 sdev <- LSI$d / sqrt(x = max(1, nrow(x = LSI$d) - 1))
 cell.embeddings <- LSI$u
 # Scale LSI loadings for each cell to mean 0 and standar deviation 1. These
 # steps are used for learning the weighting of anchors within the scATAC-seq
 # dataset.
 embed.mean <- apply(X = cell.embeddings, MARGIN = 2, FUN = mean)
 embed.sd <- apply(X = cell.embeddings, MARGIN = 2, FUN = sd)
 norm.embeddings <- t((t(cell.embeddings) - embed.mean) / embed.sd)
 # Calculate principal LSI from normalized cell embeddings
 LSI$PCs <- norm.embeddings %*% diag(LSI$d)
 rownames(LSI$PCs) <- colnames(TFIDF)
 
 # Initialize table bearing dimensionality reduction, clustering and other data
 # per cell
 cellData <- tibble::tibble(
   BARCODE = colnames(TFIDF),
   PC1 = LSI$PCs[,1],
   PC2 = LSI$PCs[,2],
   FRAGMENTS = ATACQC$peak_region_fragments,
   DROPLET.TYPE = ATACQC$droplet.type,
   INDIVIDUAL = ATACQC$individual
 )
 
 # Plot cells in LSI space, two first dimensions. Colored by number of fragments
 # in peak regions, faceted by droplet type 
 png(paste0("~/DATA/scATACseqMonocytes/lsi",sample,".png"),width=11,height=8.5,units="in",res=480)
 ggplot(data=cellData) +
   geom_point(mapping=aes(x=PC1,y=PC2,color=FRAGMENTS)) +
   facet_wrap(~DROPLET.TYPE) +
   #scale_x_continuous(trans="log10")
   theme(legend.position="bottom",
         axis.ticks=element_blank(),axis.text=element_blank()) +
   scale_color_continuous(trans="log10") +
   scale_color_viridis(name="log10(Fragments in peak regions)")
 dev.off()
 ################################################################################ 
 
 ################################################################################
 # Cell clustering and feature extraction
 # LSI1 very often captures sequencing depth (technical variation) rather than
 # biological variation. If this is the case, the component should be removed
 # from downstream analyses. To obtain a more biologically relevant custering,
 # for example.
 # Clustering in LSI space using k-nearest neighbors. Find the 20-nearest
 # neighbors for each single cell.
 k.param <- 25
 nn.ranked <- RANN::nn2(
   data=LSI$PCs[,2:16],
   k=k.param,
   searchtype="standard",
   eps=0
 )
 nn.ranked <- nn.ranked$nn.idx
 
 # Convert nn.ranked into a graph
 j <- as.numeric(t(nn.ranked))
 i <- ((1:length(j))-1) %/% k.param + 1
 nn.matrix <- as(
   Matrix::sparseMatrix(
     i=i,j=j,x=1,
     dims=c(nrow(LSI$PCs),nrow(LSI$PCs))
   ),
   Class="Graph")
 rownames(nn.matrix) <- rownames(LSI$PCs)
 colnames(nn.matrix) <- rownames(LSI$PCs)
 neighbor.graphs <- list(nn=nn.matrix)
 
 # Compute shared nearest neighbors
 prune.SNN <- 1/15
 ComputeSNN <- function(nn_ranked, prune) {
     .Call('_Seurat_ComputeSNN', PACKAGE = 'Seurat', nn_ranked, prune)
 }
 snn.matrix <- ComputeSNN(
   nn_ranked = nn.ranked,
   prune = prune.SNN
 )
 rownames(snn.matrix) <- rownames(LSI$PCs)
 colnames(snn.matrix) <- rownames(LSI$PCs)
 snn.matrix <- as.Graph(x = snn.matrix)
 neighbor.graphs[["snn"]] <- snn.matrix
 
 # Run modularity clustering. Function exported from Rcpp package. The 'Rcpp'
 # package provides R functions as well as C++ classes which offer a seamless
 # integration of R and C++.
 # Instantiate function
 RunModularityClusteringCpp <- function(SNN, modularityFunction, resolution,
 				       algorithm, nRandomStarts, nIterations,
 				       randomSeed, printOutput, edgefilename) {
     .Call('_Seurat_RunModularityClusteringCpp',
 	  PACKAGE = 'Seurat',
 	  SNN,
 	  modularityFunction, resolution, algorithm, nRandomStarts,
 	  nIterations, randomSeed, printOutput, edgefilename)
 }
 # Run modularity clustering
 clusters <- RunModularityClusteringCpp(
   SNN=neighbor.graphs$snn,
   modularityFunction=1,
   resolution=0.8,
   algorithm=3,
   nRandomStarts=10,
   nIterations=10,
   randomSeed=0,
   printOutput=T,
   edgefilename=''
 )
 names(clusters) <- colnames(neighbor.graphs$snn)
 clusters %<>% tibble::enframe(x=.,name="BARCODE",value="CLUSTER")
 
 # Add clustering data to cellData
 cellData <-
   merge(x=cellData,y=clusters,by="BARCODE") %>% dplyr::as_tibble()
 rm(clusters)
 
 # Projection from LSI dimensions to UMAP space
 set.seed(140614)
 umapATAC <- umap::umap(
   d=LSI$PCs[,2:16],method="naive",n_components=2,n_neighbors=25,
   metric="euclidean",n_epochs=200,min_dist=0.1,
   spread=1,set_op_ratio_mix_ratio=1,
   local_connectivity=1,negative_sample_rate=5,a=NA,b=NA
 )
 
 # Add UMAP projection data to cellData
 umapCoords <- tibble(
   UMAP1 = umapATAC$layout[,1],
   UMAP2 = umapATAC$layout[,2],
   BARCODE = rownames(umapATAC$layout)
 )
 cellData <-
   merge(x=cellData,y=umapCoords,by="BARCODE") %>% dplyr::as_tibble()
 rm(umapCoords)
 
 # Set color palette for clusters
 getPalette <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))
 clusterPal <-
   setNames(object=getPalette(13),as.character(0:12))
 # Set color palette for individuals
 individualPal <-
   setNames(object=RColorBrewer::brewer.pal(8,"Dark2"),unique(cellData$INDIVIDUAL))
 
 # Plot UMAP projection with points colored according to cluster belonging
 png(paste0("~/DATA/scATACseqMonocytes/umap",sample,".png"),width=11,height=8.5,units="in",res=480)
 ggplot(cellData) +
   geom_point(mapping=aes(x=UMAP1,y=UMAP2,color=as.factor(CLUSTER))) +
   scale_color_manual(name="Cluster",guide=guide_legend(nrow=2,byrow=T),values=clusterPal) +
   theme(legend.position="bottom",axis.ticks=element_blank(),axis.text=element_blank(),panel.spacing=unit(0,"mm"))
 dev.off()
 
 # Plot UMAP projections with points colored according to cluster belonging and
 # faceted by cluster
 png(paste0("~/DATA/scATACseqMonocytes/umap",sample,"Facet.png"),width=11,height=8.5,units="in",res=480)
 ggplot(cellData) +
   geom_point(data=cellData[,-7],mapping=aes(x=UMAP1,y=UMAP2),color="gray") +
   geom_point(mapping=aes(x=UMAP1,y=UMAP2,color=as.factor(CLUSTER))) +
   facet_wrap(~CLUSTER) +
   scale_color_manual(name="Cluster",guide=guide_legend(nrow=2,byrow=T),values=clusterPal) +
   theme(legend.position="bottom",axis.ticks=element_blank(),axis.text=element_blank(),panel.spacing=unit(0,"mm"))
 dev.off()
 
 # Plot UMAP projections with points colored according to droplet type and
 # faceted by cluster
 png(paste0("~/DATA/scATACseqMonocytes/umap",sample,"Droplet.png"),width=11,height=8.5,units="in",res=480)
 ggplot(cellData) +
   geom_point(data=cellData[,-7],mapping=aes(x=UMAP1,y=UMAP2),color="gray") +
   geom_point(mapping=aes(x=UMAP1,y=UMAP2,color=DROPLET.TYPE)) +
   facet_wrap(~CLUSTER) +
   scale_color_discrete(name="Droplet type",aesthetics="color",guide=guide_legend(nrow=2,byrow=T)) +
   theme(legend.position="bottom",axis.ticks=element_blank(),axis.text=element_blank(),panel.spacing=unit(0,"mm"))
 dev.off()
 
 # Plot UMAP projections with points colored according to number of fragments in
 # peak regions and faceted by cluster
 png(paste0("~/DATA/scATACseqMonocytes/umap",sample,"Fragments.png"),width=11,height=8.5,units="in",res=480)
 ggplot(cellData) +
   geom_point(data=cellData[,-7],mapping=aes(x=UMAP1,y=UMAP2),color="gray") +
   geom_point(mapping=aes(x=UMAP1,y=UMAP2,color=FRAGMENTS)) +
   facet_wrap(~CLUSTER) +
   scale_color_continuous(trans="log10") +
   scale_color_viridis(name="log10(Peak region fragments)") +
   theme(legend.position="bottom",axis.ticks=element_blank(),axis.text=element_blank(),panel.spacing=unit(0,"mm"))
 dev.off()
 
 # Plot UMAP projections with points colored according to individual in
 # peak regions and faceted by cluster
 png(paste0("~/DATA/scATACseqMonocytes/umap",sample,"Individual.png"),width=11,height=8.5,units="in",res=480)
 ggplot(cellData) +
   geom_point(data=cellData[,-7],mapping=aes(x=UMAP1,y=UMAP2),color="gray") +
   geom_point(mapping=aes(x=UMAP1,y=UMAP2,color=INDIVIDUAL)) +
   facet_wrap(~CLUSTER) +
   scale_color_manual(name="Cluster",guide=guide_legend(nrow=2,byrow=T),values=individualPal) +
   theme(legend.position="bottom",axis.ticks=element_blank(),axis.text=element_blank(),panel.spacing=unit(0,"mm"))
 dev.off()
 
 # Feature extraction: cluster-defining peaks. Use logistic regression to pull
 # out important peaks (genomic regions) that define clusters. Train a logit
 # model with the 100k accessibility peaks (features) as variables. According to
 # the score at each of the peaks, the model will predict whether each barcode
 # belongs to the cluster it was assigned to by the clustering algorithm. In each
 # prediction, each variable has a predictive weight. If the model is accurate,
 # then the top 500 variables with the highest weight in predicting belonging to
 # each cluster are the 500 genomic regions most informative for that cluster. In
 # other words, the 500 regions defining that cluster.  As normalization does not
 # necessarily improve the accuracy of the model, we use the binary accessibility
 # matrix. It is transposed, to have peaks as columns (variables). As the outcome
 # vector y, we use each barcode's cluster belonging, as defined by the
 # clustering algorithm.
 # Perform 5-fold crossvalidation to test accuracy of the model, as well as best
 # value for cost C.
 cost <- LiblineaR::LiblineaR(
   data=as.matrix(t(binaryMatrix)),target=cellData$CLUSTER,
   type=0,cross=5,
   verbose=T,findC=T
 )
 # Compute the logit model
 logit <- LiblineaR::LiblineaR(
   data=as.matrix(t(binaryMatrix)),target=cellData$CLUSTER,
   type=0,cross=0,
   verbose=T,cost=cost
 )
 
 # Matrix W in logit output bears the weights of all 100k peaks across all
 # clusters. Create a named vector as long as the number of peaks, with the
 # highest predictive value for each peak and the cluster it is associated to.
 weights <-
   sapply(X=1:ncol(logit$W),FUN=function(x){
     logit$W[,x] %>% sort(decreasing=T) %>% head(1)
   }, simplify="vector"
 )
 # Build a table associating each peak to its highest weight and the cluster it
 # is associated to. Group by cluster and only keep the 500 heaviest ones per
 # cluster.
 markerPeaks <- tibble::tibble(
   "PEAK"=colnames(logit$W),
   "CLUSTER"=names(weights),
   "WEIGHT"=weights
 ) %>% dplyr::filter(grepl("^chr",PEAK)) %>%
   dplyr::group_by(CLUSTER) %>%
   dplyr::arrange(CLUSTER,desc(WEIGHT)) %>%
   dplyr::group_modify(~head(.,500))
 
 # Plot accessibility heatmap of 500 cluster-defining peaks across all clustered
 # cells.
 # Compute accessibility z-score. First, revert SVD without PC1, to reconstruct
 # original dimensions, but without the purely technical noise captured by PC1. 
 revSVD <- t(LSI$u[,2:16] %*% diag(LSI$d[2:16]) %*% t(LSI$v[,2:16]))
 # For the heatmap, keep only cluster-defining peaks and clustered barcodes.
 barcodesOrdered <-
   cellData %>% dplyr::arrange(CLUSTER) %>%
   dplyr::pull(BARCODE)
 
 zscoreMatrix <- revSVD[unique(markerPeaks$PEAK),barcodesOrdered]
 # Compute z-score
 zscoreMatrix <- (zscoreMatrix - Matrix::rowMeans(zscoreMatrix))/matrixStats::rowSds(zscoreMatrix)
 
 vecBreaksList <- seq(-1,1,0.001)
 
 funColorScale = colorRamp2(breaks = vecBreaksList,
                      colors = viridis_pal()(length(vecBreaksList)))
 
 # Peaks will be reordered so that they match the cluster order of the barcodes,
 # and they are arranged from highest mean z-score. 
 peaksOrdered <-
   markerPeaks %>% dplyr::ungroup() %>%
     dplyr::mutate(CLUSTER = forcats::fct_relevel(CLUSTER,
						  "0","1","2","3","4","5","6","7","8","9","10","11","12")) %>%
     dplyr::mutate(meanZscore = Matrix::rowMeans(zscoreMatrix)) %>%
     dplyr::arrange(CLUSTER,desc(meanZscore)) %>% dplyr::pull(PEAK)
 
 # Define heatmap annotation
 annoCluster <- cellData$CLUSTER %>% setNames(object=.,nm=cellData$BARCODE) %>%
   .[colnames(zscoreMatrix)]
 col_fun <- rep(clusterPal,table(annoCluster)) %>% setNames(object=.,nm=annoCluster)
 col_ha <- ComplexHeatmap::HeatmapAnnotation(Cluster=annoCluster,col=list(Cluster=col_fun),show_annotation_name=F,show_legend=F)
 col_fun_row<-setNames(object=rep(clusterPal,each=500),nm=as.character(rep(seq(0,12),each=500)))
 rowAnno <- setNames(as.character(rep(seq(0,12),each=500)),as.character(rep(seq(0,12),each=500))) 
 row_ha <- ComplexHeatmap::HeatmapAnnotation(Cluster=rowAnno,col=list(Cluster=col_fun_row),show_annotation_name=F,show_legend=F,which="row")
 
 # Plot heatmap
 png(paste0("~/DATA/scATACseqMonocytes/heatmap",sample,".png"),width=11,height=8.5,units="in",res=360)
 ComplexHeatmap::Heatmap(
   matrix=zscoreMatrix[peaksOrdered,],
   name="Z score",
   top_annotation=col_ha,
   left_annotation=row_ha,
   col=funColorScale,
   cluster_rows=F,
   cluster_columns=F,
   show_row_dend=F,
   show_column_dend=F,
   show_row_names=F,
   show_column_names=F
 )
 dev.off()
 
  # Feature extraction: individual-defining peaks.
  # Compute cost function.
  cost <- LiblineaR::LiblineaR(
    data=as.matrix(t(binaryMatrix)[cellData$BARCODE,]),target=cellData$INDIVIDUAL,
    type=0,cross=5,
    verbose=T,findC=T
  )
  
  # Build model.
  logitInd <- LiblineaR::LiblineaR(
    data=as.matrix(t(binaryMatrix)[cellData$BARCODE,]),target=cellData$INDIVIDUAL,
    type=0,cross=0,
    verbose=T,cost=cost
  )
  
  # Prioritize peaks by predictive weight
  weightsInd <-
    sapply(X=1:ncol(logitInd$W),FUN=function(x){
      logitInd$W[,x] %>% sort(decreasing=T) %>% head(1)
    }, simplify="vector"
  )
  
  # Extract cluster-defining peaks
  markerPeaksInd <- tibble::tibble(
    "PEAK"=colnames(logitInd$W),
    "CLUSTER"=names(weightsInd),
    "WEIGHT"=weightsInd
  ) %>% dplyr::filter(grepl("^chr",PEAK)) %>%
    dplyr::group_by(CLUSTER) %>%
    dplyr::arrange(CLUSTER,desc(WEIGHT)) %>%
    dplyr::group_modify(~head(.,500))
  
  # Return to original dimensions, without PC1  
  revSVD <- t(LSI$u[,-1] %*% diag(LSI$d[-1]) %*% t(LSI$v[,-1]))
  
  barcodesOrdered <-
    cellData %>% dplyr::arrange(INDIVIDUAL) %>%
    dplyr::pull(BARCODE)
  
  zscoreMatrix <- revSVD[unique(markerPeaks$PEAK),barcodesOrdered]
  
  zscoreMatrix <- (zscoreMatrix - Matrix::rowMeans(zscoreMatrix))/matrixStats::rowSds(zscoreMatrix)
  
  vecBreaksList <- seq(-1,1,0.001)
  
  funColorScale = colorRamp2(breaks = vecBreaksList,
                       colors = viridis_pal()(length(vecBreaksList)))
  
  peaksOrdered <-
    markerPeaks %>% dplyr::ungroup() %>%
      dplyr::mutate(CLUSTER = forcats::fct_relevel(CLUSTER,"AFB022","AFB108","AFB131","AFB176","EUB058","EUB060","EUB061","EUB078")) %>%
      dplyr::mutate(meanZscore = Matrix::rowMeans(zscoreMatrix)) %>%
      dplyr::arrange(CLUSTER,desc(meanZscore)) %>% dplyr::pull(PEAK)
  
  annoCluster <- cellData$INDIVIDUAL %>% setNames(object=.,nm=cellData$BARCODE) %>%
    .[colnames(zscoreMatrix)]
  col_fun <- rep(individualPal,table(annoCluster)) %>% setNames(object=.,nm=annoCluster)
  col_ha <- ComplexHeatmap::HeatmapAnnotation(Cluster=annoCluster,col=list(Cluster=col_fun),show_annotation_name=F,show_legend=F)
  col_fun_row<-setNames(object=rep(individualPal,each=500),nm=as.character(rep(c("AFB022","AFB108","AFB131","AFB176","EUB058","EUB060","EUB061","EUB078"),each=500)))
  rowAnno <-
  	setNames(as.character(rep(c("AFB022","AFB108","AFB131","AFB176","EUB058","EUB060","EUB061","EUB078"),each=500)),as.character(rep(c("AFB022","AFB108","AFB131","AFB176","EUB058","EUB060","EUB061","EUB078"),each=500))) 
  row_ha <- ComplexHeatmap::HeatmapAnnotation(Cluster=rowAnno,col=list(Cluster=col_fun_row),show_annotation_name=F,show_legend=T,which="row")
  
  # Plot heatmap
  png(paste0("~/DATA/scATACseqMonocytes/heatmapIndividuals",sample,".png"),width=11,height=8.5,units="in",res=360)
    ComplexHeatmap::Heatmap(          
    matrix=zscoreMatrix[peaksOrdered,],
    name="Z score",
    top_annotation=col_ha,
    left_annotation=row_ha,
    col=funColorScale,
    cluster_rows=F,
    cluster_columns=F,
    show_row_dend=F,
    show_column_dend=F,
    show_row_names=F,
    show_column_names=F
  )
  dev.off()
###############################################################################

###############################################################################
 # Gene ontology enrichment
 # Annotate marker peaks.
 clusterAnnotation <-
 markerPeaks %>%
   dplyr::pull(PEAK) %>%
   base::strsplit(x=., split=":|-") %>% unlist() %>%
   matrix(data=.,ncol=3, byrow=T) %>% as.data.frame(stringsAsFactors=F) %>%
   magrittr::set_colnames(c("seqname","start","end")) %>%
   makeGRangesFromDataFrame() %>%
   annotatePeak(peak=.,
   TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
   level="gene",
   annoDb="org.Hs.eg.db")
 
 # Add annotation data to peaks.
 markerPeaks$Symbol <- clusterAnnotation@anno$SYMBOL
 markerPeaks$entrezId <- clusterAnnotation@anno$geneId
 markerPeaks$ensGene <- clusterAnnotation@anno$ENSEMBL
 markerPeaks$distanceToTSS <- clusterAnnotation@anno$distanceToTSS
 
 # Define background peaks. All accessibility peaks. Annotate the peaks.
 backgroundPeaks <-
   rownames(ATACQC) %>%
   base::strsplit(x=., split=":|-") %>% unlist() %>%
   matrix(data=.,ncol=3, byrow=T) %>% as.data.frame(stringsAsFactors=F) %>%
   magrittr::set_colnames(c("seqname","start","end")) %>% makeGRangesFromDataFrame()
 
 backgroundAnnotation <-
   annotatePeak(peak=backgroundPeaks,
   TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
   level="gene",
   annoDb="org.Hs.eg.db")
 
 backgroundPeaks$Symbol <- backgroundAnnotation@anno$SYMBOL
 backgroundPeaks$entrezId <- backgroundAnnotation@anno$geneId
 backgroundPeaks$ensGene <- backgroundAnnotation@anno$ENSEMBL
 backgroundPeaks$distanceToTSS <- backgroundAnnotation@anno$distanceToTSS
 
 # Define set of background genes. Genes associated to peaks for which the
 # nearest TSS is 1,000 bp away or less.
 backgroundGenes <-
   backgroundPeaks[
     elementMetadata(backgroundPeaks)[,"distanceToTSS"] %in% seq(-1000,1000)
     ]$Symbol %>% unique() %>%
   (function(x) x[which(!is.na(x))])
 
 # For loop. Perform GO enrichment analysis for each cluster. Maxime's function,
 # will comment later.
G2S=function(genes){
  GeneAnnot[match(genes,GeneAnnot[,1]),'Associated.Gene.Name']
}
odds.ratio=function(tab, alpha = 0.05){
  test=fisher.test(tab,conf.level=1-alpha)
  oframe <- data.frame(LowerCI = test$conf.int[1],
		       OR = test$est,
		       UpperCI = test$conf.int[2],
		       alpha = alpha,P=test$p.value)
  oframe
			      }
  counts2tab=function(xy,x,y,tot){matrix(c(tot-x-y+xy,x-xy,y-xy,xy),2)}

   
 for (cluster in sort(unique(markerPeaksTSS$CLUSTER))) {
   
   geneList <-
     markerPeaksTSS %>%
    # dplyr::filter(distanceToTSS %in% seq(-1000,1000)) %>%
     dplyr::filter(CLUSTER==cluster) %>%
     dplyr::pull(name) %>% unique()
 
   DE <- backgroundGenes%in%geneList
   names(DE) <- backgroundGenes
   nullP <- nullp(DE,'hg38','geneSymbol',plot.fit=F)
   res <- goseq(nullP,"hg38","geneSymbol")
   resUp=res
   resUp$FDR=p.adjust(res$over_represented_pvalue,'fdr')
   FDR <- 0.05
   resUp=resUp[resUp$FDR<FDR,]
   resUp$Term=sapply(mget(resUp[,1],GOTERM),function(x){x@Term})
   resUp$FoldEnrich=resUp[,4]/resUp[,5]/mean(DE)
   res <- resUp
   if(nrow(res)>0){
       	res$nbInGrp=rep(length(geneList),nrow(res))
       	res$nbInBckgd=rep(length(backgroundGenes),nrow(res))
       	OR=sapply(as.data.frame(t(mapply(function(xy,x,y,tot){odds.ratio(counts2tab(xy,x,y,tot))},res$numDEInCat,res$numInCat,length(geneList),res$nbInBckgd))),as.numeric)
       	if(is.matrix(OR)){
       		OR=as.data.frame(OR)
       	}else{
   			nn=names(OR)
       		OR=as.data.frame(matrix(OR,1))
       		colnames(OR)=nn
       	}
       	
       	OR$CI=paste('[',round(OR[,'LowerCI'],1),'-',round(OR[,'UpperCI'],1),']',sep='')
   		res$OR=OR$OR
   		res$CI=OR$CI
   		res$lowerCI=OR$LowerCI
       		}		 
   if (nrow(res)>0) {
   partial <-
   res %>% dplyr::as_tibble() %>%
     dplyr::filter(ontology == "BP") %>%
     dplyr::arrange(FDR,desc(lowerCI)) %>%
     dplyr::select(Term,FDR,FoldEnrich,lowerCI,category) %>%
     dplyr::mutate(CLUSTER = cluster)
   } else {
   partial <- NULL
   }
   assign(x=paste0("partial",cluster),value=partial)
 }
 
 listGO <-
   sapply(X=sort(unique(cellData$CLUSTER)),
 	 FUN=function(x) {get(paste0("partial",x))})
 names(listGO) <- paste0("Cluster ", sort(unique(cellData$CLUSTER)))
 
 # Plot top 20 GO enrichments for clusters with significantly enriched terms
 for (cluster in sort(unique(cellData$CLUSTER))) {
   if (!is.null(listGO[[as.character(cluster)]])) { 
     pdf(paste0("~/DATA/scATACseqMonocytes/GO",sample,"Cluster",cluster,".pdf"),width=11,height=8.5)
     listGO[[as.character(cluster)]] %>% dplyr::filter(CLUSTER == cluster) %>% dplyr::arrange(desc(lowerCI)) %>% dplyr::slice(1:20) %>%
     ggplot(data=.) +
     geom_col(mapping=aes(x=FoldEnrich,y=reorder(Term,sort(FoldEnrich,decreasing=T)),fill=lowerCI)) +
     scale_color_viridis(aesthetics="fill",name="Lower CI FE") +
     ylab("") +
     xlab("Fold enrichment") +
    theme(legend.position="right")
     dev.off()
   }
 }
################################################################################

################################################################################
# TFBM differential accessibility
# Create Seurat assay object with TFIDF matrix to use Signac's chromVAR
# adaptation.
ATACQC <- readRDS(file="~/DATA/scATACseqMonocytes/ATACQC.RDS")
TFIDF <- readRDS(file="~/DATA/scATACseqMonocytes/TFIDF_ATAC2.RDS")

 ATACQC[['TFIDF']] <- Seurat::CreateAssayObject(counts=TFIDF)
 Seurat::DefaultAssay(ATACQC) <- 'TFIDF'
 
 # Get a list of motif position frequency matrices from the JASPAR database
 pfm <- getMatrixSet(
   x = JASPAR2018,
   opts = list(species = 9606,
 	      all_versions = FALSE)
 )
 
 # Scan the DNA sequence of each peak for the presence of each motif
 motif.matrix <- Signac::CreateMotifMatrix(
   features = StringToGRanges(rownames(ATACQC),
 			     sep = c(":", "-")),
   pwm = pfm,
   genome = 'hg38',
   sep = c(":", "-"),
   use.counts = FALSE
 )
 
 # Create a new Mofif object to store the results
 motif <- Signac::CreateMotifObject(
   data = motif.matrix,
   pwm = pfm
 )
 
 # Add the Motif object to the assay
 ATACQC[['peaks']] <- Signac::AddMotifObject(
   object = ATACQC[['peaks']],
   motif.object = motif
 )
 
 # In order to test for overrepresented motifs, we also need to compute some
 # sequence characteristics of the peaks, such as GC content, sequence length,
 # and dinucleotide frequency. The RegionStats function computes this information
 # for us and stores the results in the feature metadata in the Seurat object.
 ATAQC <- Signac::RegionStats(
   object = ATACQC,
   genome = BSgenome.Hsapiens.UCSC.hg38,
   sep = c(":", "-")
 )
 
 # Run Signac's chromVAR adaptation
 Seurat::DefaultAssay(ATACQC) <- 'peaks'
 
 ATACQC <- Signac::RunChromVAR(
   object=ATACQC,
   genome=BSgenome.Hsapiens.UCSC.hg38
 )
 
 saveRDS(object=ATACQC,file=paste0("~/DATA/scATACseqMonocytes/ATACQCChrom_ATAC2",sample,".RDS"))
################################################################################

################################################################################
# Cicero
# cellinfo <-
  cellDataATAC2 %>% dplyr::select(BARCODE) %>% data.frame()
rownames(cellinfo) <- cellinfo$BARCODE
names(cellinfo) <- "cells"
cellinfo <- as(cellinfo,"AnnotatedDataFrame")

peakinfo <- rownames(binaryMatrixATAC2) %>%
stringr::str_split(pattern=":|-",simplify=T) %>% as.data.frame()
colnames(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name
peakinfo <- as(peakinfo,"AnnotatedDataFrame")

binaryMatrixATAC2 <- binaryMatrixATAC2[,rownames(cellinfo)]
row.names(binaryMatrixATAC2) <- row.names(peakinfo)

input_cds <-  suppressWarnings(
  newCellDataSet(cellData = binaryMatrixATAC2,
		 phenoData = cellinfo,
		 featureData = peakinfo))
input_cds <- monocle::detectGenes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
input_cds <- estimateSizeFactors(input_cds)

umap_coords <- cellDataATAC2 %>% dplyr::select(UMAP1,UMAP2) 
colnames(umap_coords) <- c("umap_ccord1","umap_coord2")
rownames(umap_coords) <- cellDataATAC2$BARCODE

cicero_cds <- make_cicero_cds(input_cds,reduced_coordinates=umap_coords)

load("~/hg38ChromLengths.RDS")

genomic_coordsHg38 <-
	data.frame(seqnames=paste0("chr",c(1:22,"X","Y")),lengths=lengths)

conns <- run_cicero(cds=cicero_cds,genomic_coords=genomic_coordsHg38)

ccansATAC2 <- generate_ccans(conns)

connsCcanATAC2 <-
	merge(x=conns,y=cansATAC2,by.x="Peak1",by.y="Peak",all.x=T,all.y=T)

connsCcanATAC2 <- connsCcanATAC2[-which(is.na(connsCcanATAC2$coaccess)),]

connsCcanATAC2$Peak2Vector <- as.vector(connsCcanATAC2$Peak2)

connsCcanATAC2Mod <-
connsCcanATAC2 %>% dplyr::mutate(Peak1Ordered = case_when(
  as.numeric(str_replace_all(pattern=c("chr[[1-9]+XY]_"="","_[0-9]+$"=""),string=Peak1))<as.numeric(str_replace_all(pattern=c("chr[[1-9]+XY]_"="","_[0-9]+$"=""),string=Peak2Vector))~Peak1,
  as.numeric(str_replace_all(pattern=c("chr[[1-9]+XY]_"="","_[0-9]+$"=""),string=Peak1))>as.numeric(str_replace_all(pattern=c("chr[[1-9]+XY]_"="","_[0-9]+$"=""),string=Peak2Vector))~Peak2Vector
)) %>%
dplyr::mutate(Peak2Ordered = case_when(
  as.numeric(str_replace_all(pattern=c("chr[[1-9]+XY]_"="","_[0-9]+$"=""),string=Peak1))<as.numeric(str_replace_all(pattern=c("chr[[1-9]+XY]_"="","_[0-9]+$"=""),string=Peak2Vector))~Peak2Vector,
  as.numeric(str_replace_all(pattern=c("chr[[1-9]+XY]_"="","_[0-9]+$"=""),string=Peak1))>as.numeric(str_replace_all(pattern=c("chr[[1-9]+XY]_"="","_[0-9]+$"=""),string=Peak2Vector))~Peak1
))

Chrom <- connsCcanATAC2Mod$Peak1Ordered %>% stringr::str_split(pattern="_",simplify=T) %>% .[,1]
Chromstart <- connsCcanATAC2Mod$Peak1Ordered %>% stringr::str_split(pattern="_",simplify=T) %>% .[,2]
Chromend <- connsCcanATAC2Mod$Peak2Ordered %>% stringr::str_split(pattern="_",simplify=T) %>% .[,3]
name <- paste0("CCAN:",connsCcanATAC2Mod$CCAN)
score <- round(range01(connsCcanATAC2Mod$coaccess)*1000,0)
value <- connsCcanATAC2Mod$coaccess
exp <- rep(".",nrow(connsCcanATAC2Mod))
color <- rep(0,nrow(connsCcanATAC2Mod))
sourceChrom <- connsCcanATAC2Mod$Peak1Ordered %>% stringr::str_split(pattern="_",simplify=T) %>% .[,1]
sourceStart <- connsCcanATAC2Mod$Peak1Ordered %>% stringr::str_split(pattern="_",simplify=T) %>% .[,2]
sourceEnd <- connsCcanATAC2Mod$Peak1Ordered %>% stringr::str_split(pattern="_",simplify=T) %>% .[,3]
sourceName <- rep(".",nrow(connsCcanATAC2Mod))
sourceStrand <- rep(".",nrow(connsCcanATAC2Mod))
targetChrom <- connsCcanATAC2Mod$Peak2Ordered %>% stringr::str_split(pattern="_",simplify=T) %>% .[,1]
targetStart <- connsCcanATAC2Mod$Peak2Ordered %>% stringr::str_split(pattern="_",simplify=T) %>% .[,2]
targetEnd <- connsCcanATAC2Mod$Peak2Ordered %>% stringr::str_split(pattern="_",simplify=T) %>% .[,3]
targetName <- rep(".",nrow(connsCcanATAC2Mod))
targetStrand <- rep(".",nrow(connsCcanATAC2Mod))

bigInteract <- data.frame(
  Chrom,
  Chromstart,
  Chromend,
  name,
  score,
  value,
  exp,
  color,
  sourceChrom,
  sourceStart,
  sourceEnd,
  sourceName,
  sourceStrand,
  targetChrom,
  targetStart,
  targetEnd,
  targetName,
  targetStrand
)

readr::write_tsv(x=bigInteract,path="~/DATA/scATACseqMonocytes/connsATAC2.inter.bb",col_names=F)
