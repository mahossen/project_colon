#Load necessary libraries
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(AnnotationHub)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(dplyr)
seed(12345)

#Define variables
analysis_parent_folder <- "./scRNA_processing/"
setwd(analysis_parent_folder)
clinicaldata<-"./metadata.csv"

#creating a Seurat object
counts <- Read10X_h5(filename = "./data/filtered_peak_bc_matrix.h5")

metadata <- read.csv(
  file = "./data/singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "./data/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

atac_data <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  project = 'ATAC',
  meta.data = metadata
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(atac_data) <- annotations

#############################################################
#Quality control & Filtering.................................
#############################################################

# compute nucleosome signal score per cell
atac_data <- NucleosomeSignal(object = atac_data)

# compute TSS enrichment score per cell
atac_data <- TSSEnrichment(object = atac_data)

# add fraction of reads in peaks
atac_data$pct_reads_in_peaks <- atac_data$peak_region_fragments / atac_data$passed_filters * 100

# add blacklist ratio
atac_data$blacklist_ratio <- FractionCountsInRegion(
					  object = atac_data, 
					  assay = 'peaks',
					  regions = blacklist_hg38_unified
					)

#Scatterplot
DensityScatter(atac_data, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

#Histogram...........
atac_data$nucleosome_group <- ifelse(atac_data$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = atac_data, group.by = 'nucleosome_group')

#VinPlot...........
VlnPlot(
  object = atac_data,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 5
)

#filter the aggregated scATACseq object using empirically-determined QC parameters
atac_data <- subset(
		  x = atac_data,
		  subset = nCount_peaks > 2000 &
			nCount_peaks < 100000 &
			pct_reads_in_peaks > 40 &
			blacklist_ratio < 0.01 &
			nucleosome_signal < 4 &
			TSS.enrichment > 4
		)
atac_data

#########################################################################################################
#Normalization and linear dimensional reduction
#########################################################################################################
#linear dimensional reduction
atac_data <- RunTFIDF(atac_data)
atac_data <- FindTopFeatures(atac_data, min.cutoff = 'q75')
atac_data <- RunSVD(atac_data)

#non-linear dimensional reduction...........
# ElbowPlot(atacAggr, ndim = 40) # select number of dimensions for UMAP embedding
atac_data <- RunUMAP(object = atac_data, reduction = 'lsi', dims = 2:30)
atac_data <- FindNeighbors(object = atac_data, reduction = 'lsi', dims = 2:30)
atac_data <- FindClusters(object = atac_data, verbose = FALSE, algorithm = 3)
DimPlot(object = atac_data, label = TRUE) + NoLegend()

##########################################################################################
#Create a gene activity matrix
#########################################################################################

#Create a gene activity matrix
gene.activities <- GeneActivity(atac_data)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
atac_data[['RNA']] <- CreateAssayObject(counts = gene.activities)
atac_data <- NormalizeData(
		  object = atac_data,
		  assay = 'RNA',
		  normalization.method = 'LogNormalize',
		  scale.factor = median(atac_data$nCount_RNA)
		)

saveRDS(atac_data, file ="atac_data.rds")

