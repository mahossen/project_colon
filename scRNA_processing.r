# Import libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(DoubletFinder)
set.seed(1)


# Define variables
sample_name <- "all_samples"
individual_qc_and_dublet_plot_location <- "./scRNA_processing/doublet_analysis/"
analysis_parent_folder <- "./scRNA_processing/"
setwd(analysis_parent_folder)
path_to_metadata <- "./rna_metadata.csv"


# Define sets and locations of files for scRNA processing
scRNA_data_path <- "./data/"
scRNA_set_names <- c("all_samples_normal", "all_samples_cancer")
scRNA_sets <- list(c("normal1","normal2","normal3","normal5"), 
			c("cancer1","cancer2","cancer3","cancer4","cancer5","cancer6","cancer7")) 

# Define functions
seurat_standard_normalize_and_scale <- function(rna_data, cluster, cluster_resolution){
	# rna_data is seurat object, 
	rna_data <- NormalizeData(rna_data, normalization.method = "LogNormalize", scale.factor = 10000)
	rna_data <- FindVariableFeatures(rna_data, selection.method = "vst", nfeatures = 2000)
	all.genes <- rownames(rna_data)
	rna_data <- ScaleData(rna_data, features = all.genes)
	rna_data <- RunPCA(rna_data, features = VariableFeatures(object = rna_data))
	if (cluster){
		rna_data <- FindNeighbors(rna_data, dims = 1:20)
		rna_data <- FindClusters(rna_data, resolution = cluster_resolution)
	}
	rna_data <- RunUMAP(rna_data, dims = 1:20)
	return(rna_data)
}

make_seurat_object_and_doublet_removal <- function(data_directory, project_name){
	# function for basic seurat based qc and doubletfinder based doublet removal
	rna_data.data <- Read10X(data.dir = data_directory)
	currentSample <- CreateSeuratObject(counts = rna_data.data, project = project_name, min.cells = 3, min.features = 40)
	currentSample[["percent.mt"]] <- PercentageFeatureSet(currentSample, pattern = "^MT-")

	# qc plot-pre filtering
	pdf(paste0("./qc_plots_", project_name, "_prefiltered.pdf"))
	print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05))
	dev.off()
	pdf(paste0("./qc_plots_", project_name, "_prefiltered_no_points.pdf"))
	print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
	dev.off()

	# filter everything to 400 unique genes/cell
	currentSample <- subset(currentSample, subset = nFeature_RNA > 400 & nFeature_RNA < 4000)
	
	# Normalize and make UMAP
	currentSample <- seurat_standard_normalize_and_scale(currentSample, FALSE)

	# Run doublet finder
	nExp_poi <- round(0.08*length(currentSample@meta.data$orig.ident)*length(currentSample@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
	seu_rna_data <- doubletFinder_v3(currentSample, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	print(head(seu_rna_data@meta.data))
	
	# rename columns
	seu_rna_data$doublet.class <- seu_rna_data[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]]
	seu_rna_data[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]] <- NULL
	pann <- grep(pattern="^pANN", x=names(seu_rna_data@meta.data), value=TRUE)
	seu_rna_data$pANN <- seu_rna_data[[pann]]
	seu_rna_data[[pann]] <- NULL

	# plot pre and post doublet finder results
	pdf(paste0("./UMAP_pre_double_removal", project_name, ".pdf"))
	print(DimPlot(seu_rna_data, reduction = "umap", group.by = "doublet.class", cols = c("#D51F26", "#272E6A")))
	dev.off()
	seu_rna_data <- subset(seu_rna_data, subset = doublet.class != "Doublet")
	pdf(paste0("./UMAP_post_double_removal", project_name, ".pdf"))
	print(DimPlot(seu_rna_data, reduction = "umap", cols = c("#D51F26")))
	dev.off()

	# Remove extra stuff and return filtered Seurat object
	seu_rna_data <- DietSeurat(seu_rna_data, counts=TRUE, data=TRUE, scale.data=FALSE, assays="RNA")
	return(seu_rna_data)
}

seurat_qc_plots <- function(rna_data, sample_name){
	# Make some basic qc plots
	pdf(paste0("./seurat_nFeature_plots_", sample_name, ".pdf"), width = 40, height = 15)
	print(VlnPlot(rna_data, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.2))
	dev.off()

	pdf(paste0("./seurat_nCount_plots_", sample_name, ".pdf"), width = 40, height = 15)
	print(VlnPlot(rna_data, features = c("nCount_RNA"), ncol = 1, pt.size = 0.2))
	dev.off()

	pdf(paste0("./seurat_pMT_plots_", sample_name, ".pdf"), width = 40, height = 15)
	print(VlnPlot(rna_data, features = c("percent.mt"), ncol = 1, pt.size = 0.2))
	dev.off()
}

##############################################################################################################################
#Create seurat objects for individual 10x runs, run doblet finder and filter most likely doublets, merge into a seurat object containing all samples#
##############################################################################################################################

	setwd(individual_qc_and_dublet_plot_location)
	for (j in 1:length(scRNA_set_names)){
		samples <- scRNA_sets[[j]]
		print(paste0(scRNA_data_path, samples[1], "/outs/filtered_feature_bc_matrix/"))
		data_directory <- paste0(scRNA_data_path, samples[1], "/outs/filtered_feature_bc_matrix/")
		sample1 <- make_seurat_object_and_doublet_removal(data_directory, samples[1])
		seu_list <- c()
		for (i in 2:length(samples)){
			data_directory <- paste0(scRNA_data_path, samples[i], "/outs/filtered_feature_bc_matrix/")
			seu_list <- c(seu_list, make_seurat_object_and_doublet_removal(data_directory, samples[i]))
		}
		current_merge <- merge(sample1, y = seu_list, add.cell.ids = samples, project = scRNA_set_names[j])
		if (j==1){
			rna_data <- current_merge
		} else if (j>1){
			rna_data <- merge(rna_data, y = current_merge, project = "full_rna_data_project")
		}
	}
	setwd(analysis_parent_folder)
	rna_data[["percent.mt"]] <- PercentageFeatureSet(rna_data, pattern = "^MT-")
	saveRDS(rna_data, "initial_full_rna_data_proj_seurat.rds")
	
 

##############################################################################################################################
# QC..........................
##############################################################################################################################

	# create and set working directory to save qc plots
	if (!dir.exists(paste0(analysis_parent_folder, "all_samples_qc_plots"))){
		dir.create(paste0(analysis_parent_folder, "all_samples_qc_plots"))
	}
	setwd(paste0(analysis_parent_folder, "all_samples_qc_plots"))

	# make the standard seurat qc plots
	seurat_qc_plots(rna_data, sample_name)

	# Now subset the project (if not done already)
	rna_data <- subset(rna_data, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 5 & nCount_RNA < 10000)

	# leave the qc directory
	setwd(analysis_parent_folder)



##############################################################################################################################
# Add metadata............................................................................................................................#
##############################################################################################################################

	metadata <- read.table(path_to_metadata, header = TRUE, sep = ",", stringsAsFactors=FALSE)
	
	# remove atac column
	metadata <- metadata[,colnames(metadata)[2:28]]
	colnames(metadata) <- c("Sample", colnames(metadata)[2:27])
	metadata <- metadata[metadata$Sample != "",]

	meta_data_types <- colnames(metadata)
	for (i in 2:length(meta_data_types)){
		identities <- rna_data[['orig.ident']]
		for (j in 1:length(metadata$Sample)){
			identities[identities==metadata$Sample[j]] <- metadata[j,meta_data_types[i]]
		}
		rna_data <- AddMetaData(rna_data, identities$orig.ident, col.name = meta_data_types[i])
	}


##############################################################################################################################
#Normalize and scale data............................................................................................................................#
##############################################################################################################################

	rna_data <- seurat_standard_normalize_and_scale(rna_data, TRUE, 1.0)


##############################################################################################################################
#Plot UMAPs............................................................................................................................#
##############################################################################################################################

	# Plot by clustering, sample, and disease state
	rna_data <- FindClusters(rna_data, resolution = 0.5)
	pdf(paste0("./UMAPclustering" , ".pdf"), onefile=F)
	print(DimPlot(rna_data, reduction = "umap", cols = paletteDiscrete(values = unique(rna_data@meta.data$seurat_clusters), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_samples.pdf"), width = 12, onefile=F)
	print(DimPlot(rna_data, reduction = "umap", group.by = "orig.ident", cols = paletteDiscrete(values = unique(rna_data@meta.data$orig.ident), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_donor.pdf"), width = 12, onefile=F)
	print(DimPlot(rna_data, reduction = "umap", group.by = "Donor", cols = paletteDiscrete(values = unique(rna_data@meta.data$Donor), set = "stallion", reverse = FALSE)))
	dev.off()

	pdf(paste0("./UMAP_disease_state.pdf"), width = 6.5, onefile=F)
	print(DimPlot(rna_data, reduction = "umap", group.by = "DiseaseState",
		cols = c("#D51F26", "#D7CEC7", "#89288F", "#D7CEC7")) + theme_ArchR())
	dev.off()

	saveRDS(rna_data, "clustered_full_rna_data_proj_seurat.rds")



##############################################################################################################################
# ID Cluster Markers..........................................................................................................
##############################################################################################################################

	# Find cluster specific markers
	rna_data.markers <- FindAllMarkers(rna_data, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.25, max.cells.per.ident = 250)
	top10 <- rna_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)



##############################################################################################################################
#Initial cluster identification...............................................................................................
##############################################################################################################################

	# Rough initial cluster identification
	new.cluster.ids <- c("Epithelial", #0
		"Epithelial",
		"Epithelial",
		"Epithelial",
		"Epithelial",
		"Epithelial",#5
		"Epithelial",
		"Epithelial",
		"Immune",
		"Stromal",
		"Immune",#10
		"Epithelial",
		"Stromal",
		"Epithelial",
		"Immune",
		"Epithelial",#15
		"Stromal",
		"Stromal",
		"Immune",
		"Tuft",
		"Epithelial",#20
		"Immune",
		"Epithelial",
		"Enteroendocrine",
		"Stromal",
		"Immune"#25
		)

	identities <- as.character(rna_data@meta.data$seurat_clusters)
	for (i in 0:length(new.cluster.ids)){
		identities[identities==as.character(i)] <- new.cluster.ids[i+1]
	}
	rna_data <- AddMetaData(rna_data, identities, col.name = "CellTypeInitial")
	pdf(paste0("./UMAP_cell_type_initial.pdf"), width = 12)
	AugmentPlot(DimPlot(rna_data, reduction = "umap", group.by = "CellTypeInitial"))
	dev.off()

	pdf(paste0("./UMAP_cell_type_initial_v2.pdf"), width = 12)
	plot = DimPlot(rna_data, reduction = "umap", group.by = "CellTypeInitial")
	plot = LabelClusters(plot = plot, id = "CellTypeInitial")
	AugmentPlot(plot)
	dev.off()

	# Save ids of immune, epithelial, and stromal cells
	immune_clusters <- c(8,10,14,18,21,25)
	stromal_clusters <- c(9,12,16,17,24)
	epithelial_clusters <- c(0:7,11,13,15,19,20,22,23)
	epithelial_specialized_clusters <- c(19,23)
	epithelial_nonspecialized_clusters <- c(0:7,12,14,20,21,25)

	# Subset to make immune, stromal, and epithelial projects
	rna_data_immune <- DietSeurat(subset(rna_data, subset = seurat_clusters %in% immune_clusters))
	rna_data_stromal <- DietSeurat(subset(rna_data, subset = seurat_clusters %in% stromal_clusters))
	rna_data_epitehlial <- DietSeurat(subset(rna_data, subset = seurat_clusters %in% epithelial_clusters))
	rna_data_specialized_epitehlial <- DietSeurat(subset(rna_data, subset = seurat_clusters %in% epithelial_specialized_clusters))
	rna_data_nonspecialized_epitehlial <- DietSeurat(subset(rna_data, subset = seurat_clusters %in% epithelial_nonspecialized_clusters))

	saveRDS(rna_data_immune, file = "rna_data_immune_all_samples_initial.rds")
	saveRDS(rna_data_stromal, file = "rna_data_stromal_all_samples_initial.rds")
	saveRDS(rna_data_epitehlial, file = "rna_data_epithelial_all_samples_initial.rds")
	saveRDS(rna_data_specialized_epitehlial, file = "rna_data_specialized_epitehlial_all_samples_initial.rds")
	saveRDS(rna_data_nonspecialized_epitehlial, file = "rna_data_nonspecialized_epitehlial_all_samples_initial.rds")




