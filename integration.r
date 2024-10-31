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
library(cowplot)
seed(12345)


#Data loading
rna_data<-readRDS("./scRNA_processing/rna_data.rds")
atac_data<-readRDS("./scRNA_processing/atac_data.rds")

#Plots of both modalities
p1 <- DimPlot(rna_data, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(atac_data, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
p1 + p2
plot <- (p1 + p2) & xlab("UMAP 1") & ylab("UMAP 2") & theme(axis.title = element_text(size = 18))
ggsave(filename = "./integration/integration.pdf", height = 6, width = 11, plot = plot, quality = 50)

#identifying Anchors between scRNA-seq and scATAC-seq
transfer.anchors <- FindTransferAnchors(
			  reference = rna_data,
			  query = atac_data,
			  reference.assay = "RNA", 
			  query.assay = "ACTIVITY"
			  reduction = 'cca'
			)
			
#Annotate scATAC-seq cells via label transfer
predicted.labels <- TransferData(
			  anchorset = transfer.anchors,
			  refdata = rna_data$celltype,
			  weight.reduction = atac_data[['lsi']],
			  dims = 2:30
			)
			
atac_data <- AddMetaData(object = atac_data, metadata = predicted.labels)

plot1 <- DimPlot(
		  object = rna_data,
		  group.by = 'celltype',
		  label = TRUE,
		  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
		  object = atac_data,
		  group.by = 'predicted.id',
		  label = TRUE,
		  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2
plot <- (plot1 + plot2) & xlab("UMAP 1") & ylab("UMAP 2") & theme(axis.title = element_text(size = 18))
ggsave(filename = "./integration/annotation_atac.pdf", height = 6, width = 11, plot = plot, quality = 50)


############################################################################################
#Co-embedding scRNA-seq and scATAC-seq datasets
#############################################################################################
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the full transcriptome if we wanted to

genes.use <- VariableFeatures(rna_data)
refdata <- GetAssayData(rna_data, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac_data[["lsi"]], dims = 2:30)
atac_data[["RNA"]] <- imputation

coembed <- merge(x = rna_data, y = atac_data)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

P<-DimPlot(coembed, group.by = c("orig.ident", "seurat_annotations"))

ggsave(filename = "./integration/embedding.pdf", height = 6, width = 11, plot = P, quality = 50)

saveRDS(coembed, "./integration/integrated_atac_rna.rds")

