# Import libraries
library(Seurat)
library(Signac)
library(tidyverse)
library(cowplot)
library(patchwork)
library(magrittr)
library(igraph)
library(reshape2)
library(ggplot2)
library(WGCNA)
library(hdWGCNA)
library(GeneOverlap)

set.seed(12345)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())


# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# load the dataset
seurat_obj_immune <- UpdateSeuratObject(sub_integrated_atac_rna)

#Checking Initial cell types.

p <- DimPlot(seurat_obj_immune, group.by='InitialCellType', label=TRUE) +
  umap_theme() + ggtitle('Major CellType ') + NoLegend()
p

#######Immune Cell Type ################
pdf(file = "Immune_CellType.pdf", width = 6, height = 5) 
p <- DimPlot(colon_immune,reduction = "umap", group.by='CellType', 
             label.size = 5) + ggtitle('Immune CellType') 
p
dev.off()


#############################################################
# Set up Seurat object for WGCNA.............................
#############################################################

seurat_obj_immune <- SetupForWGCNA(
  seurat_obj_immune,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Immune Celltype" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj_immune <- MetacellsByGroups(
  seurat_obj_immune = seurat_obj_immune,
  group.by = "InitialCellType", # specify the columns in seurat_obj_immune@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'InitialCellType' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj_immune <- NormalizeMetacells(seurat_obj_immune)

metacell_obj <- GetMetacellObject(seurat_obj_immune)

#############################################################################
#Co-expression network analysis by hdWGNCA...................................
#############################################################################

seurat_obj_immune <- SetDatExpr(
  seurat_obj_immune,
  group_name = "Immune", # the name of the group of interest in the group.by column
  group.by='InitialCellType', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

####################################################
Select soft-power threshold ........................
####################################################
# Test different soft powers:
seurat_obj_immune <- TestSoftPowers(seurat_obj_immune,  networkType = 'signed')# you can also use "unsigned" or "signed hybrid"


# plot the results:
pdf(file = "SoftPowers.pdf", width = 8, height = 5) 
plot_list <- PlotSoftPowers(seurat_obj_immune)
wrap_plots(plot_list, ncol=2)
dev.off()

#Power Table construction
power_table <- GetPowerTable(seurat_obj_immune)
head(power_table)

#########################################################################
#Construct co-expression network.........................................
#########################################################################

# construct co-expression network:
seurat_obj_immune <- ConstructNetwork(seurat_obj_immune, soft_power=5, setDatExpr=FALSE,tom_name = 'Immune') # name of the topoligical overlap matrix written to disk

#visualize the WGCNA dendrogram,
pdf(file = "Immune Dendogram.pdf", width = 8, height = 5) 
PlotDendrogram(seurat_obj_immune, main='Immune Dendrogram')
dev.off()

#The grey module should be ignored for all downstream analysis and interpretation.

#Optional: inspect the topoligcal overlap matrix (TOM)
TOM <- GetTOM(seurat_obj_immune)
head(TOM)

#########################################################################################
#Module Eigengenes and Connectivity
#########################################################################################

#Compute harmonized module eigengenes
# ScaleData first or else harmony throws an error:
seurat_obj_immune <- ScaleData(seurat_obj_immune, features=VariableFeatures(seurat_obj_immune))

# compute all MEs in the full single-cell dataset
seurat_obj_immune <- ModuleEigengenes(seurat_obj_immune, group.by.vars="DiseaseState")

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj_immune)

# module eigengenes:
MEs <- GetMEs(seurat_obj_immune, harmonized=FALSE)

#Compute module connectivity
# compute eigengene-based connectivity (kME):
seurat_obj_immune <- ModuleConnectivity(seurat_obj_immune, group.by = 'InitialCellType', group_name = 'Immune')

# rename the modules
seurat_obj_immune <- ResetModuleNames(seurat_obj_immune,new_name = "Immune-M")


# plot genes ranked by kME for each module
pdf(file = "immune_hub_gene_plot.pdf", width = 8, height = 5)
p <- PlotKMEs(seurat_obj_immune, n_hubs = 10, ncol=3)

p
dev.off()


#Getting the module assignment table
modules <- GetModules(seurat_obj_immune)

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(seurat_obj_immune, n_hubs = 25)
head(hub_df)

##This wraps up the critical analysis steps for hdWGCNA, so remember to save your output.

saveRDS(seurat_obj_immune, file='Immune_hdWGCNA_object.rds')


#Compute hub gene signature scores

# compute gene scoring for the top 25 hub genes by kME for each module with Seurat method
seurat_obj_immune <- ModuleExprScore(seurat_obj_immune,  n_genes = 25,  method='Seurat')

# compute gene scoring for the top 25 hub genes by kME for each modulewith UCell method

library(UCell)
seurat_obj_immune <- ModuleExprScore(seurat_obj_immune, n_genes = 25, method='UCell')



############################################################################
#Basic Visualization........................................................
############################################################################

#Module Feature Plots

# make a featureplot of hMEs for each module
pdf(file = "Immune_ModuleFeaturePlot.pdf", width = 6, height = 5) 
plot_list <- ModuleFeaturePlot(
  seurat_obj_immune,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
pdf(file = "Immune_ModuleFeaturePlot.pdf", width = 6, height = 5) 
wrap_plots(plot_list, ncol=3)
dev.off()


# make a featureplot of hub scores for each module
pdf(file = "Immune_ModuleFeature_score.pdf", width = 6, height = 5) 
plot_list <- ModuleFeaturePlot(
  seurat_obj_immune,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=3)
dev.off()


# plot module correlagram
pdf(file = "Immune_ModuleCorrelogram.pdf", width = 6, height = 5) 
ModuleCorrelogram(seurat_obj_immune)
dev.off()

###################################################################################
##Seurat plotting ........................................................
###################################################################################

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj_immune, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj_immune@meta.data <- cbind(seurat_obj_immune@meta.data, MEs)

# plot with Seurat's DotPlot function
pdf(file = "Immune_DotPlot.pdf", width = 10, height = 2) 
p <- DotPlot(seurat_obj_immune, features=mods, group.by = 'InitialCellType')
p
dev.off()

# flip the x/y axes, rotate the axis labels, and change color scheme:
pdf(file = "Immune_DotPlot_2.pdf", width = 8, height = 6) 
p <- DotPlot(seurat_obj_immune, features=mods, group.by = 'InitialCellType')
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')
dev.off()

# Plot INH-M4 hME using Seurat VlnPlot function
p <- VlnPlot(
  seurat_obj_immune,
  features = 'Immune-M5',
  group.by = 'InitialCellType',
  pt.size = 0 # don't show actual data points
)

# add box-and-whisker plots on top:
p <- p + geom_boxplot(width=.25, fill='white')

# change axis labels and remove legend:
p <- p + xlab('') + ylab('hME') + NoLegend()

# plot output
p

#################################################################################################
#Individual module network plots.................................................................
#################################################################################################

ModuleNetworkPlot(
  seurat_obj_immune,
  outdir = 'ModuleNetworks'
)

ModuleNetworkPlot(
  seurat_obj_immune, 
  outdir='ModuleNetworks2', # new folder name
  n_inner = 10, # number of genes in inner ring
  n_outer = 15, # number of genes in outer ring
  n_conns = Inf, # show all of the connections
  plot_size=c(10,10), # larger plotting area
  vertex.label.cex=1 # font size
)

# hubgene network
pdf(file = "hubgene network.pdf", width = 8, height = 8) 

HubGeneNetworkPlot(
  seurat_obj_immune,
  mods = "all",
  n_hubs = 5,
  n_other =20,
  sample_edges = TRUE,
  edge_prop = 0.75,
  return_graph = FALSE,
  edge.alpha = 0.5,
  vertex.label.cex = 1.0,
  hub.vertex.size = 7,
  other.vertex.size = 3,
  wgcna_name = NULL

)
dev.off()

#Applying UMAP to co-expression networks

seurat_obj_immune <- RunModuleUMAP(
  seurat_obj_immune,
  n_hubs = 5, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj_immune)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
    
  ) +
  umap_theme()

##showing the co-expression network, and labeling 5 hub genes in each module
pdf(file = "UMAP_hubgene_co-expression_network.pdf", width = 8, height = 8) 
ModuleUMAPPlot(
  seurat_obj_immune,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=5 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)
dev.off()


#extract module specific genes


# Extract module membership information
module_membership <- GetModuleMembership(seurat_obj_immune)

# Get genes assigned to a specific module (replace "module_color" with the actual color)
module_genes <- rownames(module_membership[module_membership == "blue", ])




#######################################################################################
#Module Trait Correlation.............................................................
#######################################################################################
seurat_obj_immune$DiseaseState <- as.factor(seurat_obj_immune$DiseaseState)
seurat_obj_immune$FAP <- as.factor(seurat_obj_immune$FAP)

seurat_obj_immune <- ModuleTraitCorrelation(
  seurat_obj_immune,
  traits = c("DiseaseState", "FAP"),
  group.by='InitialCellType'
)

# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(seurat_obj_immune)
names(mt_cor)
names(mt_cor$cor)

pdf(file = "immune_ModuleTraitCorrelation_FDR.pdf", width = 8, height = 5)
PlotModuleTraitCorrelation(
  seurat_obj_immune, 
  label = 'fdr',
  label_symbol = 'numeric',
  text_size = 3,
  text_digits = 3,
  text_color = 'white',
  high_color = 'red',
  mid_color = 'grey',
  low_color = 'pink',
  plot_max = 0.2,
  combine=TRUE
)

dev.off()

##############################################################################################
#Functional and pathway Enrichment analysis...................................................
##############################################################################################
library(enrichR)

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','KEGG_2021_Human')

# compute GO terms:
seurat_obj_immune <- RunEnrichr(seurat_obj_immune, dbs=dbs)

enrichr_df <- GetEnrichrTable(seurat_obj_immune) %>% subset(P.value < 0.05)
write.table(enrichr_df, quote=FALSE, sep='\t', row.names=FALSE, file=paste0("/home/malf/alim/integtation/results/immune/", 'Immune_enrichr.tsv'))


# make GO term plots:
EnrichrBarPlot(
  seurat_obj_immune,
  outdir = "./immune_wgcna/",
  n_terms = 25, plot_size = c(4,16),
  logscale=TRUE
)

# enrichr dotplot
p <- EnrichrDotPlot(
  seurat_obj_immune,
  database = dbs[1], n_terms=3,
  outdir = "/home/malf/alim/integtation/results/immune/",
  break_ties=TRUE
)
pdf( 'Immune_BP_dotplot.pdf', width=8, height=10, useDingbats=FALSE)
p
dev.off()

# inspect the enrichr table:
enrichr_df <- GetEnrichrTable(seurat_obj_immune)


# enrichr dotplot
p<-EnrichrDotPlot(
  seurat_obj_immune,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=2 # number of terms for each module
)



##############################################################################################
#Marker gene overlap analysis.................................................................
###############################################################################################

#compute cell-type marker genes with Seurat

Idents(seurat_obj_immune) <- seurat_obj_immune$InitialCellType
markers <- Seurat::FindAllMarkers(
  seurat_obj_immune,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj_immune,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)


# overlap barplot, produces a plot for each cell type
pdf( 'immune_marker_gene_overlap_barplot.pdf', width=6, height=3, useDingbats=FALSE)
plot_list <- OverlapBarPlot(overlap_df, label_size = 3)

# stitch plots with patchwork
wrap_plots(plot_list, ncol=3)
dev.off()


# plot odds ratio of the overlap as a dot plot
pdf( 'stromal_Overlap of modules & cell-type markers.pdf', width=6, height=3, useDingbats=FALSE)
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio'
) +
  ggtitle('Overlap of modules & cell-type markers')
dev.off()
