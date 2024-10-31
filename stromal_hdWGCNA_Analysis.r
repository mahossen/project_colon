#Import libraries
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

theme_set(theme_cowplot())

set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# load dataset
seurat_obj_stromal<- UpdateSeuratObject(sub_integrated_atac_rna)

#Checking Initial cell types.

pdf(file = "/home/malf/alim/integtation/InitialCellType.pdf", width = 6, height = 5) 
p <- DimPlot(seurat_obj_stromal,reduction = "umap", group.by='InitialCellType', 
             label=TRUE, label.size = 5, pt.size = 0.7,
             cols = paletteDiscrete(values = unique(seurat_obj_stromal@meta.data$InitialCellType), set = "stallion", reverse = FALSE)) +
  umap_theme() + ggtitle('Major CellType') 
p
dev.off()

#############################################################
# Set up Seurat object for WGCNA.............................
#############################################################


seurat_obj_stromal <- SetupForWGCNA(
  seurat_obj_stromal,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "stromal Celltype" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj_stromal <- MetacellsByGroups(
  seurat_obj_stromal,
  group.by = 'InitialCellType', # specify the columns in seurat_obj_stromal@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 20, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'InitialCellType' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj_stromal <- NormalizeMetacells(seurat_obj_stromal)

metacell_obj <- GetMetacellObject(seurat_obj_stromal)

#############################################################################
#Co-expression network analysis by hdWGNCA...................................
#############################################################################

seurat_obj_stromal <- SetDatExpr(
  seurat_obj_stromal,
  group_name = "stromal", # the name of the group of interest in the group.by column
  group.by='InitialCellType', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

####################################################
Select soft-power threshold ........................
####################################################

# Test different soft powers:
seurat_obj_stromal <- TestSoftPowers(seurat_obj_stromal,  networkType = 'signed')# you can also use "unsigned" or "signed hybrid"


# plot the results:
pdf(file = "SoftPowers.pdf", width = 8, height = 5) 

plot_list <- PlotSoftPowers(seurat_obj_stromal)
# assemble with patchwork
wrap_plots(plot_list, ncol=2)
dev.off()
#Power Table construction
power_table <- GetPowerTable(seurat_obj_stromal)
head(power_table)

#########################################################################
#Construct co-expression network.........................................
#########################################################################

# construct co-expression network:
seurat_obj_stromal <- ConstructNetwork(seurat_obj_stromal, soft_power=4, setDatExpr=FALSE,tom_name = 'stromal') # name of the topoligical overlap matrix written to disk


#visualize the WGCNA dendrogram,
pdf(file = "stromal Dendogram.pdf", width = 8, height = 5) 
PlotDendrogram(seurat_obj_stromal, main='Stromal Dendrogram')
dev.off()


#The grey module should be ignored for all downstream analysis and interpretation.

#########################################################################################
#Module Eigengenes and Connectivity
#########################################################################################

#Compute harmonized module eigengenes

#ScaleData first or else harmony throws an error:
seurat_obj_stromal <- ScaleData(seurat_obj_stromal, features=VariableFeatures(seurat_obj_stromal))

# compute all MEs in the full single-cell dataset
seurat_obj_stromal <- ModuleEigengenes(seurat_obj_stromal, group.by.vars="DiseaseState")

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj_stromal)

# module eigengenes:
MEs <- GetMEs(seurat_obj_stromal, harmonized=FALSE)

#Compute module connectivity
# compute eigengene-based connectivity (kME):
seurat_obj_stromal <- ModuleConnectivity(seurat_obj_stromal, group.by = 'InitialCellType', group_name = 'stromal')

# rename the modules
#seurat_obj_stromal <- ResetModuleNames(seurat_obj_stromal,new_name = "stromal-M")


# plot genes ranked by kME for each module
pdf(file = "stromal_hub_gene_plot.pdf", width = 8, height = 5) 
p <- PlotKMEs(seurat_obj_stromal,ncol = 3, n_hubs = 10)
p
dev.off()


# get the module assignment table:
modules <- GetModules(seurat_obj_stromal)

head(modules[,1:6]) # show the first 6 columns:
write.csv(modules, "moduless_stromal.csv")


# get hub genes
hub_df <- GetHubGenes(seurat_obj_stromal, n_hubs = 100)
head(hub_df)
write.csv(hub_df, "100_hub_genes_stromal.csv")
##This wraps up the critical analysis steps for hdWGCNA, so remember to save your output.

saveRDS(seurat_obj_stromal, file='stromal_hdWGCNA_object.rds')


#Compute hub gene signature scores

# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
seurat_obj_stromal <- ModuleExprScore(seurat_obj_stromal,  n_genes = 25,  method='Seurat')

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj_stromal <- ModuleExprScore(seurat_obj_stromal, n_genes = 25, method='UCell')


############################################################################
#Basic Visualization........................................................
############################################################################

#Module Feature Plots

# make a featureplot of hMEs for each module
pdf(file = "stromal_ModuleFeaturePlot.pdf", width = 6, height = 5) 
plot_list <- ModuleFeaturePlot(
  seurat_obj_stromal,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

wrap_plots(plot_list, ncol=3)
dev.off()

#We can also plot the hub gene signature score using the same function:
# make a featureplot of hub scores for each module
pdf(file = "stromal_ModuleFeature_score.pdf", width = 6, height = 5) 
plot_list <- ModuleFeaturePlot(
  seurat_obj_stromal,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)
wrap_plots(plot_list, ncol=3)
dev.off()


# plot module correlagram
pdf(file = "stromal_ModuleCorrelogram.pdf", width = 6, height = 5) 
ModuleCorrelogram(seurat_obj_stromal)
dev.off()

###################################################################################
##Seurat plotting ........................................................
###################################################################################

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj_stromal, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj_stromal@meta.data <- cbind(seurat_obj_stromal@meta.data, MEs)

# plot with Seurat's DotPlot function
pdf(file = "stromal_DotPlot.pdf", width = 10, height = 2) 
p <- DotPlot(seurat_obj_stromal, features=mods, group.by = 'InitialCellType')
p
dev.off()

# flip the x/y axes, rotate the axis labels, and change color scheme:
pdf(file = "stromal_DotPlot_2.pdf", width = 8, height = 6) 
p <- DotPlot(seurat_obj_stromal, features=mods, group.by = 'InitialCellType')
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')
p
dev.off()

# Plot  specific mudule hME using Seurat VlnPlot function
pdf(file = "stromal_Vin_Blue_module.pdf", width = 8, height = 6) 
p <- VlnPlot(
  seurat_obj_stromal,
  features = 'blue',
  group.by = 'InitialCellType',
  pt.size = 0 # don't show actual data points
)

# add box-and-whisker plots on top:
p <- p + geom_boxplot(width=.25, fill='white')
# change axis labels and remove legend:
p <- p + xlab('') + ylab('hME') + NoLegend()
p
dev.off()

#################################################################################################
#Individual module network plots.................................................................
#################################################################################################

ModuleNetworkPlot(
  seurat_obj_stromal,
  outdir = 'ModuleNetworks'
)

ModuleNetworkPlot(
  seurat_obj_stromal, 
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
  seurat_obj_stromal,
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

seurat_obj_stromal <- RunModuleUMAP(
  seurat_obj_stromal,
  n_hubs = 5, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj_stromal)

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
  seurat_obj_stromal,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=5 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)
dev.off()



#######################################################################################
#Module Trait Correlation.............................................................
#######################################################################################

seurat_obj_stromal$DiseaseState <- as.factor(seurat_obj_stromal$DiseaseState)
seurat_obj_stromal$FAP <- as.factor(seurat_obj_stromal$FAP)

seurat_obj_stromal <- ModuleTraitCorrelation(
  seurat_obj_stromal,
  traits = c("DiseaseState", "FAP"),
  group.by='InitialCellType'
)

# get the mt-correlation results
mt_cor <- GetModuleTraitCorrelation(seurat_obj_stromal)
names(mt_cor)
names(mt_cor$cor)

pdf(file = "stromal_ModuleTraitCorrelation_FDR.pdf", width = 8, height = 5)
p<-PlotModuleTraitCorrelation(
  seurat_obj_stromal, 
  label = 'fdr',
  label_symbol = 'numeric',
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
  low_color = 'purple',
  plot_max = 0.2,
  combine=TRUE
)
p
dev.off()


pdf(file = "stromal_ModuleTraitCorrelation_Pval.pdf", width = 8, height = 5)
p<-PlotModuleTraitCorrelation(
  seurat_obj_stromal, 
  label = 'pval',
  label_symbol = 'numeric',
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
  low_color = 'purple',
  plot_max = 0.2,
  combine=TRUE
)
p
dev.off()



##############################################################################################
#Functional and pathway Enrichment analysis...................................................
##############################################################################################

library(enrichR)

dbs <-c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023','KEGG_2021_Human')

# compute GO terms:
seurat_obj_stromal <- RunEnrichr(seurat_obj_stromal, dbs=dbs)

enrichr_df <- GetEnrichrTable(seurat_obj_stromal) %>% subset(P.value < 0.05)
write.table(enrichr_df, quote=FALSE, sep='\t', row.names=FALSE, file=paste0("/home/malf/alim/integtation/results/stromal/", 'stromal_enrichr.tsv'))


# make GO term plots:
EnrichrBarPlot(
  seurat_obj_stromal,
  outdir = "/home/malf/alim/integtation/results/stromal/",
  n_terms = 25, plot_size = c(4,16),
  logscale=TRUE
)



# enrichr dotplot
pdf( 'stromal_BP_dotplot.pdf', width=8, height=10, useDingbats=FALSE)
p<-EnrichrDotPlot(
  seurat_obj_stromal,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2023", # this has to be one of the lists we used above!!!
  n_terms=2 # number of terms for each module
)
p
dev.off()

# enrichr dotplot for CC
pdf( 'stromal_CC_dotplot.pdf', width=8, height=10, useDingbats=FALSE)
p <- EnrichrDotPlot(
  seurat_obj_stromal,
  database = "GO_Cellular_Component_2023", n_terms=3,
  outdir = "/home/malf/alim/integtation/results/stromal/",
  break_ties=TRUE
)
p
dev.off()

# enrichr dotplot for CC
pdf( 'stromal_MP_dotplot.pdf', width=8, height=10, useDingbats=FALSE)
p <- EnrichrDotPlot(
  seurat_obj_stromal,
  database = "GO_Molecular_Function_2023", n_terms=3,
  outdir = "/home/malf/alim/integtation/results/stromal/",
  break_ties=TRUE
)
p
dev.off()

# enrichr dotplot for CC
pdf( 'stromal_KEGG_dotplot.pdf', width=8, height=10, useDingbats=FALSE)
p <- EnrichrDotPlot(
  seurat_obj_stromal,
  database = "KEGG_2021_Human", n_terms=3,
  outdir = "/home/malf/alim/integtation/results/stromal/",
  break_ties=TRUE
)
p
dev.off()

##############################################################################################
#Marker gene overlap analysis.................................................................
###############################################################################################

#compute cell-type marker genes with Seurat

Idents(seurat_obj_stromal) <- seurat_obj_stromal$InitialCellType
markers <- Seurat::FindAllMarkers(
  seurat_obj_stromal,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj_stromal,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)


# overlap barplot, produces a plot for each cell type
pdf( 'stromal_marker_gene_overlap_barplot.pdf', width=6, height=3, useDingbats=FALSE)
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

############################################################################################
#TF motif Analysis.........................................................................
############################################################################################


# packages for TF motif analysis
library(JASPAR2020)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)

# get the pfm from JASPAR2020 using TFBSTools
pfm_core <- TFBSTools::getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# run the motif scan with these settings for the mouse dataset
seurat_obj_stromal <- MotifScan(
  seurat_obj_stromal,
  species_genome = 'hg38',
  pfm = pfm_core,
  EnsDb = EnsDb.Hsapiens.v86
)
dim(GetMotifMatrix(seurat_obj_stromal))


# TF target genes
target_genes <- GetMotifTargets(seurat_obj_stromal)

# check target genes for one TF:
head(target_genes$LAMA2)
