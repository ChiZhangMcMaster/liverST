install.packages(c("usethis", "devtools"))
usethis::create_package("spatialPipelineR")
usethis::use_git()
usethis::use_readme_md()

usethis::use_git_config(
  user.name = "Chi Zhang",
  user.email = "zhanc183@mcmaster.ca"
)
usethis::use_git()
library(spatialPipelineR)
devtools::install()
ls("package:spatialPipelineR")
usethis::use_package("Seurat")
devtools::document()
devtools::install()
library(spatialPipelineR)

liver2 <- load_spatial(
  data_dir = "/Users/yangcheng/Desktop/ExciseProject/Sam28",
  slice = "condition2a"
)

devtools::install_github("ChiZhangMcMaster/liverST")
library(liverST)


remove.packages("spatialPipelineR")
library(spatialPipelineR)
find.package("spatialPipelineR")

devtools::document()
# Remove old package if it exists
remove.packages("spatialPipelineR")

# Install fresh from local source
devtools::install()
library(liverST)
ls("package:liverST")  # Should show all your functions
usethis::use_git()
# Commit message example:
# "Rename package to liverST"
usethis::use_github()  # if not already linked


########After name changed
devtools::install_github("ChiZhangMcMaster/liverST")
library(liverST)
devtools::document()
liver2 <- load_spatial(
  data_dir = "/Users/yangcheng/Desktop/ExciseProject/Sam28",
  slice = "condition2a"
)

devtools::document()
liver2 <- setup_liver(liver2, "condition1a")
head(liver1)

devtools::document()
plot_liver_qc_feature(liver1)
plot_liver_qc_count(liver1)

##Filtering kept as it is, so they can change?

options(future.globals.maxSize = 10 * 1024^3)
devtools::document()
liver2SCT <- normalize_liver_sct(liver2)

# for liver2: show counts and features after SCTransform (same code)
## save the object at this point

### 2. Dimension reduction
### 3. Clustering
### 4. Non-linear dimensional reduction (UMAP/t-SNE)
devtools::document()
liver2SCT1 <- run_umap_clustering(liver2SCT, resolution = 1.0)
plot_umap_clustering(liver2SCT1)

liver2SCT1 <- run_tsne_clustering(liver2SCT, resolution = 1.0)
plot_tsne_clustering(liver2SCT1)

#even keep the dimheatmap???
#I think for now don't need a function for elbow plot, it is straight forward enough


## 5. Identification of cluster biomarkers
devtools::document()
liver2_markers <- identify_cluster_biomarkers(liver2SCT1, logfc.threshold = 1, top_n = 10)
biomarker_heatmap(liver2SCT1, liver2_markers) #slot = "data"


genes <- c("Oat", "Cyp2e1", "Slc1a2", "Cyp2c29")

plot_cluster_violin(liver2SCT1, genes)

plot_spatial_features(liver2SCT1, genes)

plot_umap_features(liver2SCT1, genes)


## 6. Cell type annotation
##    a. Automatic annotation using GPT

# IMPORTANT! Assign your OpenAI API key. See Vignette for details
devtools::document()
Sys.setenv(OPENAI_API_KEY = )
liver2SCT2 <- annotate_clusters_gpt(liver2SCT1, liver2_markers,
                                    output_csv = "GPTAnnotation_liv2.csv")
plot_celltype_annotation(liver2SCT2)

##  b. Annotation via deconvolution using a scRNA-seq reference
##  Methods: DWLS (need school wifi)
liver_dwls <- run_dwls_annotation(
  spatial_obj = liver2SCT2,
  reference_obj = readRDS("RefLivcombined.rds")
)

##  b. Annotation via deconvolution using a scRNA-seq reference
##  Methods: Robust Cell Type Decomposition (RCTD)
combined <- readRDS("RefLivcombined.rds")
liver_rctd <- run_rctd_annotation(
  spatial_obj = liver2SCT2,
  reference_obj = readRDS("RefLivcombined.rds")
)

## 6. Cell type annotation
##     c. Manual annotation
####function has been created before
genes <- c("Oat", "Cyp2e1", "Slc1a2", "Cyp2c29")

plot_cluster_violin(liver2SCT1, genes)

plot_spatial_features(liver2SCT1, genes)

plot_umap_features(liver2SCT1, genes)

## 7. DEG analysis, pathway enrichment analysis and Spatially Variable Genes (SVGs) analysis
## Find difference between hepatocytes (0, 1) , here we compare cluster1 to cluster0
devtools::document()
deg_markers <- run_deg_between_clusters(
  seurat_obj = liver2SCT2,
  ident_1 = "1: Hepatocytes 2",
  ident_2 = "0: Hepatocytes 1"
)

head(deg_markers)

plot_enhanced_volcano(
  deg_table = deg_markers,
  logfc_threshold = 0.5,
  padj_threshold = 0.05,
  title = "Hepatocytes2 vs Hepatocytes1"
)

plot_go_cnet(
  deg_table = deg_markers,
  logfc_threshold = 0,
  padj_threshold = 0.05,
  ont = "ALL",
  show_category = 5,
  circular = TRUE
)

# spatially variable genes (SVGs) analysis--pre-annotated anatomical regions
# Spatial heterogeneity
devtools::document()

liver2SCT3 <- annotate_spatial_regions(
  seurat_obj = liver2SCT2,
  image = "condition2a",
  region_names = c("Zone3CV", "Zone1PV")
)
###plot and selection mismatch
plot_spatial_regions(
  seurat_obj = liver2SCT3,
  group_by = "region",
  image = "condition2a",
)

####### DEG analysis between PV and CV
###Did that before



## 8. Integrative analysis across multiple samples or conditions
devtools::document()


# Integrate them using your function
liver.integrated <- integrate_ST_samples(
  object.list = liver.list,
  sample.ids = c("CR1", "EX1", "CR2", "EX2"),
  project.name = "LiverST",
  resolution = 1.2,
  verbose = TRUE
)

DimPlot(liver.integrated, reduction = "umap.harmony", label = TRUE, group.by = c("ident", "orig.ident"))
SpatialDimPlot(liver.integrated, label = TRUE, label.size = 3)


## 8. Integrative analysis across multiple samples or conditions
devtools::install_github("ChiZhangMcMaster/liverST")
library(liverST)
devtools::document()
liver1 <- load_spatial(
  data_dir = "/Users/yangcheng/Desktop/6_STsamples/Sam19",
  sample_id = "control1"
)
head(liver1)

liver2 <- load_spatial(
  data_dir = "/Users/yangcheng/Desktop/6_STsamples/Sam10",
  sample_id = "treat1"
)

liver3 <- load_spatial(
  data_dir = "/Users/yangcheng/Desktop/6_STsamples/Sam18",
  sample_id = "control2"
)

liver4 <- load_spatial(
  data_dir = "/Users/yangcheng/Desktop/6_STsamples/Sam11",
  sample_id = "treat2"
)

liver5 <- load_spatial(
  data_dir = "/Users/yangcheng/Desktop/6_STsamples/Sam39",
  sample_id = "control3"
)

liver6 <- load_spatial(
  data_dir = "/Users/yangcheng/Desktop/6_STsamples/Sam30",
  sample_id = "treat3"
)

####Filtering based on 1-99 percentile
liver1 <- quantile_filter(liver1)
liver2 <- quantile_filter(liver2)
liver3 <- quantile_filter(liver3)
liver4 <- quantile_filter(liver4)
liver5 <- quantile_filter(liver5)
liver6 <- quantile_filter(liver6)

liver1SCT <- normalize_sct(liver1)
liver2SCT <- normalize_sct(liver2)
liver3SCT <- normalize_sct(liver3)
liver4SCT <- normalize_sct(liver4)
liver5SCT <- normalize_sct(liver5)
liver6SCT <- normalize_sct(liver6)

## save the object at this point
saveRDS(liver1SCT, file = "liver1SCT_Prtc.rds")
saveRDS(liver2SCT, file = "liver2SCT_Prtc.rds")
saveRDS(liver3SCT, file = "liver3SCT_Prtc.rds")
saveRDS(liver4SCT, file = "liver4SCT_Prtc.rds")
saveRDS(liver5SCT, file = "liver5SCT_Prtc.rds")
saveRDS(liver6SCT, file = "liver6SCT_Prtc.rds")

library(liverST)
devtools::document()

liver.list <- list(liver1SCT, liver2SCT, liver3SCT, liver4SCT, liver5SCT, liver6SCT)

liver.integrated <- integrate_ST_samples(
  object.list = liver.list,
  resolution = 1.0,
  nfeatures = 3000
)


# Unlabeled UMAP
plot_umap_svg(liver.integrated_4,
              label = T,
              group.by = c("ident", "orig.ident"),
              file = "17_1-clusterMergeIntegr_DimPlotwoLabel_2samp_withlabel",
              width = 15,
              height = 5)

#####could not get cluster 0
plot_spatial_svg(liver.integrated_4,
                 label = TRUE,
                 label.size = 3,
                 file = "17_2-clusterMergeIntegr_SpatialDimPlot",
                 width = 15,
                 height = 6)


liver.list_reduced <- list(liver1SCT, liver2SCT, liver3SCT, liver4SCT)

liver.integrated_4 <- integrate_ST_samples(
  object.list = liver.list_reduced,
  resolution = 0.8,
  nfeatures = 3000
)

liver.integrated_4_markers <- identify_cluster_biomarkers(liver.integrated_4)
top10_liver.integrated_4_markers <- rank_cluster_biomarkers(liver.integrated_4_markers, logfc.threshold = 0.5, padj.threshold = 0.05, top_n = 10, file = "liver_top20_markers.csv")


saveRDS(liver.integrated_4, file = "4liver1_integrated.rds")

####Annotation
Sys.setenv(OPENAI_API_KEY = '')
liver.integrated_4_anno <- annotate_clusters_gpt(liver.integrated_4, top10_liver.integrated_4_markers, model = "gpt-5",
                                    output_csv = "GPTAnnotation_4livers.csv")

liver_markers <- read.csv("liver_top20_markers.csv")
liver.integrated_4_anno <- annotate_clusters_gpt(liver.integrated_4, liver_markers, model = "gpt-5",
                                                 output_csv = "GPTAnnotation_4livers.csv")

# Unlabeled UMAP
plot_umap_svg(liver.integrated_4_anno,
              label = T,
              group.by = c("celltype"),
              file = "GPT_annotated_4livers",
              width = 15,
              height = 5)


## Compare HSC cell features between ctl vs treatment ###########
devtools::document()
liver.integrated_4_anno1 <- prepare_celltype_condition_idents(
  liver.integrated_4_anno,
  control_conditions = c("control1", "control2")
)
head(liver.integrated_4_anno1, n = 2)

deg_markers <- run_deg_between_clusters(
  seurat_obj = liver.integrated_4_anno1,
  ident_1 = "10: Macrophages_Control",
  ident_2 = "10: Macrophages_Treat"
)

head(deg_markers)

plot_enhanced_volcano(
  deg_table = deg_markers,
  logfc_threshold = 0.5,
  padj_threshold = 0.05,
  title = "Macrophage Treat vs Control"
)

plot_go_dot(deg_markers)

####9. Pseudobulk
devtools::document()
pseudo_liv <- perform_pseudobulk(liver.integrated_4_anno1)
head(pseudo_liv@meta.data)
deg_markers <- run_pseudobulk_DE(pseudo_liv, "Macrophages") ###plot all cell types

plot_enhanced_volcano(
  deg_table = deg_markers,
  logfc_threshold = 0.5,
  padj_threshold = 0.05,
  title = "Macrophages Treat vs Control"
)
#####Mark down

## 10. Quantification of cell type composition
devtools::document()
                head(liver.integrated_4_anno1)

selected_types <- c(
  "10: Macrophages",
  "2: Hepatocytes 1",
  "6: Cycling hepatocytes"
)

df_summary <- calculate_celltype_composition(liver.integrated_4_anno1,
  selected_celltypes = selected_types
)

plot_celltype_composition(df_summary)

df_total <- calculate_celltype_composition(liver.integrated_4_anno1)
head(df_total)
plot_celltype_composition(df_total, experimental_group = "Treat")

#######11 CellChat
devtools::document()
spatial_paths <- c(
  control1 = "/Users/yangcheng/Desktop/6_STsamples/Sam19/spatial",
  treat1 = "/Users/yangcheng/Desktop/6_STsamples/Sam10/spatial",
  control2 = "/Users/yangcheng/Desktop/6_STsamples/Sam18/spatial",
  treat2 = "/Users/yangcheng/Desktop/6_STsamples/Sam11/spatial"
)

cellchat_objects <- create_cellchat(liver.integrated_4_anno1,
  spatial_paths = spatial_paths
)

cellchat_objects1 <- prepare_cellchat_calculations(
  cellchat_objects,
  species = "mouse"
)

cellchat_objects2 <- cellchat_compute_prob(
  cellchat_objects1,
  type = "truncatedMean",
  trim = 0.1,
  save_dir = "cellchat_after_CP",
  distance.use = FALSE,
  skip_existing = F
)

cellchat_objects3 <- cellchat_network_analysis(
  cellchat_objects2,
  save_dir = "cellchat_final",
  skip_existing = F
)

###Overall number of interactions and interaction strength
devtools::document()

cellchat_circle_plot(
  cellchat_objects3,
  file_format = "tif",
  save_dir = "cellchat_final/circle_plots"
)

cellchat_heatmap(
  cellchat_objects3,
  file_format = "tif",
  save_dir = "cellchat_final/heatmaps",
  measure = "count"
)

#####Pathway specific analysis (only prob)
devtools::document()
cellchat_circle_plot(
  cellchat_objects3,
  file_format = "tif",
  pathway = "MIF",
  save_dir = "cellchat_final/pathway_circle_plots"
)

cellchat_heatmap(
  cellchat_objects3,
  file_format = "tif",
  pathway = "MIF",
  save_dir = "cellchat_final/pathway_heatmaps"
)

cellchat_spatial_plot(
  cellchat_objects3,
  file_format = "tif",
  pathway = "MIF",
  save_dir = "cellchat_final/pathway_spatial_plots"
)

LR_list <- cellchat_LR_contribution_plot(
  cellchat_objects3,
  file_format = "tif",
  pathway = "MIF",
  save_dir = "cellchat_final/pathway_LR_plots"
)

