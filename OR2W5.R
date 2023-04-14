plan("multisession", workers = 10)

#OR Samples

#SS11
ss11.data<- Read10X(data.dir = "OR2W5/SS11/raw_feature_bc_matrix:/")
ss11_seurat_unfiltered<- CreateSeuratObject(counts = ss11.data, project = "SS11", min.cells = 3, min.features = 200)
ss11_seurat_unfiltered
ss11_seurat_unfiltered[["percent.mt"]] <- PercentageFeatureSet(ss11_seurat_unfiltered, pattern = "^MT-")
VlnPlot(ss11_seurat_unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(ss11_seurat_unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ss11_seurat_unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#SS12
ss12.data<- Read10X(data.dir = "OR2W5/SS12/raw_feature_bc_matrix:/")
ss12_seurat_unfiltered<-pbmc <- CreateSeuratObject(counts = ss12.data, project = "SS12", min.cells = 3, min.features = 200)
ss12_seurat_unfiltered
ss12_seurat_unfiltered[["percent.mt"]] <- PercentageFeatureSet(ss12_seurat_unfiltered, pattern = "^MT-")
VlnPlot(ss12_seurat_unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(ss1_seurat_unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ss1_seurat_unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#Merging SS11 and SS12 to see if we need integration
OR_merged<-merge(ss11_seurat_unfiltered, y = ss12_seurat_unfiltered, add.cell.ids = c("ss11", "ss12"), project = "OR2W5")
OR_merged[["percent.mt"]] <- PercentageFeatureSet(OR_merged, pattern = "^MT-")
OR_merged
VlnPlot(OR_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="orig.ident")+geom_hline(yintercept = 20)
OR_merged_filtered <- subset(OR_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
OR_merged_filtered


#Testing SCT transform (not using)
OR_SCT <- SCTransform(OR_merged_filtered, vars.to.regress = "percent.mt", verbose = FALSE)
OR_SCT <- RunPCA(OR_SCT, verbose = FALSE)
OR_SCT <- RunUMAP(OR_SCT, dims = 1:35, verbose = FALSE)

OR_SCT <- FindNeighbors(OR_SCT, dims = 1:35, verbose = FALSE)
OR_SCT <- FindClusters(OR_SCT, verbose = FALSE)
DimPlot(OR_SCT,label=TRUE) + NoLegend()
OR_SCT_Azimuth <- RunAzimuth(OR_SCT, reference = "bonemarrowref")
p1 <- DimPlot(OR_SCT_Azimuth, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) 
p1
rm(OR_SCT)
rm(OR_SCT_Azimuth)


#Regular normalization
OR_merged_normalized<-NormalizeData(OR_merged_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
OR_merged_hvg <- FindVariableFeatures(OR_merged_normalized, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(OR_merged_hvg), 10)
plot1 <- VariableFeaturePlot(OR_merged_hvg)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
plot1 + plot2


#Cell cycle assignment 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
OR_cell_cycle <- CellCycleScoring(OR_merged_hvg, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


#scaling and regressing cell cycle, MT, UMI out 

OR_cell_cycle$CC.Difference <- OR_cell_cycle$S.Score - OR_cell_cycle$G2M.Score
OR_cell_cycle <- ScaleData(OR_cell_cycle, vars.to.regress = c("CC.Difference","percent.mt","nCount_RNA"), features = rownames(OR_cell_cycle))

#PCA
OR_cell_cycle <- RunPCA(OR_cell_cycle, features = VariableFeatures(object = OR_cell_cycle))

DimPlot(OR_cell_cycle, reduction = "pca")
VlnPlot(object = OR_cell_cycle, features = "PC_1", group.by = "orig.ident", pt.size = .1)
#NO harmony data integration


OR_merged_scaled <- JackStraw(OR, dims=100,num.replicate = 100)
OR_merged_scaled <- ScoreJackStraw(OR_merged_scaled, dims = 1:50)
JackStrawPlot(OR_merged_scaled, dims = 1:50)
ElbowPlot(OR_merged_scaled,ndims=50)

#Clustering
OR_merged_scaled <- FindNeighbors(OR_merged_scaled, dims = 1:80)
OR_merged_scaled <- FindClusters(OR_merged_scaled, resolution = 0.4)

#UMAP
OR_merged_scaled <- RunUMAP(OR_merged_scaled, dims = 1:40)

DimPlot(OR_merged_scaled, reduction = "umap",label=TRUE)

OR_AZ<-RunAzimuth(OR_merged_scaled,reference="bonemarrowref")
p1 <- DimPlot(OR_AZ, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) 
p1

write_rds(OR_merged_scaled,"OR2W5/OR_SEURAT.rds")
OR<-read_rds("OR2W5/OR_SEURAT.rds")
OR


#Azimuth
OR_merged_scaled_Azimuth <- RunAzimuth(OR_merged_scaled, reference = "bonemarrowref")
p1 <- DimPlot(OR_merged_scaled_Azimuth, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) 
p1
OR_azimuth_PBMC <- RunAzimuth(OR_merged_scaled, reference = "pbmcref")
p2 <- DimPlot(OR_azimuth_PBMC, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) 
p2
p3 <- FeaturePlot(OR_merged_scaled, features = "CDK6")
p3

#DE 
OR_markers <- FindAllMarkers(OR_merged_scaled, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
OR_markers
OR_markers |>
  group_by(cluster) |>
  slice_max(n = 10, order_by = avg_log2FC) |> print(n=130)
OR_markers |>
  group_by(cluster) |>
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(OR_merged_scaled, features = top10$gene) + NoLegend()+theme(axis.text.y = element_text(size = 3))
write.csv(OR_markers,"OR2W5//OR_marekrs.csv")



#Cell cycle assignment 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
OR_cell_cycle <- CellCycleScoring(OR_merged_scaled, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(OR_cell_cycle, reduction = "umap")

























#mi-R126 High subset 
OR_miR126= AddModuleScore(
  OR_merged_scaled,
  features=list(c("FAM30A","EGFL7","GUCY1A1","ANGPT1","MSI2","CXXC5","TSC22D1","SEPT6","DDAH2","ITM2A","IGHM","AKR1C3")),
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = "miR126_high",
  seed = 1,
  search = FALSE,
)
names(x = OR_miR126[[]])
FeaturePlot(object = OR_miR126, features = 'miR126_high1')+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


#VJAY LTHSC SIGNITURE
OR_HSC= AddModuleScore(
  OR_merged_scaled,
  features=list(c("CD34", "HLF", "CRHBP")),
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = "lthsc",
  seed = 1,
  search = FALSE,
)
names(x = OR_HSC[[]])
FeaturePlot(object = OR_HSC, features = 'lthsc1')









