library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(xlsx)
library(harmony)

#Load 10x data for all 5 samples across 2 independent experiments, and filter before merging

#Sham Expt 1
htn_sham_rawdata <- Read10X("~/Documents/Hypertension_SingleCell/HTN_RawData/Sham")
htn_sham1 <- CreateSeuratObject(counts = htn_sham_rawdata, min.cells = 3, min.features = 200, project = "Sham1")

#Read in file with list of dissociation associated genes from O'Flanagan et al. 2019
dag <- read.csv("Dissociation_genes.csv", sep= ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

#find which DAGs actually overlap with genes present in our datasets
#using "gene" column which contains the 507 upregulated DAGs
overlap_htn_sham <- intersect(rownames(htn_sham1), dag$gene)
#384 overlapping genes out of the 507 total

#find the percentage of transcripts that map to either mitochondiral genes or to DAGs
htn_sham1[["htn_sham1.percent.mito"]] <- PercentageFeatureSet(htn_sham1, pattern = "^mt") 
htn_sham1[["htn_sham1.percent.dag"]] <- PercentageFeatureSet(htn_sham1, features = overlap_htn_sham) 

VlnPlot(htn_sham1, features = c("nFeature_RNA", "nCount_RNA", "htn_sham1.percent.mito", "htn_sham1.percent.dag"), ncol = 4)

plot1 <- FeatureScatter(htn_sham1, feature1 = "nCount_RNA", feature2 = "htn_sham1.percent.mito")
plot2 <- FeatureScatter(htn_sham1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(htn_sham1, feature1 = "htn_sham.percent.dag", feature2 = "htn_sham1.percent.mito")
plot4 <- FeatureScatter(htn_sham1, feature1 = "nCount_RNA", feature2 = "htn_sham1.percent.dag")
CombinePlots(plots = list(plot1, plot2, plot3, plot4))

#Filter low quality, dead or stressed cells and putative doublets
htn_sham1 <- subset(htn_sham1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & htn_sham1.percent.mito < 12 & htn_sham1.percent.dag < 20)

#4 DAY HTN
htn_4d_rawdata <- Read10X("~/Documents/Hypertension_SingleCell/HTN_RawData/Day4_AngII")
htn_4d1 <- CreateSeuratObject(counts = htn_4d_rawdata, min.cells = 3, min.features = 200, project = "Day 4_1")

overlap_htn_4d1 <- intersect(rownames(htn_4d1), dag$gene)
#395 overlapping genes out of the 507 total

htn_4d1[["htn_4d1.percent.mito"]] <- PercentageFeatureSet(htn_4d1, pattern = "^mt") 
htn_4d1[["htn_4d1.percent.dag"]] <- PercentageFeatureSet(htn_4d1, features = overlap_htn_4d1) 

VlnPlot(htn_4d1, features = c("nFeature_RNA", "nCount_RNA", "htn_4d1.percent.mito", "htn_4d1.percent.dag"), ncol = 4)

plot1 <- FeatureScatter(htn_4d1, feature1 = "nCount_RNA", feature2 = "htn_4d1.percent.mito")
plot2 <- FeatureScatter(htn_4d1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(htn_4d1, feature1 = "htn_4d1.percent.dag", feature2 = "htn_4d1.percent.mito")
plot4 <- FeatureScatter(htn_4d1, feature1 = "nCount_RNA", feature2 = "htn_4d1.percent.dag")
CombinePlots(plots = list(plot1, plot2, plot3, plot4))

htn_4d1 <- subset(htn_4d1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & htn_4d1.percent.mito < 20 & htn_4d1.percent.dag < 20)

#28 DAY HTN
htn_28d_rawdata <- Read10X("~/Documents/Hypertension_SingleCell/HTN_RawData/Day28_AngII")
htn_28d1 <- CreateSeuratObject(counts = htn_28d_rawdata, min.cells = 3, min.features = 200, project = "Day 28")

overlap_htn_28d <- intersect(rownames(htn_28d1), dag$gene)
#385 overlapping genes out of the 507 total

htn_28d1[["htn_28d1.percent.mito"]] <- PercentageFeatureSet(htn_28d1, pattern = "^mt") 
htn_28d1[["htn_28d1.percent.dag"]] <- PercentageFeatureSet(htn_28d1, features = overlap_htn_28d) 

VlnPlot(htn_28d1, features = c("nFeature_RNA", "nCount_RNA", "htn_28d.percent.mito", "htn_28d.percent.dag"), ncol = 4)

plot1 <- FeatureScatter(htn_28d1, feature1 = "nCount_RNA", feature2 = "htn_28d1.percent.mito")
plot2 <- FeatureScatter(htn_28d1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(htn_28d1, feature1 = "htn_28d1.percent.dag", feature2 = "htn_28d1.percent.mito")
plot4 <- FeatureScatter(htn_28d1, feature1 = "nCount_RNA", feature2 = "htn_28d1.percent.dag")
CombinePlots(plots = list(plot1, plot2, plot3, plot4))

htn_28d1 <- subset(htn_28d1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & htn_28d1.percent.mito < 20 & htn_28d1.percent.dag < 20)

#Load datasets from second experiment
#Sham Expt 2
htn_sham2_rawdata <- Read10X("~/Documents/HTN_Sham_2020")
htn_sham2 <- CreateSeuratObject(counts = htn_sham2_rawdata, min.cells = 3, min.features = 200, project = "Sham2")

overlap_htn_sham2 <- intersect(rownames(htn_sham), dag$gene)
#412 overlapping genes out of the 507 total

htn_sham2[["percent.mito"]] <- PercentageFeatureSet(htn_sham2, pattern = "^mt") 
htn_sham2[["percent.dag"]] <- PercentageFeatureSet(htn_sham2, features = overlap_htn_sham2) 

VlnPlot(htn_sham2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.dag"), ncol = 4)

plot1 <- FeatureScatter(htn_sham2, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(htn_sham2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(htn_sham2, feature1 = "percent.dag", feature2 = "percent.mito")
plot4 <- FeatureScatter(htn_sham2, feature1 = "nCount_RNA", feature2 = "percent.dag")
CombinePlots(plots = list(plot1, plot2, plot3, plot4))

htn_sham2 <- subset(htn_sham2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito < 30 & percent.dag < 30)

#Day 4 Expt 2
htn_day4_2_rawdata <- Read10X("~/Documents/HTN_Day4_2020")
htn_day4_2 <- CreateSeuratObject(counts = htn_day4_rawdata, min.cells = 3, min.features = 200, project = "Day 4_2")

overlap_htn_day4_2 <- intersect(rownames(htn_day4_2), dag$gene)
#413 overlapping genes out of the 507 total

htn_day4_2[["percent.mito"]] <- PercentageFeatureSet(htn_day4_2, pattern = "^mt") 
htn_day4_2[["percent.dag"]] <- PercentageFeatureSet(htn_day4_2, features = overlap_htn_day4_2) 

VlnPlot(htn_day4_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.dag"), ncol = 4)

plot1 <- FeatureScatter(htn_day4_2, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(htn_day4_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(htn_day4_2, feature1 = "percent.dag", feature2 = "percent.mito")
plot4 <- FeatureScatter(htn_day4_2, feature1 = "nCount_RNA", feature2 = "percent.dag")
CombinePlots(plots = list(plot1, plot2, plot3, plot4))

htn_day4_2 <- subset(htn_day4_2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mito < 25 & percent.dag < 30)

#specify condition in a metadata column for each dataset
htn_sham1[["condition"]] <- "Sham"
htn_sham2[["condition"]] <- "Sham"
htn_4d1[["condition"]] <- "Day 4"
htn_day4_2[["condition"]] <- "Day 4"
htn_28d1[["condition"]] <- "Day 28"

#specify experiment in a metadata column for each dataset
htn_sham1[["experiment"]] <- 1
htn_sham2[["experiment"]] <- 2
htn_4d1[["experiment"]] <- 1
htn_day4_2[["experiment"]] <- 2
htn_28d1[["experiment"]] <- 1

#specify both experiment and condition in a metadata column for each dataset
htn_sham1[["experiment_condition"]] <- "Sham_1"
htn_sham2[["experiment_condition"]] <- "Sham_2"
htn_4d1[["experiment_condition"]] <- "Day4_1"
htn_day4_2[["experiment_condition"]] <- "Day4_2"
htn_28d1[["experiment_condition"]] <- "Day28_1"

#MERGE THE FIVE TOGETHER
htn_merged <- merge(x = htn_sham1, y = c(htn_4d1, htn_28d1, htn_sham2, htn_day4_2))

htn_merged[["htn_merged.percent.mito"]] <- PercentageFeatureSet(htn_merged, pattern = "^mt") 

htn_merged <- SCTransform(htn_merged, vars.to.regress = c("nCount_RNA", "htn_merged.percent.mito"), verbose = TRUE, return.only.var.genes = F)
htn_merged <- RunPCA(htn_merged, verbose = FALSE)
ElbowPlot(htn_merged)

#Running Harmony for batch correction
htn_merged <- RunHarmony(htn_merged, group.by.vars = "experiment", theta = 4, assay.use = "SCT")
htn_merged <- FindNeighbors(object = htn_merged, dims = 1:20, reduction = "harmony")
htn_merged <- FindClusters(object = htn_merged, resolution = 1.2)
htn_merged <- RunUMAP(htn_merged, dims = 1:20, reduction = "harmony")
DimPlot(htn_merged, label = TRUE) 
DimPlot(htn_merged, split.by = "condition") 

#Example of making featureplots
FeaturePlot(htn_merged, features = c("C1qc", "Fcgr1", "Mafb", "Mertk", "Lyve1", "Timd4", "Folr2", "Igf1", "Isg20"), order = TRUE)

#Example of making violin plots
VlnPlot(htn_merged, features = c("C1qc", "Fcgr1", "Mafb", "Mertk", "Lyve1", "Timd4", "Folr2", "Igf1", "Isg20"), pt.size = F)

#Example of finding differentially expressed genes
htn_merged.markers <- FindAllMarkers(htn_merged, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
htn_merged_sig <- htn_merged.markers[which(htn_merged.markers$p_val_adj < 0.05), ]
htn_merged_sig <- htn_merged_sig[!grepl("^mt-", rownames(htn_merged_sig)), ]
htn_merged_sig <- htn_merged_sig[!grepl("^Rp", rownames(htn_merged_sig)), ]
write.xlsx(htn_merged_sig, "htn_merged.xlsx")

#Example of making a heatmap
top30 <- htn_merged_sig %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)
DoHeatmap(htn_merged, cells = WhichCells(htn_merged, downsample = 50, seed = 111), features=top30$gene) + theme(text = element_text(size = 4)) + NoLegend()

#Example of renaming clusters
new.cluster.ids <- c("Macrophages", "Ly6Chi monocytes", "Ly6Clo monocytes", "Proliferating cells") 
names(x = new.cluster.ids) <- levels(x = htn_merged)
htn_merged <- RenameIdents(object = htn_merged, new.cluster.ids) 
DimPlot(object = htn_merged)

#Subsetting data and re-clustering
MFProlif <- subset(htn_merged, idents = c("Macrophages", "Proliferating cells"))
MFProlif <- SCTransform(MFProlif, vars.to.regress = c("nCount_RNA", "htn_merged.percent.mito"), verbose = TRUE, return.only.var.genes = F)
MFProlif <- RunPCA(MFProlif, verbose = FALSE)
ElbowPlot(MFProlif)
MFProlif <- RunHarmony(MFProlif, group.by.vars = "experiment", assay.use = "SCT", theta = 4)
MFProlif <- FindNeighbors(MFProlif, dims = 1:10, verbose = FALSE, reduction = "harmony")
MFProlif <- FindClusters(MFProlif, verbose = FALSE, resolution = 0.3)
MFProlif <- RunUMAP(MFProlif, dims = 1:10, reduction = "harmony", verbose = FALSE, n.neighbors = 20, min.dist = 0.01)
DimPlot(MFProlif, label = TRUE)  

#Example of subsetting state 1 (S1) and finding differentially expressed genes across conditions
S1 <- subset(MFProlif, idents = "S1")
Idents(S1) <- "condition"

#Find DEGs between Day 28 vs sham
S1_D28vsSham.markers <- FindMarkers(S1, ident.1 = "Day 28", ident.2 = "Sham", only.pos = F, min.pct = 0.2, logfc.threshold = 0.2)
S1_D28vsSham_sig <- S1_D28vsSham.markers[which(S1_D28vsSham.markers$p_val_adj < 0.05), ]
S1_D28vsSham_sig <- S1_D28vsSham_sig[!grepl("^mt-", rownames(S1_D28vsSham_sig)), ]
S1_D28vsSham_sig <- S1_D28vsSham_sig[!grepl("^Rp", rownames(S1_D28vsSham_sig)), ]
write.xlsx(S1_D28vsSham_sig, "S1_D28vsSham.xlsx")

#Find DEGs between Day 4 vs sham
S1_D4vsSham.markers <- FindMarkers(S1_D4vsSham, ident.1 = "Day 4", ident.2 = "Sham", only.pos = F, min.pct = 0.2, logfc.threshold = 0.2)
S1_D4vsSham_sig <- S1_D4vsSham.markers[which(S1_D4vsSham.markers$p_val_adj < 0.05), ]
S1_D4vsSham_sig <- S1_D4vsSham_sig[!grepl("^mt-", rownames(S1_D4vsSham_sig)), ]
S1_D4vsSham_sig <- S1_D4vsSham_sig[!grepl("^Rp", rownames(S1_D4vsSham_sig)), ]
write.xlsx(S1_D4vsSham_sig, "S1_D4vsSham.xlsx")

#Object is MFProlif and the 9 states are in the metadata column called "states"
#Example of finding proportional frequency of each state in each condition
table(MFProlif$states, MFProlif$condition)

#Example of finding proportional frequency of each state in each experiment
table(MFProlif$states, MFProlif$condition)

#Example of finding proportional frequency of each state in each condition for both experiments
table(MFProlif$states, MFProlif$experiment_condition)





