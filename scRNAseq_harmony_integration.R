library(Seurat)
library(cowplot)
library(harmony)

data <- Read10X(data.dir = "run_cellranger_aggr/aggr_all/outs/filtered_feature_bc_matrix")

mir29ko <- CreateSeuratObject(counts = data, project = "mir29ko_mapped", min.cells = 3, min.features = 200)
mir29ko$state <- substring(names(mir29ko$orig.ident), 18, 18)
mir29ko$state = replace(mir29ko$state, mir29ko$state==1, 'Cre-')
mir29ko$state = replace(mir29ko$state, mir29ko$state==2, 'Cre+')
mir29ko$state = replace(mir29ko$state, mir29ko$state==3, 'Neo')

mir29ko[["percent.mt"]] <- PercentageFeatureSet(mir29ko, pattern = "^mt-")
#VlnPlot(mir29ko, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# filter cells with >15% mitochondria genes
mir29ko <- subset(mir29ko, subset = percent.mt < 15) 

mir29ko <- NormalizeData(mir29ko, normalization.method = "LogNormalize", scale.factor = 10000)
mir29ko <- FindVariableFeatures(mir29ko, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(mir29ko)
mir29ko <- ScaleData(mir29ko, features = all.genes)

# PCA
mir29ko <- RunPCA(mir29ko, features = VariableFeatures(object = mir29ko))


options(repr.plot.height = 2.5, repr.plot.width = 6)
mir29ko <- mir29ko %>% 
    RunHarmony("state", plot_convergence = TRUE)
    
mir29ko <- mir29ko %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()

