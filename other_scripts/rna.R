library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(harmony)
library(GenomicRanges)
library(SeuratData)
library(SeuratDisk)
library(Matrix)

setwd('/media/leon/Masha/ATAC')

data <- LoadH5Seurat('raw_merged.h5Seurat')

bed_data <- read.table("/media/leon/Masha/ATAC/aggr/outs/filtered_peak_bc_matrix/filtered_peaks.bed", header = FALSE, sep="\t", stringsAsFactors=FALSE)

peak_ids <- paste0(bed_data$V1, "-", bed_data$V2, "-", bed_data$V3)

barcodes = readLines('/media/leon/Masha/ATAC/05.02.signac_barcodes.txt')

data <- subset(data, cells = barcodes)
data <- data[peak_ids, ]

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(data) <- annotations

data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = 'q75')
data <- RunSVD(data)

DefaultAssay(data) <- 'peaks'
data <- RunUMAP(object = data, reduction = 'lsi', dims = 2:30)
data <- FindNeighbors(object = data, reduction = 'lsi', dims = 2:30)
data <- FindClusters(object = data, verbose = FALSE, algorithm = 3, resolution = 1)
DimPlot(object = data, label = TRUE, raster=FALSE) + NoLegend()


gene.activities <- GeneActivity(data)

data[['RNA']] <- CreateAssayObject(counts = gene.activities)

data <- NormalizeData(
  object = data,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(data$nCount_RNA)
)

DefaultAssay(data) <- 'RNA'

FeaturePlot(
  object = data,
  #features = c('GCG', 'INS', 'SST', 'PPY', 'REG1A', 'CFTR', 'PDGFRB', 'CLEC14A', 'CCL3'),
  features = c('GCG'),
  max.cutoff = 'q95',
  raster = F
)

mat <- GetAssayData(object = data, assay = "RNA", slot = "counts")

writeMM(mat, "full_rna_matrix.mtx")
writeLines(rownames(mat), '01.06.full_rna_genes.txt')
writeLines(colnames(mat), '01.06.rna_cells.txt')


# Extract cluster identities
clusters <- Idents(data)

# Get the unique cluster identities
unique_clusters <- unique(clusters)

# Initialize a list to store barcodes for each cluster
barcodes_per_cluster <- list()

# Loop through each cluster and get barcodes
for (cluster in unique_clusters) {
  barcodes_per_cluster[[as.character(cluster)]] <- names(clusters[clusters == cluster])
}

# Determine the maximum number of barcodes in any cluster
max_length <- max(sapply(barcodes_per_cluster, length))

# Make all lists in barcodes_per_cluster the same length by padding with NA
for (cluster in names(barcodes_per_cluster)) {
  length(barcodes_per_cluster[[cluster]]) <- max_length
}

# Convert the list to a data frame
barcodes_df <- as.data.frame(barcodes_per_cluster)

# Write the data frame to a CSV file
write.csv(barcodes_df, file = "cluster_barcodes_csv.csv", row.names = FALSE, na = "")
