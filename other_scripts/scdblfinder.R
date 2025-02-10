library(singleCellTK)
library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Signac)
suppressPackageStartupMessages(library(scDblFinder))

find_doublets <- function(sample){
  
counts <- Read10X_h5(filename = paste0('/media/leon/Masha/ATAC/CR_output/sample', sample, '/outs/filtered_peak_bc_matrix.h5'))

metadata <- read.csv(
  file = paste0('/media/leon/Masha/ATAC/CR_output/sample', sample, '/outs/singlecell.csv'),
  header = TRUE,
  row.names = 1)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = paste0('/media/leon/Masha/ATAC/CR_output/sample', sample, '/outs/fragments.tsv.gz'),
  min.cells = 0,
  min.features = 500
)

data <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

sce <- as.SingleCellExperiment(data)
sce <- scDblFinder(sce, artificialDoublets=1, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures", dbr=0.005*dim(counts)[2]/1000)

all_barcodes <- rownames(colData(sce))

doublet_barcodes <- all_barcodes[sce$scDblFinder.class == "doublet"]

writeLines(doublet_barcodes, paste0('/media/leon/Masha/ATAC/doublets/doublets_', sample, '.txt'))

}

for (i in 61:94) {
  find_doublets(as.character(i))
}
