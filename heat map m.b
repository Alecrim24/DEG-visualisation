# Load required libraries
library(DESeq2)
library(GenomicFeatures)

# Set file paths
countdata_mb_path <- "C:/Users/bop21gwh/Documents/PhD work/genomics/Hisat2/hisat2_mb.txt"
metadata_path <- "C:/Users/bop21gwh/Documents/PhD work/genomics/Hisat2/metadata.csv"
gtf_file_path <- "C:/Users/bop21gwh/Documents/PhD work/genomics/Hisat2/cleaned_m.b_annotation.gtf"

# Read count data
countdata_mb <- read.table(countdata_mb_path, header = TRUE, sep = "\t")

# Filter count data (removing columns 2-8 as per your original code)
countdata_mb_filtered <- countdata_mb[, -c(2:8)]
rownames(countdata_mb_filtered) <- countdata_mb_filtered$Geneid
countdata_mb_filtered <- countdata_mb_filtered[, -1]

# Read metadata
metadata <- read.csv(metadata_path)
metadata_filtered <- metadata[metadata$Species != "H.m", ]

# Check if sample names match between count data and metadata
metadata_samples <- metadata_filtered$Sample
countdata_mb_samples <- colnames(countdata_mb_filtered)
if (!all(metadata_samples %in% countdata_mb_samples)) stop("Not all sample names in metadata are present in countdata.")
if (!all(countdata_mb_samples %in% metadata_samples)) stop("Not all sample names in countdata are present in metadata.")

# Ensure count data columns are in the same order as metadata samples
countdata_filtered <- countdata_mb_filtered[, match(metadata_samples, countdata_mb_samples)]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(countdata_filtered),
  colData = metadata_filtered,
  design = ~ Tissue + Condition
)

# Run DESeq2
dds <- DESeq(dds)

# Perform variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Load GTF file to extract gene annotations
txdb <- makeTxDbFromGFF(gtf_file_path, format = "gtf")

# Extract results for the specified contrast (modify contrast as needed)
res <- results(dds, contrast = c('Tissue', 'Fertilizer_cavity_warty', 'Fertilizer_cavity_smooth'))

# Filter genes for heatmap (top 20 most expressed genes)
top_genes <- order(rowMeans(assay(vsd)), decreasing = TRUE)[1:20]

# Create a heatmap for selected genes using the base `heatmap()` function
heatmap(assay(vsd)[top_genes, ],
        Colv = NA,            # Don't cluster columns
        scale = "row",        # Scale gene expression across rows
        margins = c(5, 10),   # Set margins for the heatmap
        labRow = rownames(assay(vsd)[top_genes, ]),  # Gene names
        labCol = colnames(assay(vsd)),               # Sample names
        col = colorRampPalette(c("blue", "white", "red"))(100))  # Color palette
############################# pheatmap

# Assuming 'dds' is already created
library("pheatmap")

# Perform variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Select top 20 genes by normalized mean counts
select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:20]

# Extract normalized transformed counts for selected genes
mat <- assay(vsd)[select, ]

# Set the column names to the tissue names from the metadata
colnames(mat) <- colData(dds)$Tissue

pheatmap(mat,
         show_rownames = TRUE,    # Show gene names on the y-axis
         show_colnames = TRUE,    # Show the tissue names on the x-axis
         cluster_cols = FALSE,    # Disable clustering for columns if needed
         scale = "row",           # Scale the data by rows
         main = "Heatmap of Top 20 Genes by Tissue")
###############################################################

# Assuming 'dds' and 'vsd' are already created

# Extract tissue and plant labels from the metadata
plant_labels <- metadata_filtered$Plant       # Assuming Plant column contains P1, P2, P3, etc.
tissue_labels <- metadata_filtered$Tissue     # Assuming Tissue column contains tissue types

# Combine tissue and plant labels for clarity (if needed later)
combined_labels <- paste(tissue_labels, "(", plant_labels, ")", sep = "")

# Define the desired order of tissue grouping (if applicable)
desired_tissues <- c("No_ants_Tuber", "No_ants_Smooth", "No_ants_Warty", 
                     "Control_cavity_tuber", "Control_cavity_smooth", "Control_cavity_warty")  # Adjust order as needed

# Get the order of the columns based on the tissue grouping
desired_order <- unlist(lapply(desired_tissues, function(Tissue) {
  which(tissue_labels == Tissue)
}))

# Reorder the count data and combined labels according to the desired tissue grouping
countdata_ordered <- countdata_filtered[, desired_order]
combined_labels_ordered <- combined_labels[desired_order]
tissue_labels_ordered <- tissue_labels[desired_order]  # Only the tissue names for x-axis

# Select the top 20 most expressed genes (based on normalized counts)
select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:20]

# Extract normalized transformed counts for selected genes, in the desired column order
mat <- assay(vsd)[select, desired_order]

# Set the column names to only the tissue labels for the x-axis
colnames(mat) <- tissue_labels_ordered  # Only tissue labels

# Plot the heatmap with tissue names on the x-axis
pheatmap(mat,
         show_rownames = TRUE,    # Show gene names on the y-axis
         show_colnames = TRUE,    # Show the tissue names on the x-axis
         cluster_cols = FALSE,    # Disable clustering for columns if needed
         scale = "row",           # Scale the data by rows
         main = "Heatmap of Top 20 Genes by Tissue")

















