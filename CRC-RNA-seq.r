# ==============================================================================
# RNA-SEQ & MACHINE LEARNING BIOMARKER DISCOVERY PIPELINE
# Project: Colorectal Cancer (CRC) Biomarker Discovery
# Goal: End-to-end analysis from raw counts to predictive modeling 
#       using DESeq2, WGCNA, and Random Forest.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. SETUP AND LIBRARIES
# ------------------------------------------------------------------------------
message("\n[STEP 1] Initializing workspace and loading libraries...")

cran_packages <- c(
  "tidyverse", "pheatmap", "RColorBrewer", "ggplot2", "ggrepel", 
  "pROC", "randomForest", "caret", "here", "janitor", "gridExtra"
)

bioc_packages <- c(
  "DESeq2", "EnhancedVolcano", "ComplexHeatmap", "WGCNA", "mygene"
)

# Install missing CRAN packages
new_cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(new_cran)) install.packages(new_cran)

# Install missing Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
new_bioc <- setdiff(bioc_packages, rownames(installed.packages()))
if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE)

# Load core libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(EnhancedVolcano)
  library(ComplexHeatmap)
  library(WGCNA)
  library(mygene)
  library(caret)
  library(randomForest)
  library(pROC)
})

# Create output directory structure for GitHub
dirs <- c("results/qc", "results/differential", "results/figures", 
          "results/tables", "results/heatmaps", "results/wgcna", "results/ml")
invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

# ------------------------------------------------------------------------------
# 2. DATA LOADING & ANNOTATION
# ------------------------------------------------------------------------------
message("\n[STEP 2] Loading raw data and annotating via mygene...")

# Load counts and metadata (Replace with your actual file names)
counts_raw <- read.table("GSE50760_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
metadata <- read.csv("SraRunTable.csv", stringsAsFactors = FALSE)

# Fetch Gene Symbols using Entrez IDs
gene_ids <- as.character(counts_raw$GeneID)
gene_info <- queryMany(gene_ids, scopes="entrezgene", fields=c("symbol","name"), species="human")
df_genes <- as.data.frame(gene_info)
df_genes <- df_genes[!duplicated(df_genes$query), ] # Remove duplicates

# Merge annotation with counts
counts_raw$GeneID <- as.character(counts_raw$GeneID)
df_genes$query <- as.character(df_genes$query)

final_counts <- counts_raw %>%
  left_join(df_genes %>% select(query, symbol, name), by = c("GeneID" = "query")) %>%
  rename(query = GeneID)

# Create a global annotation map
gene_annotation <- final_counts %>% select(query, symbol, name) %>% rename(GeneID = query)
rownames(gene_annotation) <- gene_annotation$GeneID

# ------------------------------------------------------------------------------
# 3. MATRIX FORMATTING
# ------------------------------------------------------------------------------
message("\n[STEP 3] Formatting matrices and synchronizing metadata...")

count_matrix <- final_counts %>% select(-query, -symbol, -name)
rownames(count_matrix) <- final_counts$query

# Synchronize metadata with count matrix columns
metadata$GSM <- metadata$GEO_Accession..exp.
common_samples <- intersect(colnames(count_matrix), metadata$GSM)
count_matrix <- count_matrix[, common_samples, drop = FALSE]

metadata <- metadata %>%
  filter(GSM %in% common_samples) %>%
  arrange(match(GSM, common_samples)) %>%
  mutate(
    condition_original = case_when(
      grepl("normal", tissue, ignore.case = TRUE) ~ "Normal",
      grepl("primary", tissue, ignore.case = TRUE) ~ "Tumor",
      grepl("metast", tissue, ignore.case = TRUE) ~ "Metastasis"
    ),
    condition = factor(condition_original, levels = c("Normal", "Tumor", "Metastasis"))
  )
rownames(metadata) <- metadata$GSM

stopifnot(all(colnames(count_matrix) == rownames(metadata)))

# ------------------------------------------------------------------------------
# 4. NORMALIZATION & PCA (QUALITY CONTROL)
# ------------------------------------------------------------------------------
message("\n[STEP 4] DESeq2 Normalization and PCA...")

# Filter low count genes
keep <- rowSums(count_matrix) >= 10
counts_filt <- count_matrix[keep, , drop = FALSE]

# Initialize DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = counts_filt, colData = metadata, design = ~ condition)
dds <- estimateSizeFactors(dds)

# Variance Stabilizing Transformation for PCA and WGCNA
vst_data <- vst(dds, blind = TRUE)

# --- PCA PLOT ---
pca_df <- plotPCA(vst_data, intgroup = "condition_original", returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition_original)) +
  geom_point(size = 4, alpha = 0.85) +
  stat_ellipse(level = 0.95, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("Normal" = "#4575B4", "Tumor" = "#D6604D", "Metastasis" = "#66A61E")) +
  labs(title = "PCA - Sample Clustering by Condition", 
       x = paste0("PC1: ", percentVar[1], "% variance"), 
       y = paste0("PC2: ", percentVar[2], "% variance"))

ggsave("results/figures/PCA_plot.pdf", p_pca, width = 8, height = 6)

# ------------------------------------------------------------------------------
# 5. DIFFERENTIAL EXPRESSION & VOLCANO PLOT
# ------------------------------------------------------------------------------
message("\n[STEP 5] Differential expression analysis (Tumor vs Normal)...")

dds <- DESeq(dds)
res_TvN <- results(dds, contrast = c("condition", "Tumor", "Normal"))

results_df <- as.data.frame(res_TvN) %>% rownames_to_column("GeneID")
results_annotated <- results_df %>% left_join(gene_annotation, by = "GeneID")

# --- VOLCANO PLOT ---
volcano_plot <- EnhancedVolcano(
  results_annotated,
  lab = ifelse(is.na(results_annotated$symbol) | results_annotated$symbol == "", 
               results_annotated$GeneID, results_annotated$symbol),
  x = "log2FoldChange", y = "padj",
  title = "Tumor vs Normal Expression",
  subtitle = "Differential Expression Analysis",
  pCutoff = 0.01, FCcutoff = 1,
  pointSize = 2.5, labSize = 4,
  col = c("grey30", "forestgreen", "royalblue", "red2"),
  drawConnectors = TRUE, widthConnectors = 0.5
)

ggsave("results/figures/Volcano_plot.pdf", volcano_plot, width = 10, height = 8)

# --- DEG HEATMAP ---
top_degs <- results_annotated %>%
  filter(padj < 0.01, abs(log2FoldChange) > 2) %>%
  arrange(padj) %>% head(40)

mat_deg <- assay(vst_data)[top_degs$GeneID, ]
rownames(mat_deg) <- ifelse(is.na(top_degs$symbol) | top_degs$symbol == "", top_degs$GeneID, top_degs$symbol)
mat_deg_z <- t(scale(t(mat_deg)))

anno_col <- data.frame(Condition = metadata$condition_original)
rownames(anno_col) <- colnames(mat_deg_z)

png("results/heatmaps/DEG_Top40_Heatmap.png", width = 8, height = 10, units = "in", res = 300)
pheatmap(mat_deg_z, annotation_col = anno_col, show_colnames = FALSE, 
         main = "Top 40 Differentially Expressed Genes", clustering_distance_rows = "correlation")
dev.off()

# ------------------------------------------------------------------------------
# 6. WGCNA (CO-EXPRESSION NETWORK)
# ------------------------------------------------------------------------------
message("\n[STEP 6] Executing WGCNA...")

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Use top 5000 highly variable genes to optimize network construction
datExpr <- t(assay(vst_data))
gene_vars <- apply(datExpr, 2, var)
datExpr <- datExpr[, names(sort(gene_vars, decreasing = TRUE))[1:5000]]

# Build network automatically
net <- blockwiseModules(datExpr, power = 14, TOMType = "unsigned", 
                        minModuleSize = 30, numericLabels = FALSE, verbose = 0)

# Correlate modules with clinical traits (Tumor presence)
traits <- data.frame(Tumor = as.numeric(metadata$condition_original == "Tumor"))
rownames(traits) <- rownames(metadata)
moduleTraitCor <- cor(net$MEs, traits, use = "p")
modules_of_interest <- rownames(moduleTraitCor)[abs(moduleTraitCor[,1]) > 0.5]

# ------------------------------------------------------------------------------
# 7. MULTI-OMICS INTEGRATION & HUB SCORING
# ------------------------------------------------------------------------------
message("\n[STEP 7] Integrating DESeq2 and WGCNA outputs...")

kME <- as.data.frame(cor(datExpr, net$MEs, use = "p")) %>% rownames_to_column("gene")

# Merging network features, differential stats, and annotations
hub_candidates <- kME %>%
  pivot_longer(-gene, names_to = "module", values_to = "kME_val") %>%
  filter(module %in% modules_of_interest) %>%
  inner_join(results_df, by = c("gene" = "GeneID")) %>%
  left_join(gene_annotation, by = c("gene" = "GeneID"))

# Calculate unified biological Hub Score and extract top 20 predictive panel
final_panel <- hub_candidates %>%
  filter(!is.na(padj), padj < 0.01, abs(log2FoldChange) > 1, abs(kME_val) > 0.7) %>%
  mutate(
    hub_score = abs(kME_val) * abs(log2FoldChange) * (-log10(padj + 1e-300)),
    symbol = ifelse(is.na(symbol) | symbol == "", gene, symbol)
  ) %>%
  arrange(desc(hub_score)) %>%
  head(20)

write.csv(final_panel, "results/wgcna/Final_Biomarker_Panel.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 8. MACHINE LEARNING VALIDATION (RANDOM FOREST)
# ------------------------------------------------------------------------------
message("\n[STEP 8] Machine Learning validation (Random Forest)...")

# Prepare ML matrix with top 20 features
ml_genes <- as.character(final_panel$gene)
ml_data <- as.data.frame(t(assay(vst_data)[ml_genes, ]))

# Prevent duplication errors in column names for ML algorithms
colnames(ml_data) <- make.names(make.unique(colnames(ml_data)))
ml_data$Condition <- metadata$condition_original

# Filter for binary classification (Tumor vs Normal)
ml_data <- ml_data %>% 
  filter(Condition %in% c("Tumor", "Normal")) %>%
  mutate(Condition = factor(Condition, levels = c("Normal", "Tumor")))

# Train robust Random Forest with 5-fold Cross Validation
set.seed(123)
train_control <- trainControl(
  method = "cv", number = 5, classProbs = TRUE, 
  summaryFunction = twoClassSummary, savePredictions = "final"
)

rf_model <- train(
  Condition ~ ., data = ml_data, method = "rf", 
  metric = "ROC", trControl = train_control, importance = TRUE
)

message(sprintf("  ✓ RF Model trained successfully. CV AUC: %.3f", max(rf_model$results$ROC)))

# ------------------------------------------------------------------------------
# 9. PREDICTIVE PERFORMANCE & FINAL VISUALIZATIONS
# ------------------------------------------------------------------------------
message("\n[STEP 9] Generating Performance Metrics and Annotated Plots...")

# --- ROC CURVE ---
selected_preds <- rf_model$pred %>% filter(mtry == rf_model$bestTune$mtry)
roc_obj <- roc(selected_preds$obs, selected_preds$Tumor, levels = c("Normal", "Tumor"))

png("results/ml/ROC_Curve_Final.png", width = 6, height = 6, units = "in", res = 300)
plot(roc_obj, col = "#D6604D", lwd = 3, main = paste("ROC Curve - 20 Gene Panel\nAUC:", round(auc(roc_obj), 3)))
abline(a = 0, b = 1, lty = 2, col = "grey")
dev.off()

# --- OOB ERROR STABILITY ---
png("results/ml/RF_Error_Stability.png", width = 7, height = 5, units = "in", res = 300)
plot(rf_model$finalModel, main = "Random Forest Error Convergence (OOB Error)")
legend("topright", colnames(rf_model$finalModel$err.rate), col=1:3, fill=1:3)
dev.off()

# --- ANNOTATED VARIABLE IMPORTANCE PLOT ---
importance_rf <- varImp(rf_model, scale = FALSE)
imp_raw <- as.data.frame(importance_rf$importance)
target_col <- colnames(imp_raw)[1] 

imp_df <- imp_raw
imp_df$importance_score <- imp_df[[target_col]]
imp_df$feature <- rownames(imp_df)

# Clean R technical artifacts from feature names
imp_df$clean_id <- gsub("^X", "", imp_df$feature)
imp_df$clean_id <- gsub("\\.[0-9]+$", "", imp_df$clean_id)

# Robust Dictionary Mapping to actual Gene Symbols
anno_valid <- gene_annotation[!is.na(gene_annotation$GeneID) & gene_annotation$GeneID != "", ]
dict <- setNames(as.character(anno_valid$symbol), as.character(anno_valid$GeneID))
imp_df$display_name <- dict[imp_df$clean_id]

# Fill empty matches with the clean ID
imp_df$display_name <- ifelse(is.na(imp_df$display_name) | imp_df$display_name == "", 
                              imp_df$clean_id, as.character(imp_df$display_name))

# Resolve visual duplicates for ggplot factorization
plot_data <- imp_df %>%
  arrange(desc(importance_score)) %>%
  head(20) %>%
  mutate(display_name = make.unique(as.character(display_name))) %>%
  mutate(display_name = factor(display_name, levels = rev(display_name)))

p_imp <- ggplot(plot_data, aes(x = display_name, y = importance_score)) +
  geom_segment(aes(xend = display_name, yend = 0), color = "grey85") +
  geom_point(size = 4, color = "#D6604D") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Top 20 Predictive Biomarkers",
    subtitle = "Ranked by Random Forest Importance",
    x = "Gene Symbol", y = "Importance Score"
  ) +
  theme(axis.text.y = element_text(face = "italic", size = 11, color = "black"),
        panel.grid.minor = element_blank())

ggsave("results/ml/Variable_Importance_Annotated.png", p_imp, width = 7, height = 6, dpi = 300)

message("\n====================================================")
message("PIPELINE COMPLETED SUCCESSFULLY")
message("All plots and tables have been saved in the 'results/' directory.")
message("Ready for GitHub upload.")
# ==============================================================================