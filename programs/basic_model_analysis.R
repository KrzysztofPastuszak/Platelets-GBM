# =============================================================================
# Platelet RNA-based classification of Glioblastoma (GBM) vs Healthy Controls
# =============================================================================
#
# This script performs the TDEA (Threshold-Based Differential Expression Analysis)
# with Elastic Net regularization pipeline described in Giczewska, Pastuszak et al.
#
# The analysis proceeds in three stages:
#   1. Data loading, filtering, and sample exclusion
#   2. Differential expression analysis on the training set to identify candidate genes
#   3. Elastic Net classification (with secondary feature selection) followed by
#      Gene Ontology and Reactome pathway enrichment of the final gene panel
#
# Input data:
#   - DESeq2-normalized platelet RNA-seq counts (4,412 genes x 403 samples)
#   - Sample metadata including patient group (GBM / HC), isolation site, etc.
#   - Pre-defined train/test split IDs (60/40)
#   - Ensembl 75 gene annotations
#   - Reactome GMT gene set file (v2022.1)
#
# Output:
#   - Confusion matrix plot (PDF) for the test set
#   - Reactome and Gene Ontology enrichment bar plots (PDF)
#   - Enrichment results table (XLSX)
#   - Elastic Net model coefficients (saved as RData)
# =============================================================================

library(DESeq2)
library(edgeR)

# Flag to use the fixed train/test split from the manuscript (reproducibility).
# When TRUE, sample IDs are read from file rather than generated via createFolds.
useArticleSplit = T

# --- 1. Load normalized counts and sample metadata ---

countsAllDeseq = read.csv2("../data/rawdata/countsNormalized.tsv", sep = "\t")
colnames(countsAllDeseq) = gsub(".", "-", colnames(countsAllDeseq), fixed = T)

# Ensure all columns are numeric (read.csv2 may import as character)
for(i in 1:ncol(countsAllDeseq))
  countsAllDeseq[,i] = as.numeric(countsAllDeseq[,i])

sampleInfoAll = read.csv2("../data/rawdata/sampleInfo.tsv", sep = "\t")

# Gene annotation table (Ensembl 75): maps HGNC symbols to Ensembl IDs
load("../data/dgeGenesEnsembl75.RData")

# =============================================================================
# Helper functions
# =============================================================================

# DESeq2 variance-stabilizing normalization (standalone, without design contrast).
# Used when re-normalizing raw counts outside the main pipeline.
normalizeDESeq2NoReport = function(rawCounts) {
  library(DESeq2)
  
  d = rawCounts
  sampleIds = colnames(rawCounts)
  # Single-group dummy design (normalization only, no DE testing)
  condition <- factor(rep("A", dim(d)[2]))
  
  dds <- DESeqDataSetFromMatrix(countData = d, DataFrame(condition), design = ~ 1)
  dds <- DESeq(dds)
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  
  logs_norm_arr_ord = assay(vsd)
  return(logs_norm_arr_ord)
}

# Filter out low-abundance transcripts from a DGEList object.
# Removes genes with fewer than `minimum.read.counts` reads in > 90% of samples.
filter.for.platelet.transcriptome <- function(dge,
                                              minimum.read.counts = 30,
                                              verbose = TRUE) {
  if (missing(dge)) {
    stop("Provide DGElist object")
  }
  stopifnot(class(dge) == "DGEList")
  
  if (!is.numeric(minimum.read.counts)) {
    stop("Provide numeric value for minimum.read.counts")
  }
  
  # Identify genes that fall below the read-count threshold in >90% of samples
  bad.genes <- names(which(
    apply(dge$counts, 1, function(x) { sum(x < minimum.read.counts) }) > 0.9 * ncol(dge)
  ))
  dgeFiltered <- dge[which(!rownames(dge) %in% bad.genes), ]
  
  if (verbose == TRUE) {
    print(paste("Number of transcripts detected: ", nrow(dgeFiltered$counts), sep = ""))
  }
  
  return(dgeFiltered)
}

# Safe wrappers for Wilcoxon and t-test that return NA on failure
# (e.g., when all values in a group are identical)
my.wilcox.p.value <- function(...) {
  obj <- try(wilcox.test(...), silent = TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

my.t.test.p.value <- function(...) {
  obj <- try(t.test(...), silent = TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

# =============================================================================
# 2. Sample exclusion and cohort definition
# =============================================================================

library(matrixStats)
dataFiltered = countsAllDeseq

# Exclude one IDH-mutant sample reclassified as non-GBM under updated WHO criteria.
# TR3712-GBM-VUMC was originally labeled GBM but is IDH-wildtype negative.
id604 = grep("TR3712-GBM-VUMC", sampleInfoAll$id, fixed = T)
if(length(id604) > 0) {
  dataFiltered = dataFiltered[, -id604]
  sampleInfoAll = sampleInfoAll[-id604, ]
}

# Align gene annotation table to the genes present in filtered counts
genes = genes[which(!is.na(match(genes$hgnc_symbol, rownames(dataFiltered)))),]
rownames(genes) = genes$hgnc_symbol

# =============================================================================
# Differential expression summary table (Wilcoxon rank-sum test + FDR)
# =============================================================================
#
# For two patient groups (groupA vs groupB), computes per-gene:
#   - Mean and median expression in each group
#   - Fold change (on log2 scale, back-transformed to linear)
#   - Wilcoxon rank-sum p-value with FDR correction
#
# When `all = FALSE` (default), applies TDEA thresholds:
#   - FDR q-value < 0.1
#   - Fold change > 1.3 or < 1/1.3
#   - At least one group median > 3 (low-expression filter)

getSummaryTable = function(dataFiltered, sampleInfo, groupA, groupB,
                           all = F) {
  hcId = which(sampleInfo$ori_patientgroup == groupA)
  ocId = which(sampleInfo$ori_patientgroup == groupB)
  logsBenign = dataFiltered[, hcId]
  logsMalignant = dataFiltered[, ocId]
  
  # Group-level summary statistics
  meanBenign <- rowMeans(logsBenign)
  meanMalignant <- rowMeans(logsMalignant)
  
  medianBenign <- rowMedians(as.matrix(logsBenign))
  medianMalignant <- rowMedians(as.matrix(logsMalignant))
  
  # Fold change computed from mean log2-scale values (DESeq2-normalized)
  folds_mean_Benign_Malignant = 2^(meanBenign - meanMalignant)
  
  # Per-gene Wilcoxon rank-sum test between the two groups
  wilcox_Benign_Malignant_pvalue = sapply(1:nrow(logsBenign), function(i)
    my.wilcox.p.value(unlist(logsBenign[i, ]), unlist(logsMalignant[i, ])))
  
  genesLocal = genes[rownames(dataFiltered),]
  
  # FDR correction (Benjamini-Hochberg)
  fdr_Benign_Malignant_qvalue = p.adjust(
    wilcox_Benign_Malignant_pvalue, method = "fdr",
    n = length(wilcox_Benign_Malignant_pvalue)
  )
  
  # Assemble results table
  summary_table_Benign_Malignant = cbind.data.frame(
    genesLocal$ensembl_gene_id,
    rownames(logsBenign),
    meanBenign, meanMalignant,
    medianBenign, medianMalignant,
    folds_mean_Benign_Malignant,
    wilcox_Benign_Malignant_pvalue,
    fdr_Benign_Malignant_qvalue,
    genesLocal$description
  )
  
  groupNameA = groupA
  groupNameB = groupB
  tableCols_Benign_Malignant = cbind(
    "ENSG", "Gene",
    paste("Mean ", groupNameA, sep = ""), paste("Mean ", groupNameB, sep = ""),
    paste("Median ", groupNameA, sep = ""), paste("Median ", groupNameB, sep = ""),
    paste("Mean fold change ", groupNameA, "-", groupNameB, sep = ""),
    "Wilcoxon p-value ",
    "FDR q-value ",
    "Description"
  )
  colnames(summary_table_Benign_Malignant) = tableCols_Benign_Malignant
  
  # Sort by FDR q-value (ascending)
  summary_table_Benign_Malignant = summary_table_Benign_Malignant[order(
    summary_table_Benign_Malignant[, 9]), ]
  
  # Return all genes without filtering if requested
  if(all) {
    return(summary_table_Benign_Malignant)
  }
  
  # Apply TDEA thresholds: FDR < 0.1, fold change > 1.3, minimum expression > 3
  summaryTable = summary_table_Benign_Malignant[
    which(summary_table_Benign_Malignant$`FDR q-value ` < 0.1), ]
  summaryTable = summaryTable[c(
    which(summaryTable[, 7] > 1.3),
    which(summaryTable[, 7] < 1/1.3)
  ), ]
  
  lowExprCutOff = 3
  # NOTE: "Median OC" was the original column name when this function was used
  # with ovarian cancer data. When called with groupA="GBM", groupB="HC", the
  # actual column is "Median GBM". Because $`Median OC` returns NULL here,
  # the filter only requires HC median > 3. This is preserved for result
  # reproducibility with the published manuscript.
  highHc = which(summaryTable$`Median HC` > lowExprCutOff)
  highOc = which(summaryTable$`Median OC` > lowExprCutOff)
  highExpr = unique(c(highHc, highOc))
  summaryTable = summaryTable[highExpr, ]
  return(summaryTable)
}

# =============================================================================
# 3. Define the study cohort (GBM + HC, excluding NKI controls)
# =============================================================================

sampleInfoAll$full_ori_patientgroup = sampleInfoAll$ori_patientgroup
gbmSamples = sampleInfoAll[which(sampleInfoAll$ori_patientgroup == "GBM"),]

library(caret)

# Select healthy controls from non-NKI sites only (NKI excluded per study protocol)
hcId = which(sampleInfoAll$ori_patientgroup == "HC" &
               sampleInfoAll$isolationlocation != "NKI" &
               sampleInfoAll$group == "nonMalignant")

# Select GBM patients from non-NKI sites
gbmId = which(sampleInfoAll$ori_patientgroup == "GBM" &
                sampleInfoAll$isolationlocation != "NKI")

toInclude = c(hcId, gbmId)
sampleInfoAll = sampleInfoAll[toInclude, ]
countsAllDeseq = countsAllDeseq[, toInclude]

# Refresh indices after subsetting
hcId = which(sampleInfoAll$ori_patientgroup == "HC" &
               sampleInfoAll$isolationlocation != "NKI" &
               sampleInfoAll$group == "nonMalignant")
gbmId = which(sampleInfoAll$ori_patientgroup == "GBM" &
                sampleInfoAll$isolationlocation != "NKI")

# =============================================================================
# 4. Train/test split (60/40)
# =============================================================================

set.seed(12345)

# Default: 10-fold stratified split (folds 1-6 = train, 7-10 = test)
folds = createFolds(sampleInfoAll$ori_patientgroup[toInclude], k = 10)
trainId = toInclude[unlist(folds[1:6])]
testId = toInclude[unlist(folds[7:10])]

# Override with the fixed split used in the manuscript for reproducibility
if(useArticleSplit) {
  trainId = read.table("../data/rawdata/train_ids_final.tsv")
  trainId = trainId$V1
  trainId = trainId[which(!is.na(match(trainId, colnames(countsAllDeseq))))]
  trainId = (match(trainId, colnames(countsAllDeseq)))
  
  testId = read.table("../data/rawdata/testId_ids_final.tsv")
  testId = testId$V1
  testId = testId[which(!is.na(match(testId, colnames(countsAllDeseq))))]
  testId = (match(testId, colnames(countsAllDeseq)))
}

# =============================================================================
# 5. TDEA: identify differentially expressed candidate genes on TRAINING data only
# =============================================================================

# Run differential expression on training set (leak-free: test set not used)
gbmHc_comparisonTrimmed = getSummaryTable(
  countsAllDeseq[, trainId], sampleInfoAll[trainId, ], "GBM", "HC", F
)

# Candidate genes passing TDEA thresholds (FDR < 0.1, FC > 1.3, expression > 3)
candidates = rownames(gbmHc_comparisonTrimmed)

# Subset counts to candidate genes for modeling
countsModel = countsAllDeseq[candidates, ]
countsModel_deseq = countsModel[, toInclude]

# Record the train/test assignment for each sample
countsModel_deseq_assignment = data.frame(
  c(colnames(countsAllDeseq)[trainId], colnames(countsAllDeseq)[testId])
)
colnames(countsModel_deseq_assignment) = "Sample"
countsModel_deseq_assignment$Train = 0
countsModel_deseq_assignment$Train[1:length(trainId)] = 1
countsModel_deseq_assignment$Test = 1
countsModel_deseq_assignment$Test[1:length(trainId)] = 0
deseq_output = gbmHc_comparisonTrimmed

# =============================================================================
# 6. Elastic Net classification (first pass: all candidate genes)
# =============================================================================

# Prepare training and test matrices (samples x genes)
trainX = countsModel[, trainId]
testX = countsModel[, testId]
trainY = factor(sampleInfoAll$ori_patientgroup[trainId])
testY = factor(sampleInfoAll$ori_patientgroup[testId])

library(party)

trainX = t(trainX)
trainX = as.data.frame(trainX)
trainX$group = as.factor(trainY)

testX = t(testX)
testX = as.data.frame(testX)
testX$group = as.factor(testY)

library(dplyr)
library(glmnet)

realRefsTrain = trainX$group
realRefsTest = testX$group

# 5x5 repeated cross-validation with random hyperparameter search
control <- trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 5,
                        search = "random",
                        verboseIter = TRUE)

# Train Elastic Net (glmnet) on all TDEA-selected candidate genes
elasticModel <- train(group ~ .,
                      data = trainX,
                      method = "glmnet",
                      tuneLength = 25,
                      trControl = control)
print(elasticModel)

# --- Evaluate first-pass model ---

predictionsTrain = predict(elasticModel, trainX)
probabilitiesTrain = predict(elasticModel, trainX, type = "prob")

table(predictionsTrain, realRefsTrain)
confusionMatrixTrain = confusionMatrix(predictionsTrain, realRefsTrain)
print(confusionMatrixTrain)

predictionsTest = predict(elasticModel, testX)
probabilitiesTest = predict(elasticModel, testX, type = "prob")

if (interactive()) boxplot(probabilitiesTest[, 1] ~ realRefsTest)

# Apply a lowered probability threshold (0.25) to increase GBM sensitivity
thresholdedPredictionsTest = rep("HC", length(predictionsTest))
thresholdedPredictionsTest[which(probabilitiesTest[, 1] > 0.25)] = "GBM"

thresholdedPredictionsTrain = rep("HC", length(realRefsTrain))
thresholdedPredictionsTrain[which(probabilitiesTrain[, 1] > 0.25)] = "GBM"

confusionMatrixTrain = confusionMatrix(as.factor(thresholdedPredictionsTrain), realRefsTrain)

table(thresholdedPredictionsTest, realRefsTest)
print(confusionMatrix(as.factor(thresholdedPredictionsTest), as.factor(realRefsTest), positive = "GBM"))

table(predictionsTest, realRefsTest)
confusionMatrixTest = confusionMatrix(
  as.factor(thresholdedPredictionsTest),
  reference = as.factor(realRefsTest),
  positive = "GBM"
)
print(confusionMatrixTest)

# =============================================================================
# 7. Secondary feature selection: keep genes with |coefficient| > threshold
# =============================================================================

weightThreshold = 0.05

# Extract Elastic Net coefficients at the optimal lambda
toKeepElNet = stats::predict(
  elasticModel$finalModel, type = "coefficients",
  s = elasticModel$bestTune$lambda
)
toKeepElNetDF = as.data.frame(as.matrix(toKeepElNet))
toKeepElNetDF = toKeepElNetDF[-1, , drop = F]  # Remove intercept row
toKeepElNetDF[, 1] = as.numeric(toKeepElNetDF[, 1])  # Ensure plain numeric

# Full model: all genes with non-zero coefficients
fullModel = rownames(toKeepElNetDF)[which(abs(toKeepElNetDF[, 1]) > 0.0)]
fullModel = gsub("`", "", fullModel)

# Reduced panel: genes with |weight| above threshold
maxWeights = abs(toKeepElNetDF[, 1])
toKeepElNetDF = toKeepElNetDF[which(maxWeights > weightThreshold), , drop = F]
allSelectedEnsg = rownames(toKeepElNetDF)
allSelectedEnsg = gsub("`", "", allSelectedEnsg)

# =============================================================================
# 8. Test set AUC and confusion matrix visualization
# =============================================================================

library(cvAUC)
AUC = cvAUC(probabilitiesTest[, 1], realRefsTest, label.ordering = c("HC", "GBM"))

# Confusion matrix heatmap for the test set
cmTest <- confusionMatrixTest
plt <- as.data.frame(cmTest$table)
plt$Prediction <- factor(plt$Prediction, levels = rev(levels(plt$Prediction)))

pdf(paste0('../Report/', "Confusion_matrix_test", ".pdf"))
confPlot = ggplot(plt, aes(Prediction, Reference, fill = Freq)) +
  geom_tile() + geom_text(aes(label = Freq)) +
  scale_fill_gradient(low = "white", high = "#009194") +
  labs(x = "Reference", y = "Prediction") +
  scale_x_discrete(labels = c("HC", "GBM")) +
  scale_y_discrete(labels = c("GBM", "HC"))

print(confPlot)
dev.off()

# =============================================================================
# 9. Reactome pathway enrichment (ORA) on the final Elastic Net gene panel
# =============================================================================

library(org.Hs.eg.db)
library(clusterProfiler)
library(openxlsx)

deGenes = allSelectedEnsg         # Final panel genes (HGNC symbols)
geneUniverse = rownames(dataFiltered)  # Background: all expressed genes

# Load Reactome gene sets (MSigDB format)
c3.tf <- read.gmt("../data/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
presentFromReact = which(is.element(c3.tf$gene, geneUniverse))
c3.tf = c3.tf[presentFromReact, ]

# Run overrepresentation analysis against Reactome pathways
ans.react = enricher(allSelectedEnsg, TERM2GENE = c3.tf,
                     pAdjustMethod = "fdr", qvalueCutoff = 0.5)

# --- Prepare Reactome bar chart data ---

dfReact = ans.react@result[1:10, ]

# Parse the GeneRatio and BgRatio strings into numeric values
noGenes = strsplit(dfReact$GeneRatio, "/")
noGenesPlus = as.numeric(unlist(lapply(noGenes, `[[`, 1)))
noGenesAll = as.numeric(unlist(lapply(noGenes, `[[`, 2)))

noBg = strsplit(dfReact$BgRatio, "/")
noBgPlus = as.numeric(unlist(lapply(noBg, `[[`, 1)))
noBgAll = as.numeric(unlist(lapply(noBg, `[[`, 2)))

dfReact$noGenesPlus = noGenesPlus
dfReact$noGenesAll = noGenesAll
dfReact$noBgPlus = noBgPlus
dfReact$noBgAll = noBgAll
dfReact$Pct_selected = round(100 * dfReact$noGenesPlus / dfReact$noGenesAll, 2)
dfReact$Pct_background = round(100 * dfReact$noBgPlus / dfReact$noBgAll, 2)
dfReact$LogP = -1 * (log10(dfReact$pvalue))

# Reshape to long format for grouped bar chart
dfReactChartP = cbind(dfReact$Description, dfReact$LogP,
                      rep("p", nrow(dfReact)),
                      paste0("p=", round(dfReact$pvalue, 6)))
dfReactChartPct_selected = cbind(dfReact$Description, dfReact$Pct_selected,
                                 rep("Percentage selected", nrow(dfReact)),
                                 paste0(dfReact$Pct_selected, "%"))
dfReactChartPct_background = cbind(dfReact$Description, dfReact$Pct_background,
                                   rep("Percentage background", nrow(dfReact)),
                                   paste0(dfReact$Pct_background, "%"))

dfReactChart = rbind(dfReactChartP, dfReactChartPct_selected, dfReactChartPct_background)
dfReactChart = as.data.frame(dfReactChart)
colnames(dfReactChart) = c("Pathway", "Value", "Type", "Text")
dfReactChart$Value = as.numeric(dfReactChart$Value)

# Strip the REACTOME_ prefix for cleaner labels
dfReactChart$Pathway = gsub("REACTOME_", "", dfReactChart$Pathway)
save(dfReactChart, file = "dfReactChart_v3.RData")

# Plot top 10 Reactome pathways (grouped bars: -log10 p, % selected, % background)
pdf("../Report/reactome_top_10.pdf", width = 20, height = 8)
ggplot(data = dfReactChart, aes(x = Pathway, y = Value, fill = factor(Type))) +
  geom_bar(position = "dodge", stat = "identity") +
  coord_flip() +
  ylab("- Log10 p") +
  ggtitle("Reactome analysis") +
  geom_text(aes(label = Text),
            position = position_dodge(width = 1), size = 3) +
  theme_bw()
dev.off()

# Save significant Reactome results (p < 0.05) to Excel
toSaveX = ans.react@result[which(ans.react@result$pvalue < 0.05), ]
toSaveX = toSaveX[, -c(6, 7)]
openxlsx::write.xlsx(toSaveX, file = "Reactome_final_selection_vs_all_considered.xlsx")

# =============================================================================
# 10. Gene Ontology enrichment (ORA) on the final Elastic Net gene panel
# =============================================================================

deGenes = allSelectedEnsg

# Map HGNC symbols to Ensembl IDs, then to Entrez IDs (required by enrichGO)
ensgs = genes$ensembl_gene_id[match(deGenes, genes$hgnc_symbol)]

ensgUniverse = genes$ensembl_gene_id[match(rownames(dataFiltered), genes$hgnc_symbol)]
ensgUniverse = ensgUniverse[which(!is.na(ensgUniverse))]

deGenes = unlist(mget(ensgs[which(!is.na(ensgs))], envir = org.Hs.egENSEMBL2EG,
                      ifnotfound = NA))
symbolUniverse = unlist(mget(ensgUniverse[which(!is.na(ensgUniverse))],
                             envir = org.Hs.egENSEMBL2EG, ifnotfound = NA))
symbolUniverse = symbolUniverse[which(!is.na(symbolUniverse))]

set.seed(123)

# GO enrichment across all ontologies (for exploration)
ans.go <- enrichGO(gene = deGenes, ont = "All",
                   OrgDb = "org.Hs.eg.db",
                   universe = symbolUniverse,
                   readable = TRUE,
                   pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "fdr")
# Print top GO results (use View() interactively for full table)
head(ans.go@result, 20)

# GO enrichment restricted to Biological Process (for the manuscript figure)
ans.go3 <- enrichGO(gene = deGenes, ont = "BP",
                    OrgDb = "org.Hs.eg.db",
                    universe = symbolUniverse,
                    readable = TRUE,
                    pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "fdr")

# --- Prepare GO bar chart data (same structure as Reactome) ---

dfGo = ans.go3@result[1:10, ]

noGenes = strsplit(dfGo$GeneRatio, "/")
noGenesPlus = as.numeric(unlist(lapply(noGenes, `[[`, 1)))
noGenesAll = as.numeric(unlist(lapply(noGenes, `[[`, 2)))

noBg = strsplit(dfGo$BgRatio, "/")
noBgPlus = as.numeric(unlist(lapply(noBg, `[[`, 1)))
noBgAll = as.numeric(unlist(lapply(noBg, `[[`, 2)))

dfGo$noGenesPlus = noGenesPlus
dfGo$noGenesAll = noGenesAll
dfGo$noBgPlus = noBgPlus
dfGo$noBgAll = noBgAll
dfGo$Pct_selected = round(100 * dfGo$noGenesPlus / dfGo$noGenesAll, 2)
dfGo$Pct_background = round(100 * dfGo$noBgPlus / dfGo$noBgAll, 2)
dfGo$LogP = -1 * (log10(dfGo$pvalue))

dfGoChartP = cbind(dfGo$Description, dfGo$LogP,
                   rep("p", nrow(dfGo)),
                   paste0("p=", round(dfGo$pvalue, 6)))
dfGoChartPct_selected = cbind(dfGo$Description, dfGo$Pct_selected,
                              rep("Percentage selected", nrow(dfGo)),
                              paste0(dfGo$Pct_selected, "%"))
dfGoChartPct_background = cbind(dfGo$Description, dfGo$Pct_background,
                                rep("Percentage background", nrow(dfGo)),
                                paste0(dfGo$Pct_background, "%"))

dfGoChart = rbind(dfGoChartP, dfGoChartPct_selected, dfGoChartPct_background)
dfGoChart = as.data.frame(dfGoChart)
colnames(dfGoChart) = c("Pathway", "Value", "Type", "Text")
dfGoChart$Value = as.numeric(dfGoChart$Value)
save(dfGoChart, file = "dfGoChart.RData")

# Plot top 10 GO Biological Process terms
pdf("../Report/gene_ontology_top_10.pdf", width = 20, height = 8)
ggplot(data = dfGoChart, aes(x = Pathway, y = Value, fill = factor(Type))) +
  geom_bar(position = "dodge", stat = "identity") +
  coord_flip() +
  ylab("- Log10 p") +
  ggtitle("Gene Ontology analysis") +
  geom_text(aes(label = Text),
            position = position_dodge(width = 1), size = 3) +
  theme_bw()
dev.off()

# =============================================================================
# 11. Refit Elastic Net on reduced gene panel (|coefficient| > threshold)
# =============================================================================

selectedVars = rownames(toKeepElNetDF)
selectedVars = gsub("`", "", selectedVars)

# Subset training and test data to the reduced gene panel
trainX_subset = trainX[, c(selectedVars, "group")]
testX_subset = testX[, c(selectedVars, "group")]

# Same CV setup as before
control <- trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 5,
                        search = "random",
                        verboseIter = TRUE)

# Refit Elastic Net on the reduced feature set
elasticModel_subset <- train(group ~ .,
                             data = trainX_subset,
                             method = "glmnet",
                             tuneLength = 25,
                             trControl = control)
print(elasticModel_subset)

# --- Evaluate reduced model on training set ---
predictionsTrain = predict(elasticModel_subset, trainX_subset)
probabilitiesTrain = predict(elasticModel_subset, trainX_subset, type = "prob")

table(predictionsTrain, realRefsTrain)
confusionMatrixTrain = confusionMatrix(predictionsTrain, realRefsTrain)
print(confusionMatrixTrain)

# --- Evaluate reduced model on test set ---
predictionsTest = predict(elasticModel_subset, testX_subset)
probabilitiesTest = predict(elasticModel_subset, testX_subset, type = "prob")

if (interactive()) boxplot(probabilitiesTest[, 1] ~ realRefsTest)

# Thresholded predictions (probability > 0.25 classified as GBM)
thresholdedPredictionsTest = rep("HC", length(predictionsTest))
thresholdedPredictionsTest[which(probabilitiesTest[, 1] > 0.25)] = "GBM"

thresholdedPredictionsTrain = rep("HC", length(realRefsTrain))
thresholdedPredictionsTrain[which(probabilitiesTrain[, 1] > 0.25)] = "GBM"

# Confusion matrices for thresholded predictions
confusionMatrixTrain = confusionMatrix(as.factor(thresholdedPredictionsTrain), realRefsTrain)
print(confusionMatrixTrain)

table(thresholdedPredictionsTest, realRefsTest)
confusionMatrixTest = confusionMatrix(
  as.factor(thresholdedPredictionsTest), as.factor(realRefsTest), positive = "GBM"
)
print(confusionMatrixTest)

# Confusion matrix for default-threshold predictions
table(predictionsTest, realRefsTest)
confusionMatrixTest = confusionMatrix(
  as.factor(predictionsTest), as.factor(realRefsTest), positive = "GBM"
)
print(confusionMatrixTest)

# =============================================================================
# 12. Three-way GO Biological Process enrichment across pipeline stages
# =============================================================================
#
# Compares GO BP enrichment results across the three nested gene sets produced
# by the pipeline:
#   1. candidates         -- TDEA-selected genes (FDR < 0.1, FC > 1.3, expr > 3)
#   2. fullModel          -- all genes with non-zero Elastic Net coefficients
#   3. allSelectedEnsg    -- final panel (|coefficient| > 0.05)
#
# A Venn diagram shows overlap of significant GO terms between sets,
# produced for both FDR-adjusted (p.adjust < 0.05) and unadjusted (p < 0.05)
# significance thresholds.

library(VennDiagram)
library(grid)

# --- Convert each gene set from HGNC symbols to Entrez IDs ---
# (enrichGO requires Entrez IDs; universe already computed in section 10)

symbolsToEntrez = function(symbols) {
  ensIds = genes$ensembl_gene_id[match(symbols, genes$hgnc_symbol)]
  ensIds = ensIds[which(!is.na(ensIds))]
  entrezIds = unlist(mget(ensIds, envir = org.Hs.egENSEMBL2EG, ifnotfound = NA))
  entrezIds[which(!is.na(entrezIds))]
}

entrezCandidates  = symbolsToEntrez(candidates)       # TDEA-selected
entrezFullModel   = symbolsToEntrez(fullModel)         # non-zero coefficients
entrezFinalPanel  = symbolsToEntrez(allSelectedEnsg)   # |coef| > threshold

# --- Run GO BP enrichment for each gene set ---

set.seed(123)

goBpCandidates = enrichGO(
  gene = entrezCandidates, ont = "BP", OrgDb = "org.Hs.eg.db",
  universe = symbolUniverse, readable = TRUE,
  pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "fdr"
)

goBpFullModel = enrichGO(
  gene = entrezFullModel, ont = "BP", OrgDb = "org.Hs.eg.db",
  universe = symbolUniverse, readable = TRUE,
  pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "fdr"
)

goBpFinalPanel = enrichGO(
  gene = entrezFinalPanel, ont = "BP", OrgDb = "org.Hs.eg.db",
  universe = symbolUniverse, readable = TRUE,
  pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "fdr"
)

# Save significant results (unadjusted p < 0.05) to Excel
write.xlsx(
  goBpCandidates@result[which(goBpCandidates@result$pvalue < 0.05), -c(6, 7)],
  file = "GO_BP_candidates_vs_all.xlsx"
)
write.xlsx(
  goBpFullModel@result[which(goBpFullModel@result$pvalue < 0.05), -c(6, 7)],
  file = "GO_BP_fullModel_vs_all.xlsx"
)
write.xlsx(
  goBpFinalPanel@result[which(goBpFinalPanel@result$pvalue < 0.05), -c(6, 7)],
  file = "GO_BP_finalSelection_vs_all.xlsx"
)

# --- Extract significant pathway names at both thresholds ---

# FDR-adjusted (p.adjust < 0.05)
sigAdjCandidates  = goBpCandidates@result$Description[which(goBpCandidates@result$p.adjust < 0.05)]
sigAdjFullModel   = goBpFullModel@result$Description[which(goBpFullModel@result$p.adjust < 0.05)]
sigAdjFinalPanel  = goBpFinalPanel@result$Description[which(goBpFinalPanel@result$p.adjust < 0.05)]

# Unadjusted (p < 0.05)
sigRawCandidates  = goBpCandidates@result$Description[which(goBpCandidates@result$pvalue < 0.05)]
sigRawFullModel   = goBpFullModel@result$Description[which(goBpFullModel@result$pvalue < 0.05)]
sigRawFinalPanel  = goBpFinalPanel@result$Description[which(goBpFinalPanel@result$pvalue < 0.05)]

cat("Significant GO BP terms (FDR-adjusted p < 0.05):\n")
cat("  TDEA candidates:", length(sigAdjCandidates), "\n")
cat("  Non-zero in model:", length(sigAdjFullModel), "\n")
cat("  Final panel:", length(sigAdjFinalPanel), "\n")

cat("\nSignificant GO BP terms (unadjusted p < 0.05):\n")
cat("  TDEA candidates:", length(sigRawCandidates), "\n")
cat("  Non-zero in model:", length(sigRawFullModel), "\n")
cat("  Final panel:", length(sigRawFinalPanel), "\n")

# --- Venn diagram helper ---

vennColors = c("#4477AA", "#EE6677", "#228833")
vennSetNames = c("TDEA candidates", "Non-zero in model", "Final panel")

drawGoVenn = function(pathList, titleText, outFile) {
  pdf(outFile, width = 7, height = 7)
  vennObj = venn.diagram(
    x = setNames(pathList, vennSetNames),
    filename = NULL,
    fill = vennColors,
    alpha = 0.45,
    col = vennColors,
    lwd = 1.5,
    fontfamily = "sans",
    cat.fontfamily = "sans",
    cat.fontface = "bold",
    cat.cex = 1.1,
    cex = 1.3,
    margin = 0.1,
    main = titleText,
    main.fontfamily = "sans",
    main.cex = 1.3,
    main.fontface = "bold"
  )
  grid.draw(vennObj)
  dev.off()
}

# --- Venn diagram: FDR-adjusted significance ---

drawGoVenn(
  pathList  = list(sigAdjCandidates, sigAdjFullModel, sigAdjFinalPanel),
  titleText = "GO BP pathways (FDR-adjusted p < 0.05)",
  outFile   = "../Report/GO_BP_venn_adjusted.pdf"
)

# --- Venn diagram: unadjusted significance ---

drawGoVenn(
  pathList  = list(sigRawCandidates, sigRawFullModel, sigRawFinalPanel),
  titleText = "GO BP pathways (unadjusted p < 0.05)",
  outFile   = "../Report/GO_BP_venn_unadjusted.pdf"
)

# --- Save overlapping pathway names ---

sharedAdj = Reduce(intersect, list(sigAdjCandidates, sigAdjFullModel, sigAdjFinalPanel))
sharedRaw = Reduce(intersect, list(sigRawCandidates, sigRawFullModel, sigRawFinalPanel))

cat("\nPathways significant in all three sets (FDR-adjusted):", length(sharedAdj), "\n")
cat("Pathways significant in all three sets (unadjusted):", length(sharedRaw), "\n")

if (length(sharedAdj) > 0) {
  write.xlsx(data.frame(Pathway = sharedAdj),
             file = "GO_BP_shared_all_three_adjusted.xlsx")
}
if (length(sharedRaw) > 0) {
  write.xlsx(data.frame(Pathway = sharedRaw),
             file = "GO_BP_shared_all_three_unadjusted.xlsx")
}

# =============================================================================
# 13. Build ranked gene list for GSEA (GBM vs HC, training set)
# =============================================================================
#
# Gene ranking metric: sign(log2 FC) * -log10(Wilcoxon p-value)
# This preserves directionality (GBM-up vs GBM-down) while using the
# significance of the test to weight the magnitude of the ranking.
#
# All expressed genes are used here -- NOT just TDEA candidates --
# because GSEA operates on a full ranked list, not a pre-filtered set.
# =============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)
library(VennDiagram)
library(grid)

# Full (unfiltered) DE results on training set
gbmHc_allGenes = getSummaryTable(
  countsAllDeseq[, trainId],
  sampleInfoAll[trainId, ],
  "GBM", "HC",
  all = TRUE
)

fc_col    = which(colnames(gbmHc_allGenes) == "Mean fold change GBM-HC")
pval_col  = which(colnames(gbmHc_allGenes) == "Wilcoxon p-value ")

fc_vals   = as.numeric(gbmHc_allGenes[, fc_col])
pval_vals = as.numeric(gbmHc_allGenes[, pval_col])

# Guard against p = 0 (would produce Inf)
pval_vals[pval_vals == 0] <- .Machine$double.eps

rank_metric = sign(log2(fc_vals)) * (-log10(pval_vals))

# Map HGNC symbols -> Ensembl IDs (gseGO uses ENSEMBL keyType)
gene_symbols  = rownames(gbmHc_allGenes)
ensembl_ids   = genes$ensembl_gene_id[match(gene_symbols, genes$hgnc_symbol)]

original_gene_list_gbm = stats::setNames(rank_metric, ensembl_ids)

# Drop unmapped and duplicated entries
original_gene_list_gbm = original_gene_list_gbm[!is.na(names(original_gene_list_gbm))]
original_gene_list_gbm = original_gene_list_gbm[!duplicated(names(original_gene_list_gbm))]

gene_list_gbm = na.omit(original_gene_list_gbm)
gene_list_gbm = sort(gene_list_gbm, decreasing = TRUE)

# =============================================================================
# 14. GSEA across multiple seeds (GO Biological Process)
# =============================================================================
#
# GSEA results can vary across permutations, so 10 seeds are run and saved
# individually. The seed=1 run (pvalueCutoff = 1) is kept in memory for
# downstream visualisation and the ORA/GSEA Venn diagram.
# =============================================================================
library(enrichplot)
for (seed in 1:10) {
  set.seed(seed)
  gse_gbm_iter <- gseGO(
    geneList      = gene_list_gbm,
    ont           = "ALL",
    keyType       = "ENSEMBL",
    nPerm         = 10000,
    minGSSize     = 3,
    maxGSSize     = 800,
    pvalueCutoff  = 0.05,
    verbose       = TRUE,
    OrgDb         = org.Hs.eg.db,
    pAdjustMethod = "fdr"
  )
  write.xlsx(gse_gbm_iter,
             file = paste0("gene_set_enrichment_all_hc_vs_gbm_ref_hc3_", seed, ".xlsx"))
  save(gse_gbm_iter,
       file = paste0("gene_set_enrichment_all_hc_vs_gbm_ref_hc3_", seed, ".RData"))
}

# Re-run seed=1 with pvalueCutoff = 1 and ont = "BP" to match the ORA
# used in the Venn. Full result set needed for top-100 selection below.
set.seed(1)
gse_gbm <- gseGO(
  geneList      = gene_list_gbm,
  ont           = "BP",
  keyType       = "ENSEMBL",
  nPerm         = 10000,
  minGSSize     = 3,
  maxGSSize     = 800,
  pvalueCutoff  = 1,          # no cutoff -- needed to get >= 100 results
  verbose       = FALSE,
  OrgDb         = org.Hs.eg.db,
  pAdjustMethod = "fdr"
)

# =============================================================================
# 15. GSEA visualisation
# =============================================================================

library(DOSE)
library(ggplot2)

pdf("../Report/gsea_dotplot_BP.pdf", width = 10, height = 8)
print(
  dotplot(gse_gbm, showCategory = 10, split = ".sign") +
    facet_grid(. ~ .sign) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#F0F0F0"),
      axis.text.y      = element_text(size = 9)
    )
)
dev.off()

# emapplot requires pairwise term similarity
pdf("../Report/gsea_emapplot_BP.pdf", width = 10, height = 10)
print(emapplot(pairwise_termsim(gse_gbm), showCategory = 15))
dev.off()

# =============================================================================
# 16. Venn diagram: ORA (top 100, no FDR) vs GSEA (top 100, no FDR)
# =============================================================================
#
# Both sets are GO Biological Process terms ranked by raw p-value.
# ORA is the enrichment on the final Elastic Net gene panel (ans.go3,
# computed in section 10). GSEA is the seed=1 BP run above.
#
# "Top 100" is defined as the 100 terms with the smallest raw p-value in
# each analysis, without any FDR or adjusted-p filter.
# =============================================================================

n_top = 100

# ORA: ans.go3 was computed with pvalueCutoff = 1 so all terms are present
ora_results_sorted = ans.go3@result[order(ans.go3@result$pvalue), ]
ora_top = head(ora_results_sorted$Description, n_top)

# GSEA: gse_gbm also has pvalueCutoff = 1
gsea_results_sorted = gse_gbm@result[order(gse_gbm@result$pvalue), ]
gsea_top = head(gsea_results_sorted$Description, n_top)

n_overlap = length(intersect(ora_top, gsea_top))
cat(sprintf(
  "\nGO BP Venn summary (top %d by raw p-value):\n  ORA terms:  %d\n  GSEA terms: %d\n  Overlap:    %d\n",
  n_top, length(ora_top), length(gsea_top), n_overlap
))

# Palette consistent with the three-way Venn in section 12
venn_pal = c("#4477AA", "#EE6677")

pdf("../Report/GO_BP_venn_ORA_vs_GSEA.pdf", width = 6, height = 6)
venn_obj = venn.diagram(
  x          = list(
    `ORA`  = ora_top,
    `GSEA` = gsea_top
  ),
  filename        = NULL,
  fill            = venn_pal,
  alpha           = 0.45,
  col             = venn_pal,
  lwd             = 1.8,
  fontfamily      = "sans",
  cat.fontfamily  = "sans",
  cat.fontface    = "bold",
  cat.cex         = 1.15,
  cex             = 1.4,
  margin          = 0.12,
  main            = sprintf("GO BP top %d pathways: ORA vs GSEA", n_top),
  main.fontfamily = "sans",
  main.cex        = 1.15,
  main.fontface   = "bold"
)
grid.newpage()
grid.draw(venn_obj)
dev.off()

# Save the overlapping pathway names
shared_ora_gsea = intersect(ora_top, gsea_top)
if (length(shared_ora_gsea) > 0) {
  write.xlsx(
    data.frame(
      Pathway         = shared_ora_gsea,
      ORA_pvalue      = ora_results_sorted$pvalue[
        match(shared_ora_gsea, ora_results_sorted$Description)],
      ORA_padj        = ora_results_sorted$p.adjust[
        match(shared_ora_gsea, ora_results_sorted$Description)],
      GSEA_pvalue     = gsea_results_sorted$pvalue[
        match(shared_ora_gsea, gsea_results_sorted$Description)],
      GSEA_padj       = gsea_results_sorted$p.adjust[
        match(shared_ora_gsea, gsea_results_sorted$Description)],
      GSEA_NES        = gsea_results_sorted$NES[
        match(shared_ora_gsea, gsea_results_sorted$Description)]
    ),
    file = "GO_BP_shared_ORA_GSEA_top100.xlsx"
  )
}

