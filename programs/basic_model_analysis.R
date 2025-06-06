 
library(DESeq2)
library(edgeR)
useArticleSplit = T
countsAllDeseq = read.csv2("../data/rawdata/countsNormalized.tsv", sep = "\t")
colnames(countsAllDeseq) = gsub(".", "-", colnames(countsAllDeseq), fixed = T)
useArticleSplit = T
for(i in 1:ncol(countsAllDeseq))
  countsAllDeseq[,i] = as.numeric(countsAllDeseq[,i])
sampleInfoAll = read.csv2("../data/rawdata/sampleInfo.tsv", sep = "\t")
load("../data/dgeGenesEnsembl75.RData")
# rownames(epilepsyRaw) = epilepsyRaw$X1 

normalizeDESeq2NoReport = function ( rawCounts )
{
  library(DESeq2)
  
  d = rawCounts
  sampleIds = colnames(rawCounts)
  condition <- factor(rep("A",dim(d)[2]))
  
  dds <- DESeqDataSetFromMatrix(countData = d, DataFrame(condition), design = ~ 1)
  dds <- DESeq(dds)
  vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
  
  logs_norm_arr_ord = assay(vsd)
  return(logs_norm_arr_ord) 
}

filter.for.platelet.transcriptome <- function(dge, 
                                              minimum.read.counts = 30,
                                              verbose = TRUE){
  # Filters the DGE-object containing the raw data for low-abundant RNAs.
  # 
  # Args:
  #   dge: DGEList outputted by the prepare.dge.object-function, contains raw count table,
  #       sample info and gene info.
  #   minimum.read.counts: Numeric-value containing the minimum number of gene counts to 
  #       be detected in at least 90% of the samples.
  #   verbose:  Whether or not to show function output.
  #
  # Returns:
  #   DGEList with filtered count table, sample info table and gene info table.
  
  if (missing(dge)){
    stop("Provide DGElist object")
  }
  stopifnot(class(dge) == "DGEList")
  
  if (!is.numeric(minimum.read.counts)){
    stop("Provide numeric value for minimum.read.counts")
  }
  
  # filter for low-abundant transcripts, 
  # i.e. those transcripts with less than minimum.read.counts reads in more 
  # than 90% of all s amples
  bad.genes <- names(which(apply(dge$counts, 1, function(x){sum(x < minimum.read.counts)}) > 0.9 * ncol(dge)))
  dgeFiltered <- dge[which(!rownames(dge) %in% bad.genes), ]
  
  if (verbose == TRUE){
    print(paste("Number of transcripts detected: ", nrow(dgeFiltered$counts), sep = ""))
  }
  
  # return DGEList-object
  return(dgeFiltered)
}





  


my.wilcox.p.value <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
} 

my.t.test.p.value <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
} 

library(matrixStats)
dataFiltered = countsAllDeseq
# not IDH wildtype, excluded as non GBM according to new classification
id604 = grep("TR3712-GBM-VUMC", sampleInfoAll$id, fixed = T)  
if(length(id604) > 0)
{ 
  dataFiltered = dataFiltered[, -id604]
  sampleInfoAll = sampleInfoAll[-id604, ]
} 

genes = genes[which(!is.na(match(genes$hgnc_symbol, rownames(dataFiltered)))),] 
rownames(genes) = genes$hgnc_symbol
library(readxl)  
getSummaryTable = function(dataFiltered, sampleInfo, groupA, groupB,
                           all = F)
{ 
  hcId = which(sampleInfo$ori_patientgroup == groupA)
  ocId =  which(sampleInfo$ori_patientgroup == groupB)
  logsBenign = dataFiltered[, hcId]
  logsMalignant = dataFiltered[,ocId]
   
  meanBenign <- rowMeans(logsBenign)
  meanMalignant <- rowMeans(logsMalignant)
  
  medianBenign <- rowMedians(as.matrix(logsBenign))
  medianMalignant <-  rowMedians(as.matrix(logsMalignant))
  folds_mean_Benign_Malignant = 2^(meanBenign - meanMalignant)  
  ttest_Benign_Malignant_pvalue = sapply(1:nrow(logsBenign), function(i)
    my.t.test.p.value(logsBenign[i, ], logsMalignant[i, ]))
  ttest_Benign_Malignant_pvalue = sapply(1:nrow(logsBenign), function(i)
    my.wilcox.p.value(unlist(logsBenign[i, ]), unlist(logsMalignant[i, ]))) 
  genesLocal = genes[rownames(dataFiltered),] 
  fdrt_Benign_Malignant_qvalue = p.adjust(ttest_Benign_Malignant_pvalue, method = "fdr", n = length(ttest_Benign_Malignant_pvalue))
  summary_table_Benign_Malignant = cbind.data.frame( genesLocal$ensembl_gene_id
                                                     ,rownames(logsBenign),
                                                    meanBenign, meanMalignant,
                                                    medianBenign, medianMalignant, folds_mean_Benign_Malignant,
                                                    ttest_Benign_Malignant_pvalue, fdrt_Benign_Malignant_qvalue, 
                                                    genesLocal$description 
  )
  
  groupNameA = groupA
  groupNameB = groupB
  tableCols_Benign_Malignant=  cbind( "ENSG","Gene",paste("Mean ", groupNameA, sep = ""), paste("Mean ", groupNameB, sep = ""),
                                      paste("Median ", groupNameA, sep = ""), paste("Median ", groupNameB, sep = ""),
                                      paste("Mean fold change ", groupNameA, "-", groupNameB, sep = ""),
                                      "T-test p-value ",
                                      "FDR q-value ", 
                                      "Description" 
  ) 
  colnames(summary_table_Benign_Malignant) = tableCols_Benign_Malignant
  summary_table_Benign_Malignant = summary_table_Benign_Malignant[order(
    summary_table_Benign_Malignant[,9]),]
  if(all)
  {
    cNames = colnames(summary_table_Benign_Malignant)  
    return(summary_table_Benign_Malignant)
  }
  summaryTable = summary_table_Benign_Malignant[which(summary_table_Benign_Malignant$`FDR q-value ` < 0.1),]
  summaryTable = summaryTable[c(which(summaryTable[,7] > 1.3),
                                which(summaryTable[,7] < 1/1.3)), ]
  
  lowExprCutOff = 3
  highHc = which(summaryTable$`Median HC` > lowExprCutOff)
  highOc = which(summaryTable$`Median OC` > lowExprCutOff)
  highExpr = unique(c(highHc, highOc))
  summaryTable = summaryTable[highExpr,] 
  return(summaryTable)
  
}
 
sampleInfoAll$full_ori_patientgroup = sampleInfoAll$ori_patientgroup
gbmSamples = sampleInfoAll[which(sampleInfoAll$ori_patientgroup == "GBM"),] 
library(caret) 
hcId = which(sampleInfoAll$ori_patientgroup == "HC" & sampleInfoAll$isolationlocation != "NKI"  & sampleInfoAll$group == "nonMalignant")

gbmId = which(sampleInfoAll$ori_patientgroup == "GBM" & sampleInfoAll$isolationlocation != "NKI")
toInclude = c(hcId, gbmId)
sampleInfoAll = sampleInfoAll[toInclude,]
countsAllDeseq = countsAllDeseq[, toInclude]
hcId = which(sampleInfoAll$ori_patientgroup == "HC" & sampleInfoAll$isolationlocation != "NKI"  & sampleInfoAll$group == "nonMalignant")
gbmId = which(sampleInfoAll$ori_patientgroup == "GBM" & sampleInfoAll$isolationlocation != "NKI")

set.seed(12345)

folds = createFolds(sampleInfoAll$ori_patientgroup[toInclude], k = 10)
trainId = toInclude[unlist(folds[1:6])]
testId =toInclude[ unlist(folds[7:10])]
if(useArticleSplit)
{
  trainId = read.table("../data/rawdata/train_ids_final.tsv")
  trainId = trainId$V1
  trainId = trainId[which(!is.na(match(trainId, colnames(countsAllDeseq))))]
  trainId = (match(trainId, colnames(countsAllDeseq)))
  testId = read.table("../data/rawdata/testId_ids_final.tsv")
  testId = testId$V1
  testId = testId[which(!is.na(match(testId, colnames(countsAllDeseq))))]
  testId = (match(testId, colnames(countsAllDeseq)))
}
gbmHc_comparisonTrimmed  = getSummaryTable(countsAllDeseq[,trainId], sampleInfoAll[trainId,], "GBM", "HC", F) 
# candidates
candidates = rownames(gbmHc_comparisonTrimmed) 
countsModel = countsAllDeseq[candidates, ] 
countsModel_deseq = countsModel[, toInclude]
countsModel_deseq_assignment = data.frame(c(colnames(countsAllDeseq)[trainId], colnames(countsAllDeseq)[testId]))
colnames(countsModel_deseq_assignment) = "Sample"
countsModel_deseq_assignment$Train = 0
countsModel_deseq_assignment$Train[1:length(trainId)] = 1
countsModel_deseq_assignment$Test = 1
countsModel_deseq_assignment$Test[1:length(trainId)] = 0
deseq_output = gbmHc_comparisonTrimmed   

trainX = countsModel[, trainId]

testX = countsModel[, testId] 
trainY = factor(sampleInfoAll$ori_patientgroup[trainId])
testY = factor(sampleInfoAll$ori_patientgroup[testId]) 
library(party) 
trainX = t(trainX ) 
trainX = as.data.frame(trainX) 
trainX$group = as.factor(trainY)

testX = t(testX ) 
testX = as.data.frame(testX) 
testX$group = as.factor(testY)
library(dplyr) 
library(glmnet)
library(caret) 
realRefsTrain = trainX$group
realRefsTest = testX$group
control <- trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 5,
                        search = "random",
                        verboseIter = TRUE)
elasticModel <- train(group ~ .,
                       data = trainX,
                       method = "glmnet",
                       #preProcess = c("center", "scale"),
                       tuneLength = 25,
                       trControl = control)
elasticModel

predictionsTrain = predict(elasticModel, trainX)
probabilitiesTrain = predict(elasticModel, trainX, type = "prob")
# save(elasticModel, file = "regressionModel.RData")

table(predictionsTrain, realRefsTrain)
confusionMatrixTrain = confusionMatrix(predictionsTrain, realRefsTrain)
confusionMatrixTrain

predictionsTest = predict(elasticModel, testX)
probabilitiesTest = predict(elasticModel, testX, type = "prob")

boxplot(probabilitiesTest[,1] ~ realRefsTest)

thresholdedPredictionsTest = rep("HC", length(predictionsTest))
thresholdedPredictionsTest[which(probabilitiesTest[,1] > 0.25)] = "GBM"

thresholdedPredictionsTrain = rep("HC", length(realRefsTrain))
thresholdedPredictionsTrain[which(probabilitiesTrain[,1] > 0.25)] = "GBM"

confusionMatrixTrain = confusionMatrix(as.factor(thresholdedPredictionsTrain), realRefsTrain) 

table(thresholdedPredictionsTest, realRefsTest)
confusionMatrix(as.factor(thresholdedPredictionsTest), as.factor(realRefsTest), positive = "GBM")

table(predictionsTest, realRefsTest)
confusionMatrixTest = confusionMatrix(as.factor(thresholdedPredictionsTest), reference = as.factor(realRefsTest), positive = "GBM")
confusionMatrixTest

# save.image(paste0("afterElasticNet", Sys.Date(), ".RData")) 
weightThreshold = 0.05 
toKeepElNet =  stats::predict(elasticModel$finalModel, type = "coefficients", s = elasticModel$bestTune$lambda) 
toKeepElNetDF = as.matrix(toKeepElNet)
toKeepElNetDF = as.data.frame(toKeepElNetDF)
toKeepElNetDF = toKeepElNetDF[-1,, drop = F]
fullModel = rownames(toKeepElNetDF)[which(abs(toKeepElNetDF$s1) > 0.0)]
fullModel = gsub("`", "", fullModel)
# write.table(fullModel, file = "fullModel_final.tsv", row.names = F, col.names = F, sep = "\t", quote = F) 
maxWeights = rowMaxs(as.matrix(abs(toKeepElNetDF)))
toKeepElNetDF = toKeepElNetDF[which(maxWeights > weightThreshold),, drop =F] 
allSelectedEnsg = rownames(toKeepElNetDF) 
allSelectedEnsg = gsub("`", "", allSelectedEnsg) 
library(cvAUC)
cvAUC( probabilitiesTest[,1], realRefsTest, label.ordering = c("HC", "GBM"))
AUC =  cvAUC( probabilitiesTest[,1], realRefsTest, label.ordering = c("HC", "GBM")) 
cmTest <- confusionMatrixTest
plt <- as.data.frame(cmTest$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

pdf(paste0('../Report/',
  "Confusion_matrix_test",  ".pdf" ) )
confPlot = ggplot(plt, aes(Prediction,Reference, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#009194") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels=c("HC","GBM")) +
  scale_y_discrete(labels=c("GBM","HC"))


print(confPlot)
dev.off()
library(ggforce)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)

library(org.Hs.eg.db)


library(ggforce)
library(dplyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(openxlsx)
  
deGenes = allSelectedEnsg 
geneUniverse = rownames(dataFiltered) 
library(openxlsx)


c3.tf <- read.gmt("../data/c2.cp.reactome.v2022.1.Hs.symbols.gmt") 
presentFromReact = which(is.element(c3.tf$gene,geneUniverse))# rownames(dataFiltered)))
c3.tf = c3.tf[presentFromReact,]
 
ans.react = enricher(allSelectedEnsg, TERM2GENE=c3.tf , pAdjustMethod = "fdr", qvalueCutoff = 0.5)
dfReact = ans.react@result[1:10,]
noGenes = strsplit(dfReact$GeneRatio, "/")
noGenesPlus = as.numeric(unlist( lapply(noGenes, `[[`, 1)))
noGenesAll = as.numeric(unlist( lapply(noGenes, `[[`, 2)))

noBg = strsplit(dfReact$BgRatio, "/")
noBgPlus = as.numeric(unlist( lapply(noBg, `[[`, 1)))
noBgAll = as.numeric(unlist( lapply(noBg, `[[`, 2)))
dfReact$noGenesPlus = noGenesPlus
dfReact$noGenesAll = noGenesAll
dfReact$noBgPlus = noBgPlus
dfReact$noBgAll = noBgAll
dfReact$Pct_selected = round(100*dfReact$noGenesPlus/dfReact$noGenesAll,2)
dfReact$Pct_background = round(100*dfReact$noBgPlus/dfReact$noBgAll,2)
dfReact$LogP = -1*(log10(dfReact$pvalue))

dfReactChartP = cbind(dfReact$Description, dfReact$LogP, rep("p", nrow(dfReact)), paste0("p=", round(dfReact$pvalue,6)))
dfReactChartPct_selected = cbind(dfReact$Description, dfReact$Pct_selected, rep("Percentage selected", nrow(dfReact)), paste0( dfReact$Pct_selected, "%") )
dfReactChartPct_background = cbind(dfReact$Description, dfReact$Pct_background, rep("Percentage background", nrow(dfReact)), paste0( dfReact$Pct_background, "%") )
dfReactChart = rbind(dfReactChartP, dfReactChartPct_selected,  dfReactChartPct_background)
dfReactChart = as.data.frame(dfReactChart)
colnames(dfReactChart) = c("Pathway", "Value", "Type", "Text")
dfReactChart$Value = as.numeric(dfReactChart$Value)
dfReactChart$Pathway = gsub("REACTOME_", "", dfReactChart$Pathway)
save(dfReactChart, file = "dfReactChart_v3.RData")
pdf("../Report/reactome_top_10.pdf", width = 20, height = 8)
ggplot(data=dfReactChart, aes(x=Pathway,y=Value,fill=factor(Type))) +
  geom_bar(position="dodge",stat="identity") +
  coord_flip() +
  ylab("- Log10 p") +
  ggtitle("Reactome analysis") + 
  geom_text(aes(label = Text),
            position = position_dodge(width= 1), size = 3) + 
  theme_bw()
dev.off()
toSaveX = ans.react@result[which(ans.react@result$pvalue < 0.05),]
toSaveX = toSaveX[, -c(6,7)]
openxlsx::write.xlsx(toSaveX, file = "Reactome_final_selection_vs_all_considered.xlsx") 

deGenes = allSelectedEnsg 


ensgs = genes$ensembl_gene_id[match(deGenes, genes$hgnc_symbol)]
 
ensgUniverse = genes$ensembl_gene_id[match(rownames(dataFiltered), genes$hgnc_symbol)]

ensgUniverse = ensgUniverse[which(!is.na(ensgUniverse))]
deGenes =   unlist(mget(ensgs[which(!is.na(ensgs))], envir=org.Hs.egENSEMBL2EG,
                        ifnotfound = NA))
symbolUniverse = unlist(mget(ensgUniverse[which(!is.na(ensgUniverse))], envir=org.Hs.egENSEMBL2EG,
                             ifnotfound = NA))
symbolUniverse = symbolUniverse[which(!is.na(symbolUniverse))]
set.seed(123)

ans.go <- enrichGO(gene = deGenes, ont = "All",
                   OrgDb ="org.Hs.eg.db",
                   universe = symbolUniverse,
                   readable=TRUE,
                   pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "fdr")

View(ans.go@result) 
ans.go3 <- enrichGO(gene = deGenes, ont = "BP",
                    OrgDb ="org.Hs.eg.db",
                    universe = symbolUniverse,
                    readable=TRUE,
                    pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "fdr")

dfGo = ans.go3@result[1:10,]
noGenes = strsplit(dfGo$GeneRatio, "/")
noGenesPlus = as.numeric(unlist( lapply(noGenes, `[[`, 1)))
noGenesAll = as.numeric(unlist( lapply(noGenes, `[[`, 2)))

noBg = strsplit(dfGo$BgRatio, "/")
noBgPlus = as.numeric(unlist( lapply(noBg, `[[`, 1)))
noBgAll = as.numeric(unlist( lapply(noBg, `[[`, 2)))
dfGo$noGenesPlus = noGenesPlus
dfGo$noGenesAll = noGenesAll
dfGo$noBgPlus = noBgPlus
dfGo$noBgAll = noBgAll
dfGo$Pct_selected = round(100*dfGo$noGenesPlus/dfGo$noGenesAll,2)
dfGo$Pct_background = round(100*dfGo$noBgPlus/dfGo$noBgAll,2)
dfGo$LogP = -1*(log10(dfGo$pvalue))

dfGoChartP = cbind(dfGo$Description, dfGo$LogP, rep("p", nrow(dfGo)), paste0("p=", round(dfGo$pvalue,6)))
dfGoChartPct_selected = cbind(dfGo$Description, dfGo$Pct_selected, rep("Percentage selected", nrow(dfGo)), paste0( dfGo$Pct_selected, "%") )
dfGoChartPct_background = cbind(dfGo$Description, dfGo$Pct_background, rep("Percentage background", nrow(dfGo)), paste0( dfGo$Pct_background, "%") )
dfGoChart = rbind(dfGoChartP, dfGoChartPct_selected,  dfGoChartPct_background)
dfGoChart = as.data.frame(dfGoChart)
colnames(dfGoChart) = c("Pathway", "Value", "Type", "Text")
dfGoChart$Value = as.numeric(dfGoChart$Value)
save(dfGoChart, file = "dfGoChart.RData")
pdf("../Report/gene_ontology_top_10.pdf", width = 20, height = 8)
ggplot(data=dfGoChart, aes(x=Pathway,y=Value,fill=factor(Type))) +
  geom_bar(position="dodge",stat="identity") +
  coord_flip() +
  ylab("- Log10 p") +
  ggtitle("Gene Ontology analysis") + 
  geom_text(aes(label = Text),
            position = position_dodge(width= 1), size = 3) + 
  theme_bw()
dev.off()

selectedVars = rownames(toKeepElNetDF)[which(abs(toKeepElNetDF$s1) > weightThreshold)]
selectedVars = gsub("`", "", selectedVars)  # Clean variable names

# Subset trainX and testX to include only selected variables (plus the group column)
trainX_subset = trainX[, c(selectedVars, "group")]
testX_subset = testX[, c(selectedVars, "group")]

# Define training control
control <- trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 5,
                        search = "random",
                        verboseIter = TRUE)

# Train new elastic net model with selected variables
elasticModel_subset <- train(group ~ .,
                             data = trainX_subset,
                             method = "glmnet",
                             tuneLength = 25,
                             trControl = control)
elasticModel_subset

# Predictions and probabilities on training set
predictionsTrain = predict(elasticModel_subset, trainX_subset)
probabilitiesTrain = predict(elasticModel_subset, trainX_subset, type = "prob")

# Confusion matrix for training set
table(predictionsTrain, realRefsTrain)
confusionMatrixTrain = confusionMatrix(predictionsTrain, realRefsTrain)
confusionMatrixTrain

# Predictions and probabilities on test set
predictionsTest = predict(elasticModel_subset, testX_subset)
probabilitiesTest = predict(elasticModel_subset, testX_subset, type = "prob")

# Boxplot of probabilities
boxplot(probabilitiesTest[,1] ~ realRefsTest)

# Thresholded predictions for test set
thresholdedPredictionsTest = rep("HC", length(predictionsTest))
thresholdedPredictionsTest[which(probabilitiesTest[,1] > 0.25)] = "GBM"

# Thresholded predictions for training set
thresholdedPredictionsTrain = rep("HC", length(realRefsTrain))
thresholdedPredictionsTrain[which(probabilitiesTrain[,1] > 0.25)] = "GBM"

# Confusion matrix for thresholded predictions (training)
confusionMatrixTrain = confusionMatrix(as.factor(thresholdedPredictionsTrain), realRefsTrain) 
confusionMatrixTrain

# Confusion matrix for thresholded predictions (test)
table(thresholdedPredictionsTest, realRefsTest)
confusionMatrixTest = confusionMatrix(as.factor(thresholdedPredictionsTest), as.factor(realRefsTest), positive = "GBM")
confusionMatrixTest

# Confusion matrix for raw predictions (test)
table(predictionsTest, realRefsTest)
confusionMatrixTest = confusionMatrix(as.factor(predictionsTest), as.factor(realRefsTest), positive = "GBM")
confusionMatrixTest
