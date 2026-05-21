library(data.table)
library(matrixStats)
library(tidyr)
library(ggplot2)
library(dplyr)
library(QFeatures)
library(limma)
library(missBayes)
library(MSstats)
library(DEqMS)
library(proDA)
library(DreamAI)
library(limpa)
library(msqrob2)
library(Pirat)
library(SummarizedExperiment)
library(Matrix)
library(iq)
source("ND_imputation.R")

# Import data ----
pg<-fread("CQE/report.pg_matrix.tsv")
pr<-fread("CQE/report.pr_matrix.tsv")
raw <- fread("CQE/report.tsv")
t <- table(unique(pr[,c('Protein.Group','Stripped.Sequence')])$Protein.Group)
pg$Peptide.Count <- t[match(pg$Protein.Group,names(t))]
metadata <- fread("CQE/metadata.txt", header = TRUE)
#Filter by protein groups with <2 peptides, remove contaminants, and log2 transform the data
df.prot <- filter(pg, Peptide.Count > 1)
df.prot2 <- df.prot[!grepl("Cont_",df.prot$Protein.Group),]
df.prot.all<- df.prot2[, 5:29]
protein.matrix_all <- log2(as.matrix(df.prot.all))
log2all.df <- as.data.frame(protein.matrix_all)
rownames(log2all.df) <- df.prot2$Protein.Names

# 1. MSstats ----
# convert DIANN output to MSstats raw data
MSstats_data <- DIANNtoMSstatsFormat(raw, metadata, removeFewMeasurements = FALSE, 
                                     quantificationColumn = "FragmentQuantRaw" )
timepoint1 <- Sys.time()
# data pre-processing and summarization
MSstats_results <- dataProcess(MSstats_data, normalization = "equalizeMedians", summaryMethod = "TMP", MBimpute = TRUE)
timepoint2 <- Sys.time()
MSstasts_summarization_time <- timepoint2 - timepoint1
# differential abundance analysis for group comparison
MSstats_results$ProteinLevelData$GROUP <- paste0("G", regmatches(MSstats_results$ProteinLevelData$originalRUN, 
                                                                 regexpr("\\d+_\\d+_\\d+", MSstats_results$ProteinLevelData$originalRUN)))
MSstats_results$FeatureLevelData$GROUP <- paste0("G", regmatches(MSstats_results$FeatureLevelData$originalRUN, 
                                                                 regexpr("\\d+_\\d+_\\d+", MSstats_results$FeatureLevelData$originalRUN)))

levels <- levels(as.factor(MSstats_results$ProteinLevelData$GROUP))
comparison <- matrix(0, nrow = 1, ncol = length(levels))
colnames(comparison) <- levels
rownames(comparison) <- "G70_10_20 - G70_20_10"
comparison[1, "G70_10_20"] <- 1
comparison[1, "G70_20_10"] <- -1

timepoint3 <- Sys.time()
testResultOneComparison <- groupComparison(contrast.matrix= comparison, data= MSstats_results)
timepoint4 <- Sys.time()
MSstats_testing_time <- timepoint4 - timepoint3
MSstats.results <- testResultOneComparison$ComparisonResult

# 2. MSqRob2 ----
raw[["File.Name"]] <- raw[["Run"]]
raw_wide <- raw[,c(2:6, 14:16, 28)] %>% pivot_wider(names_from = Run, values_from = Precursor.Normalised, values_fn = sum)

metadata[["Run"]] <- sort(unique(raw[["File.Name"]]))
colnames(metadata)[1] <- "quantCols"

qf <- readQFeatures(raw_wide, colData = metadata, quantCols = metadata$quantCols, name = "quant")
qf <- logTransform(qf, base = 2, i = "quant", name = "log2quant")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$quantCols
metadata_aligned <- metadata[rownames(colData(qf)), ]
colData(qf)$Condition <- as.factor(metadata_aligned$Condition)
colData(qf)$Bio.Rep <- metadata_aligned$Bio.Rep
colData(qf)$Run <- metadata_aligned$Run
# normalization
qf<- normalize(qf, i = "log2quant", name = "log2norm", method = "center.median")
plotDensities(assay(qf[["log2norm"]]), legend = FALSE)

timepoint5 <- Sys.time()
#Fit ridge regression peptide-level model 
qf <- msqrobAggregate(qf, i = "log2norm", fcol = "Protein.Names", formula = ~Condition + (1|Modified.Sequence), ridge = TRUE)
#Fit Hurdle model (to summarized protein quant) 
qf_agg <- aggregateFeatures(qf, i = "log2norm", fcol = "Protein.Names", name = "log2proteinQuant")
qf_hurdle <- msqrobHurdle(qf_agg, i = "log2proteinQuant", formula = ~Condition)
timepoint6 <- Sys.time()
msqrob2_summarization_time <- timepoint6 - timepoint5


# Check names in the Intensity component
coef <- getCoef(rowData(qf_hurdle[["log2proteinQuant"]])$msqrobHurdleIntensity[[1]])
contrast <- makeContrast(contrasts = "ConditionG70_10_20 - ConditionG70_20_10 = 0", parameterNames = names(coef))

timepoint7 <- Sys.time()
msqrob2_result <- hypothesisTestHurdle(object = qf_hurdle, contrast = contrast, i = "log2proteinQuant")
timepoint8 <- Sys.time()
msqrob2_testing_time <- timepoint8 - timepoint7

msqRob2.results_701020 <- rowData(msqrob2_result[["log2proteinQuant"]])[["hurdle_ConditionG70_10_20 - ConditionG70_20_10"]]

# 3. DEqMS ----
# median normalization
pep.count.table <- data.frame(count = as.vector(pg$Peptide.Count),
                              row.names = pg$Protein.Names)

log2norm.df <- normalizeMedianValues(log2all.df)
class <- as.factor(metadata$Condition)
design <- model.matrix(~0+class) # fitting without intercept
colnames(design) <- sub("^class", "", colnames(design))

timepoint9 <- Sys.time()
fit1 <- lmFit(log2norm.df,design = design)
cont <- makeContrasts(G70_10_20-G70_20_10, levels = levels(class))
fit2 <- contrasts.fit(fit1,contrasts = cont)
fit3 <- eBayes(fit2)

fit3$count <- pep.count.table[rownames(fit3$coefficients),"count"]

#check the values in the vector fit3$count
#if min(fit3$count) return NA or 0, troubleshoot the error first
min(fit3$count)
fit4 <- spectraCounteBayes(fit3)
timepoint10 <- Sys.time()
DEqMS_testing_time <- timepoint10 - timepoint9

DEqMS.results <- outputResult(fit4,coef_col = 1)

# 4. proDA ----
timepoint11 <- Sys.time()
fit <- proDA(log2norm.df, design = design) # fit dropout model
proDA.results_701020 <- test_diff(fit, `G70_10_20` - `G70_20_10`)
timepoint12 <- Sys.time()
proDA_testing_time <- timepoint12 - timepoint11

# 5. limpa ----
y.prec <- readDIANN(file = "CQE/report.tsv",
                    run.column = "File.Name", q.columns = c("Q.Value", "Protein.Q.Value"))

timepoint13 <- Sys.time()
dpcest <- dpcCN(y.prec) # estimate dpc 
y <- dpcQuant(y.prec, dpc = dpcest) # summarization
timepoint14 <- Sys.time()
limpa_summarization_time <- timepoint14 - timepoint13

Group <- as.factor(c(rep("G35_25_40", 3), rep("G35_40_25", 2), rep("G70_10_20", 4), rep("G70_20_10", 3), rep("G70_30_0", 3), rep("G35_25_40",2), rep("G35_40_25", 3), "G70_10_20", rep("G70_20_10", 2), rep("G70_30_0", 2)))
design <- model.matrix(~0+Group)
colnames(design) <- sub("^Group", "", colnames(design))
timepoint15 <- Sys.time()
fit <- dpcDE(y, design, plot = TRUE)
fit <- eBayes(fit)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)
timepoint16 <- Sys.time()
limpa_testing_time <- timepoint16 - timepoint15

limpa_result_701020 <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")

# 7. msImpute+limma ----
# Generate MaxLFQ peptide-level intensities
pep_maxlfq <- diann_maxlfq(raw[raw$Q.Value <= 0.01 & raw$PG.Q.Value <= 0.01,], 
                           group.header="Stripped.Sequence", 
                           id.header = "Precursor.Id", 
                           quantity.header = "Precursor.Normalised")

# Convert to log2 
pep_log2 <- log2(pep_maxlfq)
pep_log2[is.infinite(pep_log2)] <- NA  # Ensure 0s become NAs
# median normalization
pep_log2 <- normalizeMedianValues(pep_log2)
# msImpute
design <- model.matrix( ~0+Condition, metadata )
timepoint17 <- Sys.time()
msImpute.pep <- msImpute(pep_log2, method = "v2-mnar", design = design, relax_min_obs = TRUE)
timepoint18 <- Sys.time()
msImpute_imputation_time <- timepoint18 - timepoint17

# robust summarization
# 1. Create a SummarizedExperiment for the peptide data
msImpute.pep <- msImpute.pep[, rownames(colData(qf))]
se_pep <- SummarizedExperiment(assays = list(log2quant = msImpute.pep),
                               colData = colData(qf)) # use the colData from msqorob2
# 2. Wrap it in a QFeatures object
qf_imputed <- QFeatures(list(pep_imputed = se_pep))
# 3. Create a mapping 
# Get unique mappings between sequences and protein groups
pep_to_prot <- unique(raw[, .(Stripped.Sequence, Protein.Group)])
setDF(pep_to_prot)
rownames(pep_to_prot) <- pep_to_prot$Stripped.Sequence

# 4. Add this mapping to the rowData of new QFeatures object
# This ensures that for every row in the matrix, R knows the 'Protein.Group'
rowData(qf_imputed[["pep_imputed"]]) <- pep_to_prot[rownames(msImpute.pep), ]
# Ensure the metadata row names match the colData row names
# metadata should be the original data frame imported from CSV
rownames(metadata) <- metadata$quantCols 
# Re-align and assign the actual character names
# This ensures Sample 'X' gets Condition 'Y' regardless of order
colData(qf_imputed)$Condition <- metadata[rownames(colData(qf_imputed)), "Condition"]
# Convert to factor 
colData(qf_imputed)$Condition <- as.factor(colData(qf_imputed)$Condition)
# Verify 
colData(qf_imputed)$Condition
# 5. Aggregate to Protein Level
timepoint19 <- Sys.time()
qf_imputed <- aggregateFeatures(qf_imputed, 
                                i = "pep_imputed", 
                                fcol = "Protein.Group", 
                                name = "protein_imputed", 
                                fun = MsCoreUtils::robustSummary)
timepoint20 <- Sys.time()
msImpute_summarization_time <- timepoint20 - timepoint19

# extract the aggregated protein-level data
msImpute.results <- assay(qf_imputed[["protein_imputed"]])
conditions_msImpute <- sub(".*HYB_(\\d+_\\d+_\\d+)_.*", "G\\1", colnames(msImpute.results))
design_msImpute <- model.matrix(~0 + conditions_msImpute)
# statistical inference 
timepoint21 <- Sys.time()
fit <- lmFit(msImpute.results, design_msImpute)
contrast_matrix <- makeContrasts(conditions_msImputeG70_10_20 - conditions_msImputeG70_20_10, levels = design_msImpute)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- treat(fit2, lfc=log2(1.15), trend = TRUE)
timepoint22 <- Sys.time()
msImpute_testing_time <- timepoint22 - timepoint21

msImpute.results_701020 <- topTreat(fit2, coef = "conditions_msImputeG70_10_20 - conditions_msImputeG70_20_10", adjust.method = "BH", number = Inf, confint = TRUE)
# replace row names with protein names
# Match row names of df1 to Protein.Group in df2, then pull the corresponding Protein.Names
rownames(msImpute.results_701020) <- raw$Protein.Names[match(rownames(msImpute.results_701020), raw$Protein.Group)]

# 8. DreamAI+limma ----
timepoint23 <- Sys.time()
DreamAI_imputed <- DreamAI(log2norm.df)
timepoint24 <- Sys.time()
DreamAI_imputation_time <- timepoint24 - timepoint23

DreamAI_imputed <- DreamAI_imputed$Ensemble

# limma inference
timepoint25 <- Sys.time()
fit <- lmFit(DreamAI_imputed, design_limma)
contrast_matrix <- makeContrasts(G70_10_20 - G70_20_10, levels = design_limma)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- treat(fit2, lfc=log2(1.15), trend = TRUE)
timepoint26 <- Sys.time()
DreamAI_testing_time <- timepoint26 - timepoint25

DreamAI.results <- topTreat(fit2, coef = "G70_10_20 - G70_20_10", adjust.method = "BH", number = Inf, confint = TRUE)

# 9. mice-pmm+limma ----
colnames(log2norm.df) <- metadata$rename
timepoint27 <- Sys.time()
mice_imputed <- mice(log2norm.df, method = "pmm", seed = 456)
mice_pooled <- Reduce("+", lapply(1:5, function(i) complete(mice_imputed, i))) / 5
mice_imputed <- mice_pooled
timepoint28 <- Sys.time()
mice_imputation_time <- timepoint28 - timepoint27

# limma inference
timepoint29 <- Sys.time()
fit <- lmFit(mice_imputed, design_limma)
contrast_matrix <- makeContrasts(G70_10_20 - G70_20_10, levels = design_limma)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- treat(fit2, lfc=log2(1.15), trend = TRUE)
timepoint30 <- Sys.time()
mice_testing_time <- timepoint30 - timepoint29

mice.results <- topTreat(fit2, coef = "G70_10_20 - G70_20_10", adjust.method = "BH", number = Inf, confint = TRUE)

# 10. missBayes ----
group <- as.factor(metadata$Condition)
comparison <- makeContrasts("G70_10_20 - G70_20_10", levels = levels(group))
set.seed(746)
timepoint31 <- Sys.time()
missBayes.results_701020 <- BayesMissingModel(values = log2norm.df, groups = group, comparisons = comparison, parallel = TRUE,
                             n.adapt = 500, burn.in = 500, n.iter = 10000, n.chains = 2, mcmcDiag = TRUE)
timepoint32 <- Sys.time()
missBayes_testing_time <- timepoint32 - timepoint31

# 11. ND+limma ----
timepoint33 <- Sys.time()
ND_imputed <- impute_normal(log2norm.df)
timepoint34 <- Sys.time()
ND_imputation_time <- timepoint34 - timepoint33

# limma inference
timepoint35 <- Sys.time()
fit <- lmFit(ND_imputed, design_limma)
contrast_matrix <- makeContrasts(G70_10_20 - G70_20_10, levels = design_limma)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- treat(fit2, lfc=log2(1.15), trend = TRUE)
timepoint36 <- Sys.time()
ND_testing_time <- timepoint36 - timepoint35

ND.results <- topTreat(fit2, coef = "G70_10_20 - G70_20_10", adjust.method = "BH", number = Inf, confint = TRUE)


# 12. limma ----
design_limma <- model.matrix(~ 0 + Condition, data = metadata)
colnames(design_limma) <- gsub("Condition", "", colnames(design_limma))
timepoint37 <- Sys.time()
fit <- lmFit(log2norm.df, design_limma)
contrast_matrix <- makeContrasts(G70_10_20 - G70_20_10, levels = design_limma)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- treat(fit2, lfc=log2(1.15), trend = TRUE)
timepoint38 <- Sys.time()
limma_testing_time <- timepoint38 - timepoint37

limma.results <- topTreat(fit2, coef = "G70_10_20 - G70_20_10", adjust.method = "BH", number = Inf, confint = TRUE)

# Keep proteins present in all methods for fair comparison
rownames(MSstats.results) <- MSstats.results$Protein
rownames(limpa_result_701020) <- limpa_result_701020$Protein.Names
rownames(proDA.results_701020) <- proDA.results_701020$name

common_proteins <- Reduce(intersect, list(rownames(MSstats.results_701020), rownames(msqRob2.results_701020), rownames(DEqMS.results_701020), rownames(proDA.results_701020), rownames(limpa_result_701020), 
                                          rownames(DreamAI.results_701020), rownames(mice.results_701020), rownames(msImpute.results_701020), rownames(ND.results_701020), rownames(missBayes.results_701020[[1]]), rownames(limma.results_701020)))
                          
# 13. RMSE calculation ----
# rename logFC columns before calculation
missBayes.results_701020[[1]] <- missBayes.results_701020[[1]] %>% rename(logFC = Median)
MSstats.results <- MSstats.results %>% rename(logFC = log2FC)
proDA.results_701020 <- proDA.results_701020 %>% rename(logFC = diff)

method_list <- list(
  MSstats = MSstats.results_701020[common_proteins, ],
  MSqRob2 = msqRob2.results_701020[common_proteins, ],
  DEqMS = DEqMS.results_701020[common_proteins, ],
  proDA = proDA.results_701020[common_proteins, ],
  limpa = limpa_result_701020[common_proteins, ],
  DreamAI_limma = DreamAI.results_701020[common_proteins, ],
  mice_limma = mice.results_701020[common_proteins, ],
  msImpute_limma = msImpute.results_701020[common_proteins, ],
  ND_limma = ND.results_701020[common_proteins, ],
  missBayes = missBayes.results_701020[[1]][common_proteins, ],
  limma = limma.results_701020[common_proteins, ]
)

RMSE_data <- data.frame(method = names(method_list), RMSE = NA)

# Loop through each named element in the list
for (i in seq_along(method_list)) {
  method_name <- names(method_list)[i]  
  method_df <- method_list[[i]]        
  
  # Add the true value column based on protein names
  method_df$trueVal <- ifelse(grepl("HUMAN", rownames(method_df)), 0,
                              ifelse(grepl("ECOLI", rownames(method_df)), 1,
                                     ifelse(grepl("YEAST", rownames(method_df)), -1, NA)))
  method_df$logFC[!is.finite(method_df$logFC)] <- NA
  
  # Calculate RMSE for this method
  rmse <- sqrt(mean((method_df$trueVal - method_df$logFC)^2, na.rm = TRUE))
  
  # Assign the result to the correct row in RMSE_data
  RMSE_data[RMSE_data$method == method_name, "RMSE"] <- rmse
}

# View the results
print(RMSE_data)

# 14. TPR-FPP curves ----
# set FDR to three levels : 0.01, 0.05 and 0.1
pr_results_list <- list()
FDR_thresholds <- c(0, 0.01, 0.05, 0.1, 1.0)
pval_cols <- c("adj.P.Val", "fisherAdjPval", "adj.pvalue", "adj_pval")


for (i in seq_along(method_list)) {
  method_name <- names(method_list)[i]  
  method_df <- method_list[[i]]
  recall <- precision <- c()
  for (percent in FDR_thresholds) {
    
    # Define conditions for each species
    is_human <- grepl("HUMAN", rownames(method_df))
    is_ecoli <- grepl("ECOLI", rownames(method_df))
    is_yeast <- grepl("YEAST", rownames(method_df))
    
    
    # Significant results
    if (method_name == "missBayes") {
      significant <- (method_df$pLtROPE >= (100 - 100 * percent)) | 
        (method_df$pGtROPE >= (100 - 100 * percent))
    } else {
      pcol <- pval_cols[pval_cols %in% colnames(method_df)][1]
      
      if (is.na(pcol)) {
        stop(paste("No p-value column found for", method_name))
      }
      
      significant <- method_df[[pcol]] <= percent
    }
    
    # fold change
    right_direction <- (grepl("ECOLI", rownames(method_df)) & method_df$logFC > 0) |
      (grepl("YEAST", rownames(method_df)) & method_df$logFC < 0)
    
    # Compute TP, FP, TN, FN
    TP <- sum((is_ecoli  & significant & right_direction) |
                (is_yeast  & significant & right_direction), na.rm = TRUE)
    TN <- sum(is_human & !significant, na.rm = TRUE)
    FP <- sum((is_human & significant) | (is_yeast  & significant & !right_direction) | (is_ecoli  & significant & !right_direction), na.rm = TRUE)
    FN <- sum((is_ecoli & !significant) | (is_yeast & !significant), na.rm = TRUE)
    
    # Metrics
    rec <- TP / (TP + FN)
    prec <- TP / (TP + FP)
    
    if (is.nan(prec)) prec <- 1
    
    recall <- c(recall, rec)
    precision <- c(precision, prec)
  }
  pr_data <- data.frame(
    precision = precision,
    recall = recall,
    Threshold = FDR_thresholds
  )
  
  pr_results_list[[method_name]] <- pr_data
}
# plot curve, show FDR thresholds on the plot
pr_combined <- do.call(rbind, lapply(names(pr_results_list), function(name) {
  df <- pr_results_list[[name]]
  df$method <- name
  return(df)
}))
pr_combined$FPP <- 1 - pr_combined$precision

# get label positions
label_data <- pr_combined %>%
  group_by(method) %>%
  arrange(FPP) %>%
  dplyr::slice(n() - 1)
label_data <- label_data %>%
  mutate(FPP_label = FPP + 0.02)

pr_combined$method <- factor(pr_combined$method, levels = c(
  "MSstats",
  "MSqRob2",
  "missBayes",
  "DEqMS",
  "proDA",
  "limpa",
  "DreamAI_limma",
  "mice_limma",
  "msImpute_limma",
  "ND_limma",
  "limma"
))


ggplot(pr_combined, aes(x = FPP, y = recall, color = method, group = method)) +
  
  # step curve (important change)
  geom_step(direction = "hv", linewidth = 1.0) +
  
  # optional points
  geom_point(size = 2) +
  
  # direct labels instead of legend
  geom_text_repel(
    data = label_data,
    aes(x = FPP_label, y = recall, label = method),
    direction = "y",              # keep vertical alignment clean
    hjust = 0,                   # left align text
    size = 4,
    segment.color = "grey50",
    segment.size = 0.6,
    min.segment.length = 0,
    box.padding = 0.3,
    point.padding = 0.2,
    force = 10,
    show.legend = FALSE
  ) +
  
  # reference FDR thresholds
  geom_vline(xintercept = c(0.01, 0.05, 0.1),
             linetype = "dashed",
             color = "grey60") +
  
  # zoom into region of interest
  coord_cartesian(xlim = c(0, 0.15), ylim = c(0, 1)) +
  
  labs(
    x = "False Positive Proportion",
    y = "True Positive Rate"
  ) +
  
  scale_color_manual(
    values = c(
      MSstats = "#1b9e77",
      MSqRob2 = "#d95f02",
      missBayes = "#7570b3",
      DEqMS = "#e7298a",
      proDA = "#66a61e",
      limpa = "#e6ab02",
      DreamAI_limma = "#a6761d",
      mice_limma = "#666666",
      msImpute_limma = "#1f78b4",
      ND_limma = "#33a02c",
      limma = "#fb9a99"
    ),
    labels = c(
      MSstats = "MSstats (15.0%, 25.3%, 31.2%)",
      MSqRob2 = "MSqRob2 (3.3%, 12.4%, 20.5%)",
      missBayes = "missBayes (0.7%, 2.7%, 5.3%)",
      DEqMS = "DEqMS (32.6%, 41.6%, 45.8%)",
      proDA = "proDA (30.9%, 39.3%, 42.9%)",
      limpa = "limpa (0.4%, 3.0%, 6.0%)",
      DreamAI_limma = "DreamAI_limma (0.7%, 2.0%, 3.4%)",
      mice_limma = "mice_limma (6.7%, 9.1%, 10.5%)",
      msImpute_limma = "msImpute_limma (1.0%, 3.4%, 5.7%)",
      ND_limma = "ND_limma (0.4%, 1.2%, 2.4%)",
      limma = "limma (0.6%, 1.7%, 3.0%)"
    )
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 1.0),
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(0.25, "cm"),
    legend.position = "right"   
  )


# 15. box plots ----
ecoli_logFC_list <- list()
for (method_name in names(method_list)) {
  
  method_df <- method_list[[method_name]]
  
  # species
  is_ecoli <- grepl("ECOLI", rownames(method_df))
  
  
  # subset
  selected <- is_ecoli
  
  # extract logFC
  ecoli_logFC_list[[method_name]] <- method_df$logFC[selected]
}

ecoli_logFC_df <- do.call(rbind, lapply(names(ecoli_logFC_list), function(m) {
  data.frame(
    method = m,
    logFC = ecoli_logFC_list[[m]]
  )
}))

# Compute medians
medians <- ecoli_logFC_df %>%
  group_by(method) %>%
  summarise(median_logFC = median(logFC, na.rm = TRUE))

ggplot(ecoli_logFC_df, aes(x = method, y = logFC, fill = method)) +
  geom_boxplot(width = 0.9, alpha = 0.8) +
  coord_cartesian(ylim = c(-3, 5)) +
  geom_text(
    data = medians,
    aes(x = method, y = median_logFC, label = round(median_logFC, 2)),
    y = 4.9,                      # move label slightly above median line
    fontface = "bold",
    size = 4
  ) +
  labs(x = "Methods", y = "logFC") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  ) +
  scale_fill_brewer(palette = "Set3") +
  geom_hline(yintercept = 1, color = "red", linewidth = 1, linetype = "dashed")

# 16. Unique proteins ----
# index of unique proteins
index_unique_B <- rownames(log2norm.df)[
  xor(
    apply(log2norm.df[, 1:5], 1, function(x) all(is.na(x))),
    apply(log2norm.df[, 6:10], 1, function(x) all(is.na(x)))
  )
]

index_unique_A <- rownames(log2norm.df)[
  xor(
    apply(log2norm.df[, 11:15], 1, function(x) all(is.na(x))),
    apply(log2norm.df[, 16:20], 1, function(x) all(is.na(x)))
  )
]
# calculate true positive rate, precision and false sign rate for unique proteins
method_list <- list(
  MSqRob2 = msqRob2.results_701020,
  proDA = proDA.results_701020,
  limpa = limpa_result_701020,
  DreamAI_limma = DreamAI.results_701020,
  mice_limma = mice.results_701020,
  msImpute_limma = msImpute.results_701020,
  ND_limma = ND.results_701020,
  missBayes = missBayes.results_701020[[1]]
)

threshold <- 0.05
pval_cols <- c("adj.P.Val", "fisherAdjPval", "adj.pvalue", "adj_pval")

unique_metrics_list <- list()

for (method_name in names(method_list)) {
  
  method_df <- method_list[[method_name]]
  
  # restrict to unique proteins
  method_df <- method_df[index_unique_A, , drop = FALSE]
  
  # species
  is_human <- grepl("HUMAN", rownames(method_df))
  is_ecoli <- grepl("ECOLI", rownames(method_df))
  is_yeast <- grepl("YEAST", rownames(method_df))
  
  positives <- is_ecoli | is_yeast
  
  # significance
  if (method_name == "missBayes") {
    significant <- (method_df$pLtROPE >= (100 - 100 * threshold)) | 
      (method_df$pGtROPE >= (100 - 100 * threshold))
  } else {
    pcol <- pval_cols[pval_cols %in% colnames(method_df)][1]
    
    if (is.na(pcol)) {
      stop(paste("No p-value column found for", method_name))
    }
    
    significant <- method_df[[pcol]] <= threshold
  }
  
  # direction
  right_direction <- (is_ecoli & method_df$logFC > 0) |
    (is_yeast & method_df$logFC < 0)
  
  wrong_direction <- (is_ecoli & method_df$logFC < 0) |
    (is_yeast & method_df$logFC > 0)
  
  # counts
  TP <- sum(positives & significant & right_direction, na.rm = TRUE)
  FP <- sum(is_human & significant, na.rm = TRUE)
  FN <- sum(positives & !significant, na.rm = TRUE)
  
  # metrics
  recall <- TP / (TP + FN)
  precision <- TP / (TP + FP)
  
  if (is.nan(precision)) precision <- 1
  
  # FSR 
  FSR <- sum(positives & significant & wrong_direction, na.rm = TRUE) /
    sum(positives, na.rm = TRUE)
  
  # store
  unique_metrics_list[[method_name]] <- data.frame(
    method = method_name,
    recall = recall,
    precision = precision,
    FSR = FSR
  )
}

unique_metrics_df <- do.call(rbind, unique_metrics_list)
unique_metrics_df$precision <- unique_metrics_df$precision * 100
unique_metrics_df$recall <- unique_metrics_df$recall * 100
unique_metrics_df$FSR <- unique_metrics_df$FSR * 100

# draw bar plot
plot_df <- unique_metrics_df %>%
  pivot_longer(
    cols = c(recall, precision, FSR),
    names_to = "metric",
    values_to = "value"
  )
ggplot(plot_df, aes(x = method, y = value, fill = metric)) +
  
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.3
  ) +
  
  # nicer colors
  scale_fill_manual(
    values = c(
      recall = "#1b9e77",
      precision = "#d95f02",
      FSR = "#7570b3"
    ),
    labels = c(
      recall = "Recall",
      precision = "Precision",
      FSR = "False Sign Rate"
    )
  ) +
  
  # axis labels
  labs(
    x = NULL,
    y = "Percentage (%)",
    fill = NULL
  ) +
  
  # y-axis range
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  
  # clean theme
  theme_classic(base_size = 13) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11,
      face = "bold"
    ),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13, face = "bold"),
    
    legend.position = "top",
    legend.text = element_text(size = 11),
    
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_line(linewidth = 1)
  )


# RMSE calculation for all unique proteins
RMSE_data <- data.frame(method = names(method_list), RMSE = NA)

# Loop through each named element in the list
for (i in seq_along(method_list)) {
  method_name <- names(method_list)[i]  
  method_df <- method_list[[i]]
  method_df <- method_df[index_unique_A, , drop = FALSE]
  
  # Add the true value column based on protein names
  method_df$trueVal <- ifelse(grepl("HUMAN", rownames(method_df)), 0,
                              ifelse(grepl("ECOLI", rownames(method_df)), 1,
                                     ifelse(grepl("YEAST", rownames(method_df)), -1, NA)))
  method_df$logFC[!is.finite(method_df$logFC)] <- NA
  
  # Calculate RMSE for this method
  rmse <- sqrt(mean((method_df$trueVal - method_df$logFC)^2, na.rm = TRUE))
  
  # Assign the result to the correct row in RMSE_data
  RMSE_data[RMSE_data$method == method_name, "RMSE"] <- rmse
}
RMSE_data$RMSE <- round(RMSE_data$RMSE, 3)
# View the results
print(RMSE_data)

# 17. how logFC distribute with different number of peptide counts ----
pep_count_map <- pg[, c("Protein.Names", "Peptide.Count")]
pep_count_map$Peptide.Count <- as.numeric(pep_count_map$Peptide.Count)
pep_count_vec <- setNames(pep_count_map$Peptide.Count, pep_count_map$Protein.Names)

logFC_list <- list()

for (method_name in names(method_list)) {
  
  method_df <- method_list[[method_name]]
  
  # species
  is_ecoli <- grepl("ECOLI", rownames(method_df))
  is_yeast <- grepl("YEAST", rownames(method_df))
  
  
  # map peptide count using rownames
  pep_count <- pep_count_vec[rownames(method_df)]
  
  df <- data.frame(
    method = method_name,
    logFC = method_df$logFC,
    pep_count = pep_count,
    species = ifelse(is_ecoli, "ECOLI",
                     ifelse(is_yeast, "YEAST", NA))
  )
  
  # keep only ecoli + yeast + significant
  df <- df[is_yeast, ]
  
  logFC_list[[method_name]] <- df
}

logFC_df <- bind_rows(logFC_list)

table(logFC_df$pep_count)


peptide_counts <- c(2:10)

plot_list <- list()

for (pc in peptide_counts) {
  
  df_sub <- logFC_df[logFC_df$pep_count == pc, ]
  medians <- df_sub %>%
    group_by(method) %>%
    summarise(median_logFC = median(logFC, na.rm = TRUE))
  
  p <- ggplot(df_sub, aes(x = method, y = logFC, fill = method)) +
    
    geom_boxplot(width = 0.9, alpha = 0.8) +
    
    coord_cartesian(ylim = c(-5, 5)) +
    
    geom_text(
      data = medians,
      aes(x = method, y = median_logFC, label = round(median_logFC, 2)),
      y = 4.6,                      # move label slightly above median line
      fontface = "bold",
      size = 4
    ) +
    
    # horizontal reference lines
    geom_hline(yintercept = -1, color = "red", linewidth = 1, linetype = "dashed") +
    
    labs(
      title = paste("Peptide count =", pc),
      x = "Methods",
      y = "logFC"
    ) +
    
    theme_minimal() +
    
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
    ) +
    
    scale_fill_brewer(palette = "Set3")
  
  plot_list[[as.character(pc)]] <- p
}
