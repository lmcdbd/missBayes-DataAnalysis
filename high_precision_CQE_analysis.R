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
pg<-fread("singlecell_level/report.pg_matrix.tsv")
pr<-fread("singlecell_level/report.pr_matrix.tsv")
raw <- fread("singlecell_level/report.tsv")
t <- table(unique(pr[,c('Protein.Group','Stripped.Sequence')])$Protein.Group)
pg$Peptide.Count <- t[match(pg$Protein.Group,names(t))]
metadata <- fread("singlecell_level/metadata.txt", header = TRUE)
#Filter by protein groups with <2 peptides, remove contaminants, and log2 transform the data
df.prot <- filter(pg, Peptide.Count > 1)
df.prot2 <- df.prot[!grepl("Cont_",df.prot$Protein.Group),]
df.prot.all<- df.prot2[, 7:16]
protein.matrix_all <- log2(as.matrix(df.prot.all))
log2all.df <- as.data.frame(protein.matrix_all)
rownames(log2all.df) <- df.prot2$Protein.Names
# normalization
log2norm.df <- normalizeMedianValues(log2all.df)
plotDensities(log2norm.df, legend = FALSE)

# missing proportion and number of unique proteins
sum(is.na(log2all.df)) / nrow(log2all.df) / ncol(log2all.df)

# 1. MSstats ----
# convert DIANN output to MSstats raw data
# remove blanks
metadata_msstats <- metadata
colnames(metadata_msstats)[2] <- "Run"
raw <- raw[!grepl("BlankPost", raw$File.Name), ]
MSstats_data <- DIANNtoMSstatsFormat(raw, metadata_msstats, removeFewMeasurements = FALSE, 
                                     quantificationColumn = "FragmentQuantRaw" )

# data pre-processing and summarization
MSstats_results <- dataProcess(MSstats_data, normalization = "equalizeMedians", summaryMethod = "TMP", MBimpute = TRUE)

MSstats_results$ProteinLevelData <- MSstats_results$ProteinLevelData[!is.na(MSstats_results$ProteinLevelData$GROUP), ]
MSstats_results$FeatureLevelData <- MSstats_results$FeatureLevelData[!is.na(MSstats_results$FeatureLevelData$GROUP), ]

# differential abundance analysis for group comparison
levels <- levels(as.factor(MSstats_results$ProteinLevelData$GROUP))
comparison <- matrix(0, nrow = 1, ncol = length(levels))
colnames(comparison) <- levels
rownames(comparison) <- "G70_10_20 - G70_20_10"
comparison[1, "G70_10_20"] <- 1
comparison[1, "G70_20_10"] <- -1
testResultOneComparison <- groupComparison(contrast.matrix= comparison, data= MSstats_results)
MSstats.results <- testResultOneComparison$ComparisonResult

# 2. MSqRob2 ----
# remove blanks
raw <- raw[!grepl("BlankPost", raw$File.Name), ]
raw[["File.Name"]] <- raw[["Run"]]

#Assumes report has already been filtered according to Q-values, removing contaminants, shared peptides, etc as required
raw_wide <- raw[,c(2:6, 14:16, 28)] %>% pivot_wider(names_from = Run, values_from = Precursor.Normalised, values_fn = sum)

metadata[["Run"]] <- sort(unique(raw[["File.Name"]]))
colnames(metadata)[2] <- "quantCols"

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

#Fit ridge regression peptide-level model (works but may error on single peptide proteins or proteins with lots of missing values)
qf <- msqrobAggregate(qf, i = "log2norm", fcol = "Protein.Names", formula = ~Condition + (1|Modified.Sequence), ridge = TRUE)
#Fit Hurdle model (to summarized protein quant) 
qf_agg <- aggregateFeatures(qf, i = "log2norm", fcol = "Protein.Names", name = "log2proteinQuant")
qf_hurdle <- msqrobHurdle(qf_agg, i = "log2proteinQuant", formula = ~Condition)
# Check names in the Intensity component
coef <- getCoef(rowData(qf_hurdle[["log2proteinQuant"]])$msqrobHurdleIntensity[[1]])
contrast <- makeContrast(contrasts = "-ConditionG70_20_10 = 0", parameterNames = names(coef))

msqrob2_result <- hypothesisTestHurdle(object = qf_hurdle, contrast = contrast, i = "log2proteinQuant")

msqRob2.results <- rowData(msqrob2_result[["log2proteinQuant"]])[["hurdle_-ConditionG70_20_10"]]

# 3. DEqMS ----
# median normalization
pep.count.table <- data.frame(count = as.vector(pg$Peptide.Count),
                              row.names = pg$Protein.Names)
class <- as.factor(metadata$Condition)
design <- model.matrix(~0+class) # fitting without intercept
colnames(design) <- sub("^class", "", colnames(design))
fit1 <- lmFit(log2norm.df,design = design)
cont <- makeContrasts(G70_10_20-G70_20_10, levels = levels(class))
fit2 <- contrasts.fit(fit1,contrasts = cont)
fit3 <- eBayes(fit2)

fit3$count <- pep.count.table[rownames(fit3$coefficients),"count"]

#check the values in the vector fit3$count
#if min(fit3$count) return NA or 0, troubleshoot the error first
min(fit3$count)
fit4 <- spectraCounteBayes(fit3)
DEqMS.results <- outputResult(fit4,coef_col = 1)

# 4. proDA ----
fit <- proDA(log2norm.df, design = design) # fit dropout model
proDA.results <- test_diff(fit, `G70_10_20` - `G70_20_10`)

# 5. limpa ----
y.prec <- readDIANN(file = "singlecell_level/report.tsv",
                    run.column = "File.Name", q.columns = c("Q.Value", "Protein.Q.Value"))
dpcest <- dpcCN(y.prec) # estimate dpc 
plotDPC(dpcest)
y <- dpcQuant(y.prec, dpc = dpcest) # summarization
fit <- dpcDE(y, design, plot = TRUE)
fit <- eBayes(fit)
cont <- makeContrasts(G70_10_20 - G70_20_10, levels = design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)
limpa.result <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
rownames(limpa.result) <- limpa.result$Protein.Names


# 7. msImpute+limma ----
# Convert to log2 
pep_log2 <- as.matrix(log2(pr[, 13:22]))
rownames(pep_log2) <- pr$Stripped.Sequence

pep_log2[is.infinite(pep_log2)] <- NA  # Ensure 0s become NAs
# median normalization
pep_log2 <- normalizeMedianValues(pep_log2)
# msImpute
msImpute.pep <- msImpute(pep_log2, method = "v2-mnar", design = design, relax_min_obs = TRUE)
# robust summarization
# 1. Create a SummarizedExperiment for the peptide data
colnames(msImpute.pep) <- metadata$quantCols

se_pep <- SummarizedExperiment(assays = list(log2quant = msImpute.pep),
                               colData = metadata) 
rownames(se_pep) <- make.unique(rownames(se_pep))
# 2. Wrap it in a QFeatures object
qf_imputed <- QFeatures(list(pep_imputed = se_pep))

# 3. Create a mapping (using original 'raw' DIA-NN data)
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
# Convert to factor so msqrob knows they are groups, not numbers
colData(qf_imputed)$Condition <- as.factor(colData(qf_imputed)$Condition)
# Verify - you should now see the G-names, not 1, 2, 3...
colData(qf_imputed)$Condition
# 5. Aggregate to Protein Level

cd <- colData(qf_imputed)

cd$Condition <- as.character(cd$Condition)

colData(qf_imputed) <- cd
colData(qf_imputed[["pep_imputed"]])$Condition <- cd$Condition

qf_imputed <- aggregateFeatures(qf_imputed, 
                                i = "pep_imputed", 
                                fcol = "Protein.Group", 
                                name = "protein_imputed", 
                                fun = MsCoreUtils::robustSummary)

# extract the aggregated protein-level data
msImpute.results <- assay(qf_imputed[["protein_imputed"]])
# statistical inference 
fit <- lmFit(msImpute.results, design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- treat(fit2, lfc=log2(1.15), trend = TRUE)
msImpute.results <- topTreat(fit2, coef = "G70_10_20 - G70_20_10", adjust.method = "BH", number = Inf, confint = TRUE)
# replace row names with protein names
# Match row names of df1 to Protein.Group in df2, then pull the corresponding Protein.Names
rownames(msImpute.results) <- raw$Protein.Names[match(rownames(msImpute.results), raw$Protein.Group)]

# 8. DreamAI+limma ----
DreamAI_imputed <- DreamAI(log2norm.df)
DreamAI_imputed <- DreamAI_imputed$Ensemble

# limma inference
fit <- lmFit(DreamAI_imputed, design)
fit2 <- contrasts.fit(fit, cont)
keep <- !is.na(fit2$Amean)
fit2 <- fit2[keep, ]
fit2 <- treat(fit2, lfc = log2(1.15), trend = TRUE)
DreamAI.results <- topTreat(fit2, coef = "G70_10_20 - G70_20_10", adjust.method = "BH", number = Inf, confint = TRUE)

# 9. mice-pmm+limma ----
colnames(log2norm.df) <- metadata$quantCols
mice_imputed <- mice(log2norm.df, method = "pmm", seed = 456)
mice_pooled <- Reduce("+", lapply(1:5, function(i) complete(mice_imputed, i))) / 5
mice_imputed <- mice_pooled
# limma inference
fit <- lmFit(mice_imputed, design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- treat(fit2, lfc=log2(1.15), trend = TRUE)
mice.results <- topTreat(fit2, coef = "G70_10_20 - G70_20_10", adjust.method = "BH", number = Inf, confint = TRUE)

# 10. missBayes ----
group <- as.factor(metadata$Condition)
comparison <- makeContrasts("G70_10_20 - G70_20_10", levels = levels(group))
set.seed(746)
missBayes.results <- BayesMissingModel(values = log2norm.df, groups = group, comparisons = comparison, parallel = TRUE,
                             n.adapt = 500, burn.in = 500, n.iter = 5000, n.chains = 2, mcmcDiag = TRUE)

# 11. ND+limma ----
ND_imputed <- impute_normal(log2norm.df)
# limma inference
fit <- lmFit(ND_imputed, design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- treat(fit2, lfc=log2(1.15), trend = TRUE)
ND.results <- topTreat(fit2, coef = "G70_10_20 - G70_20_10", adjust.method = "BH", number = Inf, confint = TRUE)


# 12. limma ----
fit <- lmFit(log2norm.df, design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- treat(fit2, lfc=log2(1.15), trend = TRUE)
limma.results <- topTreat(fit2, coef = "G70_10_20 - G70_20_10", adjust.method = "BH", number = Inf)

# Keep proteins present in all methods for fair comparison
rownames(MSstats.results) <- MSstats.results$Protein
rownames(limpa.result) <- limpa.result$Protein.Names
rownames(proDA.results) <- proDA.results$name

MSstats.results <- MSstats.results[!is.infinite(MSstats.results$log2FC), ]

common_proteins <- Reduce(intersect, list(rownames(MSstats.results), rownames(msqRob2.results), rownames(DEqMS.results), rownames(proDA.results), rownames(limpa.result), 
                                          rownames(DreamAI.results), rownames(mice.results), rownames(msImpute.results), rownames(ND.results), rownames(missBayes.results[[1]]), rownames(limma.results)))
                          
# 13. RMSE calculation ----
# rename logFC columns before calculation
missBayes.results[[1]] <- missBayes.results[[1]] %>% dplyr::rename(logFC = Median)
MSstats.results <- MSstats.results %>% dplyr::rename(logFC = log2FC)
proDA.results <- proDA.results%>% dplyr::rename(logFC = diff)


method_list <- list(
  MSstats = MSstats.results[common_proteins, ],
  MSqRob2 = msqRob2.results[common_proteins, ],
  DEqMS = DEqMS.results[common_proteins, ],
  proDA = proDA.results[common_proteins, ],
  limpa = limpa.result[common_proteins, ],
  DreamAI_limma = DreamAI.results[common_proteins, ],
  mice_limma = mice.results[common_proteins, ],
  msImpute_limma = msImpute.results[common_proteins, ],
  ND_limma = ND.results[common_proteins, ],
  missBayes = missBayes.results[[1]][common_proteins, ],
  limma = limma.results[common_proteins, ]
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

# 14. TPR_FPP curves ----
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
  dplyr::slice(n() - 2)
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
    force = 15,
    show.legend = FALSE
  ) +
  
  # reference FDR thresholds
  geom_vline(xintercept = c(0.01, 0.05, 0.1),
             linetype = "dashed",
             color = "grey60") +
  
  # zoom into region of interest
  coord_cartesian(xlim = c(0, 0.3), ylim = c(0, 0.5)) +
  
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
      MSstats = "MSstats (0.0%, 34.9%, 58.3%)",
      MSqRob2 = "MSqRob2 (7.8%, 13.2%, 27.4%)",
      missBayes = "missBayes (2.6%, 9.4%, 15.0%)",
      DEqMS = "DEqMS (5.5%, 17.1%, 28.7%)",
      proDA = "proDA (6.9%, 12.3%, 23.7%)",
      limpa = "limpa (3.7%, 16.6%, 35.3%)",
      DreamAI_limma = "DreamAI_limma (0.0%, 4.9%, 9.4%)",
      mice_limma = "mice_limma (3.9%, 10.2%, 16.7%)",
      msImpute_limma = "msImpute_limma (4.4%, 11.2%, 15.0%)",
      ND_limma = "ND_limma (2.5%, 5.7%, 12.3%)",
      limma = "limma (4.0%, 7.6%, 12.5%)"
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

ggsave("/Users/mengchunli/Documents/work/Paper revision/paper documents_new/figures/PR_curve_CQE_701020.png", width = 6, height = 4)

# 15. box plots ----
yeast_logFC_list <- list()
for (method_name in names(method_list)) {
  
  method_df <- method_list[[method_name]]
  
  # species
  is_yeast <- grepl("YEAST", rownames(method_df))
  
  
  # subset
  selected <- is_yeast
  
  # extract logFC
  yeast_logFC_list[[method_name]] <- method_df$logFC[selected]
}

yeast_logFC_df <- do.call(rbind, lapply(names(yeast_logFC_list), function(m) {
  data.frame(
    method = m,
    logFC = yeast_logFC_list[[m]]
  )
}))

# Compute medians
medians <- yeast_logFC_df %>%
  group_by(method) %>%
  summarise(median_logFC = median(logFC, na.rm = TRUE))

ggplot(yeast_logFC_df, aes(x = method, y = logFC, fill = method)) +
  geom_boxplot(width = 0.9, alpha = 0.8) +
  coord_cartesian(ylim = c(-3.5, 2)) +
  geom_text(
    data = medians,
    aes(x = method, y = median_logFC, label = round(median_logFC, 2)),
    y = 1.9,                      # move label slightly above median line
    fontface = "bold",
    size = 4
  ) +
  labs(x = "Methods", y = "logFC") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold", margin = ggplot2::margin(t = 10)),
    axis.title.y = element_text(size = 14, face = "bold", margin = ggplot2::margin(r = 10)),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  ) +
  scale_fill_brewer(palette = "Set3") +
  geom_hline(yintercept = -1, color = "red", linewidth = 1, linetype = "dashed")

# 16. Unique proteins ----
# index of unique proteins
index_unique <- rownames(log2norm.df)[
  xor(
    apply(log2norm.df[, 1:5], 1, function(x) all(is.na(x))),
    apply(log2norm.df[, 6:10], 1, function(x) all(is.na(x)))
  )
]

# calculate true positive number, false positive number for unique proteins
method_list <- list(
  MSqRob2 = msqRob2.results,
  proDA = proDA.results,
  limpa = limpa.result,
  DreamAI_limma = DreamAI.results,
  mice_limma = mice.results,
  msImpute_limma = msImpute.results,
  ND_limma = ND.results,
  missBayes = missBayes.results[[1]]
)

threshold <- 0.05
pval_cols <- c("adj.P.Val", "fisherAdjPval", "adj.pvalue", "adj_pval")

unique_metrics_list <- list()

for (method_name in names(method_list)) {
  
  method_df <- method_list[[method_name]]
  
  # restrict to unique proteins
  method_df <- method_df[index_unique, , drop = FALSE]
  
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
    FSR = FSR,
    TP = TP,
    FP = FP
  )
}

unique_metrics_df <- do.call(rbind, unique_metrics_list)
unique_metrics_df$precision <- unique_metrics_df$precision * 100
unique_metrics_df$recall <- unique_metrics_df$recall * 100
unique_metrics_df$FSR <- unique_metrics_df$FSR * 100
unique_metrics_df <- unique_metrics_df[, c(1, 5:6)]
# draw bar plot
plot_df <- unique_metrics_df %>%
  pivot_longer(
    cols = c(TP,FP),
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
      TP = "#1b9e77",
      FP = "#d95f02"
    ),
    labels = c(
      TP = "True positive",
      FP = "False positive"
    )
  ) +
  
  # axis labels
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  
  # y-axis range
  scale_y_continuous(
    limits = c(0, 15),
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
    axis.ticks = element_line(linewidth = 1),
    plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 20)
  )


# RMSE calculation for all unique proteins
RMSE_data <- data.frame(method = names(method_list), RMSE = NA)

# Loop through each named element in the list
for (i in seq_along(method_list)) {
  method_name <- names(method_list)[i]  
  method_df <- method_list[[i]]
  method_df <- method_df[index_unique, , drop = FALSE]
  
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

