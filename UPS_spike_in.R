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
pg<-fread("UPS_spikein/proteinGroups.txt")
pr<-fread("UPS_spikein/peptides.txt")
t <- table(unique(pr[,c('Proteins','Sequence')])$'Proteins')
pg$Peptide.Count <- t[match(pg$'Protein IDs',names(t))]
metadata <- fread("UPS_spikein/metadata.txt", header = TRUE)
#Filter by protein groups with <2 peptides, remove contaminants, and log2 transform the data
df.prot <- filter(pg, Peptide.Count > 1)
df.prot2 <- df.prot[!grepl("^CON_",df.prot$'Protein IDs'),]
df.prot.all<- df.prot2[, 305:344]
log2all.df <- as.matrix(log2(df.prot.all))
log2all.df[is.infinite(log2all.df) & log2all.df < 0] <- NA

log2norm.df <- normalizeMedianValues(log2all.df)
log2norm.df <- as.data.frame(log2norm.df)
rownames(log2norm.df) <- df.prot2$`Protein IDs`
log2norm.df <- log2norm.df[rowSums(!is.na(log2norm.df)) > 0, ]

colnames(metadata) <- c("Run", "rename", "Condition")
metadata$BioReplicate <- c(1:40)

# comparisons to be made
contrasts <- c("50amol - 500amol", "100amol - 500amol", "250amol - 500amol",
                  "1fmol - 500amol", "5fmol - 500amol")

# 1. MSstats ----
# convert MaxQuant output to MSstats raw data
infile <- read.table("UPS_spikein/evidence.txt", 
                     sep="\t", header=TRUE, fill = TRUE)

run_mapping <- unique(infile[, c("Experiment", "Raw.file")])
metadata <- merge(metadata, run_mapping, 
                  by.x = "rename", 
                  by.y = "Experiment", 
                  all.x = TRUE)
metadata$Run <- metadata$Raw.file

raw <- MaxQtoMSstatsFormat(evidence=infile, 
                           annotation=metadata, 
                           proteinGroups=pg,)
# data pre-processing and summarization
MSstats_results <- dataProcess(raw, normalization = "equalizeMedians", summaryMethod = "TMP", MBimpute = TRUE)

# differential abundance analysis for group comparison
MSstats.results <- list()
levels <- levels(as.factor(MSstats_results$ProteinLevelData$GROUP))
for (contrast in contrasts) {
  # Seperate the contrast string into two group names
  parts <- strsplit(contrast, " - ")[[1]]
  group1 <- trimws(parts[1])  
  group2 <- trimws(parts[2]) 
  
  # Build contrast matrix
  comparison <- matrix(0, nrow = 1, ncol = length(levels))
  colnames(comparison) <- levels
  rownames(comparison) <- contrast
  
  comparison[1, group1] <- 1
  comparison[1, group2] <- -1
  
  # Run comparison
  testResultOneComparison <- groupComparison(contrast.matrix = comparison, 
                                             data = MSstats_results)
  
  MSstats.results[[contrast]] <- testResultOneComparison$ComparisonResult
}

# 2. MSqRob2 ----
pr <- as.data.frame(pr)
quantCols <- grep("^Intensity", colnames(pr), value = TRUE)[2:41] # Remove sum intensity column
metadata$quantCols <- quantCols
pr[quantCols] <- lapply(pr[quantCols], as.numeric)
qf <- readQFeatures(
  assayData = pr,
  fnames = "Sequence",
  colData = metadata,
  quantCols = quantCols,
  name = "peptideRaw"
)

qf <- logTransform(qf, base = 2, i = "peptideRaw", name = "log2quant")
Protein_filter <- rowData(qf[["log2quant"]])$Proteins %in% smallestUniqueGroups(rowData(qf[["log2quant"]])$Proteins)
qf <- qf[Protein_filter,,]
qf <- filterFeatures(qf, ~ Reverse != "+")
qf <- filterFeatures(qf, ~ Potential.contaminant != "+")
# normalization
qf<- normalize(qf, i = "log2quant", name = "log2norm", method = "center.median")
assay(qf[["log2norm"]])[!is.finite(assay(qf[["log2norm"]]))] <- NA # Convert -Inf to NA after log transformation
limma::plotDensities(assay(qf[["log2quant"]]), legend = FALSE)

# Filter out proteins with less than 2 values, and then filter out proteins with more than 90% missing values
keep <- rowSums(!is.na(assay(qf[["log2norm"]]))) >= 2
qf <- qf[keep, , ]
qf <- filterNA(
  qf,
  i = "log2norm",
  pNA = 0.9
)

#Fit ridge regression peptide-level model (works but may error on single peptide proteins or proteins with lots of missing values)
qf <- msqrobAggregate(qf, i = "log2norm", fcol = "Proteins", formula = ~Condition + (1|Sequence), ridge = TRUE)

#Fit Hurdle model (to summarized protein quant) 
qf_agg <- aggregateFeatures(qf, i = "log2norm", fcol = "Protein.Names", name = "log2proteinQuant")
qf_hurdle <- msqrobHurdle(qf_agg, i = "log2proteinQuant", formula = ~Condition)

# loop over contrasts and perform hypothesis testing
msqRob2.results <- list()
coef <- getCoef(rowData(qf_hurdle[["log2proteinQuant"]])$msqrobHurdleIntensity[[1]])
ref_level <- "100amol"  # the missing one from names(coef)

for (contrast in contrasts) {
  parts <- strsplit(contrast, " - ")[[1]]
  group1 <- trimws(parts[1])
  group2 <- trimws(parts[2])
  
  # Build contrast string handling reference level
  if (group1 == ref_level) {
    # e.g. 100amol - 500amol = 0 - Condition500amol
    contrast_str <- paste0("-Condition", group2, " = 0")
    contrast_name <- paste0("hurdle_-Condition", group2)
  } else if (group2 == ref_level) {
    # e.g. 50amol - 100amol = Condition50amol - 0
    contrast_str <- paste0("Condition", group1, " = 0")
    contrast_name <- paste0("hurdle_Condition", group1)
  } else {
    # Neither group is reference, standard contrast
    contrast_str <- paste0("Condition", group1, " - Condition", group2, " = 0")
    contrast_name <- paste0("hurdle_Condition", group1, " - Condition", group2)
  }
  
  contrast_matrix <- makeContrast(contrasts = contrast_str,
                                  parameterNames = names(coef))
  
  msqrob2_result <- hypothesisTestHurdle(object = qf_hurdle,
                                         contrast = contrast_matrix,
                                         i = "log2proteinQuant")
  
  msqRob2.results[[contrast]] <- rowData(msqrob2_result[["log2proteinQuant"]])[[contrast_name]]
}


# 3. DEqMS ----
# median normalization
pep.count.table <- data.frame(count = as.vector(pg$Peptide.Count),
                              row.names = pg$`Protein IDs`)
class <- as.factor(metadata$Condition)
design <- model.matrix(~0+class) # fitting without intercept
colnames(design) <- sub("^class", "", colnames(design))
fit1 <- lmFit(log2norm.df,design = design)
DEqMS.results <- list()

for (contrast in contrasts) {
  parts <- strsplit(contrast, " - ")[[1]]
  group1 <- trimws(parts[1])
  group2 <- trimws(parts[2])
  
  # Make syntactically valid R names
  group1_valid <- make.names(group1)  # e.g. "50amol" -> "X50amol"
  group2_valid <- make.names(group2)  # e.g. "500amol" -> "X500amol"
  
  # Build contrast using valid names
  cont <- makeContrasts(
    contrasts = paste0(group1_valid, " - ", group2_valid),
    levels = make.names(levels(class))
  )
  
  # Fit and eBayes
  fit2 <- contrasts.fit(fit1, contrasts = cont)
  fit3 <- eBayes(fit2)
  
  # Add peptide counts
  fit3$count <- pep.count.table[rownames(fit3$coefficients), "count"]
  
  # Check counts
  if (is.na(min(fit3$count)) || min(fit3$count) == 0) {
    warning(paste("Count issue in contrast:", contrast))
    next
  }
  
  # DEqMS
  fit4 <- spectraCounteBayes(fit3)
  
  # Store result
  DEqMS.results[[contrast]] <- outputResult(fit4, coef_col = 1)
}

# 4. proDA ----
proDA.results <- list()
fit <- proDA(as.matrix(log2norm.df), design = design) # fit dropout model
for (contrast in contrasts) {
  parts <- strsplit(contrast, " - ")[[1]]
  group1 <- trimws(parts[1])
  group2 <- trimws(parts[2])
  
  # Build contrast expression as string and parse it
  contrast_expr <- paste0("`", group1, "` - `", group2, "`")
  
  # Store result
  proDA.results[[contrast]] <- test_diff(fit, contrast = eval(parse(text = contrast_expr)))
}


# 5. limpa ----

y.prec <- readMaxQuant(file = "UPS_spikein/peptides.txt")

dpcest <- dpcCN(y.prec) # estimate dpc 
plotDPC(dpcest)

y <- dpcQuant(y.prec, dpc = dpcest, protein.id = "Proteins") # summarization
design <- model.matrix(~0+Condition, data = metadata)
colnames(design) <- sub("^Condition", "", colnames(design))
fit <- dpcDE(y, design, plot = TRUE)
fit <- eBayes(fit)
limpa.results <- list()

for (contrast in contrasts) {
  parts <- strsplit(contrast, " - ")[[1]]
  group1 <- trimws(parts[1])
  group2 <- trimws(parts[2])
  
  # Make valid R names
  group1_valid <- make.names(group1)
  group2_valid <- make.names(group2)
  
  # Build contrast matrix
  cont <- makeContrasts(
    contrasts = paste0(group1_valid, " - ", group2_valid),
    levels = make.names(colnames(design))
  )
  
  # Fit contrasts and eBayes
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- eBayes(fit2)
  
  # Store result
  limpa.results[[contrast]] <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
}



# 7. msImpute+limma ----
# MaxLFQ peptide-level intensities
pep <- pr[, grepl("^Intensity", colnames(pr))]
rownames(pep) <- pr$Sequence

# Convert to log2 
pep_log2 <- as.matrix(log2(pep))
pep_log2[is.infinite(pep_log2)] <- NA  

# msImpute
design <- model.matrix( ~0+Condition, metadata )
pep_log2 <- normalizeMedianValues(pep_log2)
# median normalization after imputation
msImpute.pep <- normalizeMedianValues(msImpute.pep)

# robust summarization
# 1. Create a SummarizedExperiment for the peptide data
msImpute.pep <- msImpute.pep[, rownames(colData(qf))]
se_pep <- SummarizedExperiment(assays = list(log2quant = msImpute.pep),
                               colData = colData(qf)) # colData from MSqRob2 QFeature
# 2. Wrap it in a QFeatures object
qf_imputed <- QFeatures(list(pep_imputed = se_pep))

# Get unique mappings between sequences and protein groups
pep_to_prot <- unique(pr[, c("Sequence", "Proteins")])
setDF(pep_to_prot)
rownames(pep_to_prot) <- pep_to_prot$Sequence

# 4. Add this mapping to the rowData of new QFeatures object
rowData(qf_imputed[["pep_imputed"]]) <- pep_to_prot[rownames(msImpute.pep), ]

# Ensure the metadata row names match the colData row names
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$quantCols 

# Re-align and assign the actual character names
colData(qf_imputed)$Condition <- metadata[rownames(colData(qf_imputed)), "Condition"]
# Convert to factor 
colData(qf_imputed)$Condition <- as.factor(colData(qf_imputed)$Condition)

# 5. Aggregate to Protein Level
colData(qf_imputed[["pep_imputed"]]) <- colData(qf_imputed)

rd <- rowData(qf_imputed[["pep_imputed"]])
rd$Proteins <- sapply(strsplit(rd$Proteins, ";"), `[`, 1) # Take the first protein group for each peptide
rowData(qf_imputed[["pep_imputed"]]) <- rd

qf_imputed <- aggregateFeatures(qf_imputed, 
                                i = "pep_imputed", 
                                fcol = "Proteins", 
                                name = "protein_imputed", 
                                fun = MsCoreUtils::robustSummary)


# extract the aggregated protein-level data
msImpute.results <- assay(qf_imputed[["protein_imputed"]])

# statistical inference 
fit <- lmFit(msImpute.results, design)
msImpute.results <- list()

for (contrast in contrasts) {
  parts <- strsplit(contrast, " - ")[[1]]
  group1 <- trimws(parts[1])
  group2 <- trimws(parts[2])
  
  # Build contrast matrix using Condition prefix
  cont <- makeContrasts(
    contrasts = paste0("Condition", group1, " - Condition", group2),
    levels = design
  )
  
  # Fit contrasts
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- treat(fit2, lfc = log2(1.15), trend = TRUE)
  
  # Extract coef name (matches what treat() generates)
  coef_name <- paste0("Condition", group1, " - Condition", group2)
  
  # Get results
  result <- topTreat(fit2, coef = coef_name, adjust.method = "BH", 
                     number = Inf, confint = TRUE)
  
  # Store result
  msImpute.results[[contrast]] <- result
}


# 8. DreamAI+limma ----
DreamAI_imputed <- DreamAI(log2norm.df)
DreamAI_imputed <- DreamAI_imputed$Ensemble
# Remove rows with any NAs
DreamAI_imputed_clean <- DreamAI_imputed[complete.cases(DreamAI_imputed), ]

# Refit
fit <- lmFit(DreamAI_imputed_clean, design)

DreamAI.results <- list()

for (contrast in contrasts) {
  parts <- strsplit(contrast, " - ")[[1]]
  group1 <- trimws(parts[1])
  group2 <- trimws(parts[2])
  
  coef_name <- paste0("Condition", group1, " - Condition", group2)
  
  cont <- makeContrasts(
    contrasts = coef_name,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- treat(fit2, lfc = log2(1.15), trend = TRUE)
  
  result <- topTreat(fit2, coef = coef_name, adjust.method = "BH",
                     number = Inf, confint = TRUE)
  
  DreamAI.results[[contrast]] <- result
}


# 9. mice-pmm+limma ----
colnames(log2norm.df) <- metadata$rename
mice_imputed <- mice(log2norm.df, method = "pmm", seed = 456)
mice_pooled <- Reduce("+", lapply(1:5, function(i) complete(mice_imputed, i))) / 5
mice_imputed <- mice_pooled
# limma inference
fit <- lmFit(mice_imputed, design)


mice.results <- list()

for (contrast in contrasts) {
  parts <- strsplit(contrast, " - ")[[1]]
  group1 <- trimws(parts[1])
  group2 <- trimws(parts[2])
  
  coef_name <- paste0("Condition", group1, " - Condition", group2)
  
  cont <- makeContrasts(
    contrasts = coef_name,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- treat(fit2, lfc = log2(1.15), trend = TRUE)
  
  result <- topTreat(fit2, coef = coef_name, adjust.method = "BH",
                     number = Inf, confint = TRUE)
  
  mice.results[[contrast]] <- result
}


# 10. missBayes ----

group <- as.factor(make.names(metadata$Condition))

comparison <- makeContrasts(
  "X50amol - X500amol",
  "X100amol - X500amol",
  "X250amol - X500amol",
  "X1fmol - X500amol",
  "X5fmol - X500amol",
  levels = levels(group)
)

set.seed(746)
missBayes.results <- BayesMissingModel(values = log2norm.df, groups = group, comparisons = comparison, parallel = TRUE,
                             n.adapt = 500, burn.in = 500, n.iter = 5000, n.chains = 2, mcmcDiag = TRUE)


# 11. ND+limma ----
ND_imputed <- impute_normal(log2norm.df)

# limma inference
fit <- lmFit(ND_imputed, design)
ND.results <- list()

for (contrast in contrasts) {
  parts <- strsplit(contrast, " - ")[[1]]
  group1 <- trimws(parts[1])
  group2 <- trimws(parts[2])
  
  coef_name <- paste0("Condition", group1, " - Condition", group2)
  
  cont <- makeContrasts(
    contrasts = coef_name,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- treat(fit2, lfc = log2(1.15), trend = TRUE)
  
  result <- topTreat(fit2, coef = coef_name, adjust.method = "BH",
                     number = Inf, confint = TRUE)
  
  ND.results[[contrast]] <- result
}


# 12. limma ----

fit <- lmFit(log2norm.df, design)
limma.results <- list()

for (contrast in contrasts) {
  parts <- strsplit(contrast, " - ")[[1]]
  group1 <- trimws(parts[1])
  group2 <- trimws(parts[2])
  
  coef_name <- paste0("Condition", group1, " - Condition", group2)
  
  cont <- makeContrasts(
    contrasts = coef_name,
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, cont)
  fit2 <- treat(fit2, lfc = log2(1.15), trend = TRUE)
  
  result <- topTreat(fit2, coef = coef_name, adjust.method = "BH",
                     number = Inf, confint = TRUE)
  
  limma.results[[contrast]] <- result
}

# Keep proteins present in all methods for fair comparison
rownames(MSstats.results) <- MSstats.results$Protein
rownames(limpa_result_701020) <- limpa_result_701020$Protein.Names
rownames(proDA.results_701020) <- proDA.results_701020$name

common_proteins <- Reduce(intersect, list(rownames(MSstats.results_701020), rownames(msqRob2.results_701020), rownames(DEqMS.results_701020), rownames(proDA.results_701020), rownames(limpa_result_701020), 
                                          rownames(DreamAI.results_701020), rownames(mice.results_701020), rownames(msImpute.results_701020), rownames(ND.results_701020), rownames(missBayes.results_701020[[1]]), rownames(limma.results_701020)))
                          
# 13. Number of TP and FP ----
# set FDR to three levels : 0.01, 0.05 and 0.1
for (i in 1:5){
  missBayes.results[[i]] <- missBayes.results[[i]] %>% dplyr::rename(logFC = Median)
  MSstats.results[[i]] <- MSstats.results[[i]] %>% dplyr::rename(logFC = log2FC)
  proDA.results[[i]] <- proDA.results[[i]] %>% dplyr::rename(logFC = diff)
  rownames(MSstats.results[[i]]) <- MSstats.results[[i]]$Protein
  rownames(proDA.results[[i]]) <- proDA.results[[i]]$name
}
# go through 1 to 5 comparison. Note: change the logFC direction for comparison 4 and 5
method_list <- list(
  MSstats = MSstats.results[[2]],
  MSqRob2 = msqRob2.results[[2]],
  DEqMS = DEqMS.results[[2]],
  proDA = proDA.results[[2]],
  limpa = limpa.results[[2]],
  DreamAI_limma = DreamAI.results[[2]],
  mice_limma = mice.results[[2]],
  msImpute_limma = msImpute.results[[2]],
  ND_limma = ND.results[[2]],
  missBayes = missBayes.results[[2]],
  limma = limma.results[[2]]
)

results_list <- list()

FDR_thresholds <- c(0.01, 0.05, 0.1)
pval_cols <- c("adj.P.Val", "fisherAdjPval", "adj.pvalue", "adj_pval")

for (method_name in names(method_list)) {
  
  method_df <- method_list[[method_name]]
  metric_list <- vector("list", length(FDR_thresholds))
  
  is_ups <- grepl("ups", rownames(method_df))
  
  for (j in seq_along(FDR_thresholds)) {
    
    threshold <- FDR_thresholds[j]
    
    # Determine significance
    if (method_name == "missBayes") {
      significant <- (method_df$pLtROPE >= (100 - 100 * threshold)) |
        (method_df$pGtROPE >= (100 - 100 * threshold))
    } else {
      pcol <- pval_cols[pval_cols %in% colnames(method_df)][1]
      
      if (is.na(pcol)) {
        stop(sprintf("No adjusted p-value column found for %s", method_name))
      }
      
      significant <- method_df[[pcol]] <= threshold
    }
    
    # Direction check
    right_direction <- method_df$logFC > 0
    
    # Confusion matrix
    TP <- sum(is_ups & significant & right_direction, na.rm = TRUE)
    FP <- sum(!is_ups & significant, na.rm = TRUE)
    FN <- sum(is_ups & (!significant | !right_direction), na.rm = TRUE)
    TN <- sum(!is_ups & !significant, na.rm = TRUE)
    
    metric_list[[j]] <- data.frame(
      Threshold = threshold,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN
    )
  }
  
  results_list[[method_name]] <- do.call(rbind, metric_list)
}

# Combine all methods into one data frame
confusion_results <- do.call(
  rbind,
  lapply(names(results_list), function(method) {
    df <- results_list[[method]]
    df$Method <- method
    df
  })
)
rownames(confusion_results) <- NULL


plot_df <- confusion_results %>%
  mutate(
    Threshold = factor(
      Threshold,
      levels = c(0.01, 0.05, 0.10),
      labels = c("0.01", "0.05", "0.10")
    )
  )

ggplot(plot_df, aes(x = FP, y = TP)) +
  geom_point(
    aes(colour = Method, shape = Threshold),
    size = 3.8,
    alpha = 0.9
  ) +
  geom_text_repel(
    aes(label = Method, colour = Method),
    size = 3.2,
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.linewidth = 0.4
  ) +
  facet_wrap(~ Threshold, nrow = 1) +
  labs(
    x = "False Positives",
    y = "True Positives",
    colour = "Method",
    shape = "FDR Threshold"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1.3
    ),
    panel.grid.major = element_line(linewidth = 0.25),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.background = element_blank(),
    legend.key = element_blank()
  )

# 14. log2FC distributions ----
# Prepare UPS protein logFC data
ups_logFC_df <- bind_rows(
  lapply(names(method_list), function(method) {
    df <- method_list[[method]]
    
    data.frame(
      Method = method,
      Protein = rownames(df),
      logFC = df$logFC
    )
  })
) %>%
  filter(
    grepl("ups", Protein, ignore.case = TRUE),
    !is.na(logFC),
    is.finite(logFC)
  )

# Preserve original method order
ups_logFC_df$Method <- factor(
  ups_logFC_df$Method,
  levels = names(method_list)
)

# Compute medians
medians <- ups_logFC_df %>%
  group_by(Method) %>%
  summarise(
    median_logFC = median(logFC, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(ups_logFC_df, aes(x = Method, y = logFC, fill = Method)) +
  geom_boxplot(
    width = 0.85,
    alpha = 0.8,
    linewidth = 0.5,
    outlier.size = 0.7
  ) +
  geom_hline(
    yintercept = 3.3,
    colour = "red",
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_text(
    data = medians,
    aes(
      x = Method,
      y = 9.6,
      label = sprintf("%.2f", median_logFC)
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.8
  ) +
  coord_cartesian(ylim = c(0, 10)) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    x = "Methods",
    y = expression(logFC)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1.3
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 11
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 12
    ),
    axis.title.x = element_text(
      face = "bold",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      face = "bold",
      margin = margin(r = 10)
    ),
    legend.position = "none"
  )

# 15. Count detected number of ups/all proteins(coverage) ----
method_list <- list(
  MSstats = MSstats.results[[1]],
  MSqRob2 = msqRob2.results[[1]],
  DEqMS = DEqMS.results[[1]],
  proDA = proDA.results[[1]],
  limpa = limpa.results[[1]],
  DreamAI_limma = DreamAI.results[[1]],
  mice_limma = mice.results[[1]],
  msImpute_limma = msImpute.results[[1]],
  ND_limma = ND.results[[1]],
  missBayes = missBayes.results[[1]],
  limma = limma.results[[1]]
)

summary_detected <- bind_rows(lapply(names(method_list), function(method) {
  
  df <- method_list[[method]]
  
  # ensure logFC exists
  if (!"logFC" %in% colnames(df)) {
    stop(paste("logFC column missing in", method))
  }
  
  total_detected <- sum(!is.na(df$logFC) & is.finite(df$logFC))
  
  ups_detected <- sum(
    grepl("ups", rownames(df), ignore.case = TRUE) &
      !is.na(df$logFC) &
      is.finite(df$logFC)
  )
  
  data.frame(
    Method = method,
    Total_detected = total_detected,
    UPS_detected = ups_detected
  )
}))


