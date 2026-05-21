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
library(SummarizedExperiment)
library(Matrix)
library(iq)
library(purrr)
source("ND_imputation.R")

# Import data ----
pg<-fread("phagoFACS/report.pg_matrix.tsv")
pr<-fread("phagoFACS/report.pr_matrix.tsv")
raw <- fread("phagoFACS/report.tsv")
t <- table(unique(pr[,c('Protein.Group','Stripped.Sequence')])$Protein.Group)
pg$Peptide.Count <- t[match(pg$Protein.Group,names(t))]
metadata <- fread("phagoFACS/metadata.txt", header = TRUE)
#Filter by protein groups with <2 peptides, remove contaminants, and log2 transform the data
df.prot <- filter(pg, Peptide.Count > 1)
df.prot2 <- df.prot[!grepl("Cont_",df.prot$Protein.Group),]
df.prot.all<- df.prot2[, 5:40]
protein.matrix_all <- log2(as.matrix(df.prot.all))
log2all.df <- as.data.frame(protein.matrix_all)
rownames(log2all.df) <- make.unique(df.prot2$Protein.Names)
log2norm.df <- normalizeMedianValues(log2all.df)


# generate 20 splits
# Store all results and metadata splits
MSstats.results.list <- list()
metadata_split.list <- list()

# Generate all possible splits
library(combinat)

# Get replicate indices for one condition (same structure for all)
all_combinations <- combn(1:6, 3, simplify = FALSE)

# specify contrasts
contrasts <- c("PFA_4h_A - PFA_4h_B", "PFA_UP_A - PFA_UP_B", "PhoP_4h_A - PhoP_4h_B",
               "PhoP_UP_A - PhoP_UP_B", "WT_4H_A - WT_4H_B", "WT_UP_A - WT_UP_B")

# 1. MSstats ----
for (i in seq_along(all_combinations)) {
  i = 18
  cat("Running split", i, "of", length(all_combinations), "\n")

  # Split metadata using fixed combination instead of random
  metadata_split <- metadata %>%
    group_by(Condition) %>%
    mutate(
      idx = row_number(),
      split = ifelse(idx %in% all_combinations[[i]], "A", "B"),
      Condition = paste0(Condition, "_", split)
    ) %>%
    select(-idx, -split) %>%
    ungroup()

  metadata_split$BioReplicate <- 1:nrow(metadata_split)

  # Store metadata split
  metadata_split.list[[i]] <- metadata_split

  # Run MSstats
  tryCatch({
    MSstats_data <- DIANNtoMSstatsFormat(raw, metadata_split,
                                         removeFewMeasurements = FALSE,
                                         quantificationColumn = "FragmentQuantRaw")

    MSstats_results <- dataProcess(MSstats_data,
                                   normalization = "equalizeMedians",
                                   summaryMethod = "TMP",
                                   MBimpute = TRUE)

    levels <- levels(as.factor(MSstats_results$ProteinLevelData$GROUP))

    # Build contrast matrix
    comparison <- matrix(0, nrow = length(contrasts), ncol = length(levels))
    colnames(comparison) <- levels
    rownames(comparison) <- contrasts

    for (j in seq_along(contrasts)) {
      parts <- trimws(strsplit(contrasts[j], "-")[[1]])
      comparison[j, parts[1]] <- 1
      comparison[j, parts[2]] <- -1
    }

    # Run group comparison
    testResultAllComparisons <- groupComparison(contrast.matrix = comparison,
                                                data = MSstats_results)
    result <- testResultAllComparisons$ComparisonResult
    result <- result[!is.infinite(result$log2FC), ]

    MSstats.results.list[[i]] <- result

  }, error = function(e) {
    cat("Error in split", i, ":", conditionMessage(e), "\n")
    MSstats.results.list[[i]] <<- NULL
  })
  rm(MSstats_data, MSstats_results, result, comparison)
  gc()
}


# 2. MSqRob2 ----
raw[["File.Name"]] <- raw[["Run"]]

#Assumes report has already been filtered according to Q-values, removing contaminants, shared peptides, etc as required
raw_wide <- raw[,c(2:6, 14:16, 28)] %>% pivot_wider(names_from = Run, values_from = Precursor.Normalised, values_fn = sum)

msqRob2.results.list <- list()

for (i in seq_along(metadata_split.list)) {
  cat("Running msqrob2 split", i, "of", length(metadata_split.list), "\n")

  tryCatch({
    # Prepare metadata
    meta <- as.data.frame(metadata_split.list[[i]])
    meta[["Run"]] <- sort(unique(raw[["File.Name"]]))
    colnames(meta)[1] <- "quantCols"
    rownames(meta) <- meta$quantCols

    # Build QFeatures object
    qf <- readQFeatures(raw_wide, colData = meta,
                        quantCols = meta$quantCols, name = "quant")
    qf <- logTransform(qf, base = 2, i = "quant", name = "log2quant")

    # Align metadata
    metadata_aligned <- meta[rownames(colData(qf)), ]
    colData(qf)$Condition <- as.factor(metadata_aligned$Condition)

    # Normalization
    qf <- normalize(qf, i = "log2quant", name = "log2norm", method = "center.median")

    # Aggregate and fit hurdle model
    qf_agg <- aggregateFeatures(qf, i = "log2norm", fcol = "Protein.Names",
                                name = "log2proteinQuant")
    qf_hurdle <- msqrobHurdle(qf_agg, i = "log2proteinQuant", formula = ~Condition)

    # Get coefficients
    coef <- getCoef(rowData(qf_hurdle[["log2proteinQuant"]])$msqrobHurdleIntensity[[1]])
    ref_level <- levels(as.factor(metadata_aligned$Condition))[1]  # reference level

    # Run all contrasts
    split_results <- list()

    for (contrast in contrasts) {
      parts <- strsplit(contrast, " - ")[[1]]
      group1 <- trimws(parts[1])
      group2 <- trimws(parts[2])

      # Handle reference level
      if (group1 == ref_level) {
        contrast_str <- paste0("-Condition", group2, " = 0")
        contrast_name <- paste0("hurdle_-Condition", group2)
      } else if (group2 == ref_level) {
        contrast_str <- paste0("Condition", group1, " = 0")
        contrast_name <- paste0("hurdle_Condition", group1)
      } else {
        contrast_str <- paste0("Condition", group1, " - Condition", group2, " = 0")
        contrast_name <- paste0("hurdle_Condition", group1, " - Condition", group2)
      }

      contrast_matrix <- makeContrast(contrasts = contrast_str,
                                      parameterNames = names(coef))

      msqrob2_result <- hypothesisTestHurdle(object = qf_hurdle,
                                             contrast = contrast_matrix,
                                             i = "log2proteinQuant")

      split_results[[contrast]] <- rowData(msqrob2_result[["log2proteinQuant"]])[[contrast_name]]
    }

    msqRob2.results.list[[i]] <- split_results

  }, error = function(e) {
    cat("Error in split", i, ":", conditionMessage(e), "\n")
    msqRob2.results.list[[i]] <<- NULL
  })
  rm(qf, qf_agg, qf_hurdle, msqrob2_result, contrast_matrix)
  gc()
}

# 3. DEqMS ----
run_DEqMS_split <- function(meta, log2norm.df, pep.count.table, contrasts) {
  
  # Build design matrix
  class <- as.factor(meta$Condition)
  design <- model.matrix(~0 + class)
  colnames(design) <- sub("^class", "", colnames(design))
  
  # Fit model once per split
  fit1 <- lmFit(log2norm.df, design = design)
  
  split_results <- list()
  
  for (contrast in contrasts) {
    parts <- trimws(strsplit(contrast, " - ")[[1]])
    group1 <- parts[1]
    group2 <- parts[2]
    
    coef_name <- paste0(group1, " - ", group2)
    
    cont <- makeContrasts(
      contrasts = coef_name,
      levels = colnames(design)
    )
    
    fit2 <- contrasts.fit(fit1, contrasts = cont)
    fit3 <- eBayes(fit2)
    
    # Add peptide counts
    fit3$count <- pep.count.table[rownames(fit3$coefficients), "count"]
    
    if (is.na(min(fit3$count)) || min(fit3$count) == 0) {
      warning(paste("Count issue in contrast:", contrast))
      next
    }
    
    fit4 <- spectraCounteBayes(fit3)
    split_results[[contrast]] <- outputResult(fit4, coef_col = 1)
  }
  
  return(split_results)
}

# Compute pep.count.table once outside loop
pep.count.table <- data.frame(
  count = as.vector(pg$Peptide.Count),
  row.names = make.unique(pg$Protein.Names)
)

# Main loop
DEqMS.results.list <- list()

for (i in seq_along(metadata_split.list)) {
  cat("Running DEqMS split", i, "of", length(metadata_split.list), "\n")
  
  tryCatch({
    DEqMS.results.list[[i]] <- run_DEqMS_split(
      meta = as.data.frame(metadata_split.list[[i]]),
      log2norm.df = log2norm.df,
      pep.count.table = pep.count.table,
      contrasts = contrasts
    )
  }, error = function(e) {
    cat("Error in split", i, ":", conditionMessage(e), "\n")
    DEqMS.results.list[[i]] <<- NULL
  })
  
  gc()
  saveRDS(DEqMS.results.list, "DEqMS_results_checkpoint.rds")
}

# 4. proDA ----
proDA.results.list <- list()

run_proDA_split <- function(meta, log2norm.df, contrasts) {

  # Build design matrix
  design <- model.matrix(~0 + Condition, data = meta)
  colnames(design) <- sub("^Condition", "", colnames(design))

  # Fit dropout model
  fit <- proDA(log2norm.df, design = design)

  split_results <- list()

  for (contrast in contrasts) {
    parts <- strsplit(contrast, " - ")[[1]]
    group1 <- trimws(parts[1])
    group2 <- trimws(parts[2])

    contrast_expr <- paste0("`", group1, "` - `", group2, "`")

    split_results[[contrast]] <- test_diff(fit,
                                           contrast = eval(parse(text = contrast_expr)))
  }

  return(split_results)
}

for (i in seq_along(metadata_split.list)) {
  cat("Running proDA split", i, "of", length(metadata_split.list), "\n")

  tryCatch({
    proDA.results.list[[i]] <- run_proDA_split(
      meta = as.data.frame(metadata_split.list[[i]]),
      log2norm.df = log2norm.df,
      contrasts = contrasts
    )
  }, error = function(e) {
    cat("Error in split", i, ":", conditionMessage(e), "\n")
    proDA.results.list[[i]] <<- NULL
  })

  gc()
}



# 5. limpa ----
run_limpa_split <- function(meta, y, contrasts) {

  # Match metadata to column order in y$E
  metadata_ordered <- meta[match(colnames(y$E), meta$Run), ]
  stopifnot(all(metadata_ordered$Run == colnames(y$E)))

  # Create design matrix
  group <- as.factor(metadata_ordered$Condition)
  design <- model.matrix(~0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  rownames(design) <- colnames(y$E)

  # Fit model (once per split)
  fit <- dpcDE(y, design, plot = FALSE)
  fit <- eBayes(fit)

  split_results <- list()

  for (contrast in contrasts) {
    parts <- trimws(strsplit(contrast, "-")[[1]])
    group_A <- parts[1]
    group_B <- parts[2]

    contrast_matrix <- makeContrasts(
      contrasts = paste0(group_A, " - ", group_B),
      levels = colnames(design)
    )

    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)

    split_results[[contrast]] <- topTable(fit2, coef = 1,
                                          number = Inf,
                                          adjust.method = "BH")
  }

  return(split_results)
}

y.prec <- readDIANN(file = "P:/Trost-group/Mengchun/8) MIP/Data/real proteomics dataset/report.tsv",
                    run.column = "Run", q.columns = c("Q.Value", "Protein.Q.Value"))
dpcest <- dpcCN(y.prec) # estimate dpc
plotDPC(dpcest)
y <- dpcQuant(y.prec, dpc = dpcest) # summarization
# Main loop
limpa.results.list <- list()

for (i in seq_along(metadata_split.list)) {
  cat("Running limpa split", i, "of", length(metadata_split.list), "\n")

  tryCatch({
    limpa.results.list[[i]] <- run_limpa_split(
      meta = as.data.frame(metadata_split.list[[i]]),
      y = y,
      contrasts = contrasts
    )
  }, error = function(e) {
    cat("Error in split", i, ":", conditionMessage(e), "\n")
    limpa.results.list[[i]] <<- NULL
  })

  gc()
}


# 7. msImpute+limma ----
# Convert to log2
raw <- fread("phagoFACS/report.tsv")
pep_log2 <- as.matrix(log2(pr[, 11:46]))

pep_log2[is.infinite(pep_log2)] <- NA
# median normalization
pep_log2 <- normalizeMedianValues(pep_log2)
rownames(pep_log2) <- pr$Stripped.Sequence

run_msImpute_split <- function(meta, pep_log2, raw, pep_to_prot, contrasts) {

  # Check what pep_log2 colnames look like vs meta$Run
  cat("pep_log2 colnames example:", head(colnames(pep_log2)), "\n")
  cat("meta$Run example:", head(meta$Run), "\n")

  # Reorder pep_log2 columns to match meta$Run order BEFORE imputation
  # Try direct match first
  if (all(meta$Run %in% colnames(pep_log2))) {
    pep_log2_ordered <- pep_log2[, meta$Run]
  } else {
    # Try matching by basename (strip path/extension differences)
    col_base <- tools::file_path_sans_ext(basename(colnames(pep_log2)))
    run_base <- tools::file_path_sans_ext(basename(meta$Run))
    idx <- match(run_base, col_base)
    cat("Match via basename:", sum(!is.na(idx)), "of", length(run_base), "\n")
    pep_log2_ordered <- pep_log2[, idx]
    colnames(pep_log2_ordered) <- meta$Run  # rename to match meta$Run exactly
  }

  # msImpute imputation on reordered matrix
  design_imp <- model.matrix(~0 + Condition, data = meta)
  msImpute.pep <- msImpute(pep_log2_ordered, method = "v2-mnar",
                           design = design_imp, relax_min_obs = TRUE)

  # Now colnames should match meta$Run
  cat("Columns match after reorder:", identical(colnames(msImpute.pep), meta$Run), "\n")

  # Create SummarizedExperiment
  se_pep <- SummarizedExperiment(
    assays = list(log2quant = msImpute.pep),
    colData = DataFrame(meta, row.names = meta$Run)
  )
  rownames(se_pep) <- make.unique(rownames(se_pep))

  # rest of function unchanged ...

  # Wrap in QFeatures
  qf_imputed <- QFeatures(list(pep_imputed = se_pep))

  # Add peptide-to-protein mapping
  rowData(qf_imputed[["pep_imputed"]]) <- pep_to_prot[rownames(msImpute.pep), ]

  # Assign condition from metadata
  colData(qf_imputed)$Condition <- as.factor(
    meta$Condition[match(rownames(colData(qf_imputed)), meta$Run)]
  )

  # Fix colData for aggregation
  cd <- colData(qf_imputed)
  cd$Condition <- as.character(cd$Condition)
  colData(qf_imputed) <- cd
  colData(qf_imputed[["pep_imputed"]])$Condition <- cd$Condition

  # Aggregate to protein level
  qf_imputed <- aggregateFeatures(qf_imputed,
                                  i = "pep_imputed",
                                  fcol = "Protein.Group",
                                  name = "protein_imputed",
                                  fun = MsCoreUtils::robustSummary)

  # Extract protein-level matrix
  msImpute.mat <- assay(qf_imputed[["protein_imputed"]])

  # Build conditions using Run to match original column order
  conditions <- as.factor(meta$Condition[match(colnames(msImpute.mat), meta$Run)])
  cat("Conditions found:", levels(conditions), "\n")
  cat("Any NA conditions:", any(is.na(conditions)), "\n")

  design_lm <- model.matrix(~0 + conditions)
  colnames(design_lm) <- gsub("conditions", "", colnames(design_lm))

  fit <- lmFit(msImpute.mat, design_lm)

  split_results <- list()

  for (contrast in contrasts) {
    parts <- trimws(strsplit(contrast, " - ")[[1]])
    group1 <- parts[1]
    group2 <- parts[2]

    coef_name <- paste0(group1, " - ", group2)

    cont <- makeContrasts(
      contrasts = coef_name,
      levels = colnames(design_lm)
    )

    fit2 <- contrasts.fit(fit, cont)
    fit2 <- treat(fit2, lfc = log2(1.15), trend = TRUE)

    result <- topTreat(fit2, coef = coef_name, adjust.method = "BH",
                       number = Inf, confint = TRUE)

    # Replace rownames with protein names
    new_rownames <- raw$Protein.Names[match(rownames(result), raw$Protein.Group)]
    rownames(result) <- make.unique(as.character(new_rownames))

    split_results[[contrast]] <- result
  }

  return(split_results)
}


pep_to_prot <- unique(raw[, .(Stripped.Sequence, Protein.Group)])
setDF(pep_to_prot)
rownames(pep_to_prot) <- pep_to_prot$Stripped.Sequence

# Main loop
msImpute.results.list <- list()

for (i in seq_along(metadata_split.list)) {
  cat("Running msImpute split", i, "of", length(metadata_split.list), "\n")

  tryCatch({
    msImpute.results.list[[i]] <- run_msImpute_split(
      meta = as.data.frame(metadata_split.list[[i]]),
      pep_log2 = pep_log2,
      raw = raw,
      pep_to_prot = pep_to_prot,
      contrasts = contrasts
    )
  }, error = function(e) {
    cat("Error in split", i, ":", conditionMessage(e), "\n")
    msImpute.results.list[[i]] <<- NULL
  })

  gc()
}

# 8. DreamAI+limma ----
DreamAI_imputed <- DreamAI(log2norm.df)
DreamAI_imputed <- DreamAI_imputed$Ensemble
run_DreamAI_split <- function(meta, DreamAI_imputed, contrasts) {
  
  # Remove rows that are entirely or mostly NA
  na_count <- rowSums(is.na(DreamAI_imputed))
  DreamAI_clean <- DreamAI_imputed[na_count == 0, ]
  cat("Removed", nrow(DreamAI_imputed) - nrow(DreamAI_clean), "rows with NAs\n")
  
  # Build design matrix
  design_lm <- model.matrix(~0 + Condition, data = meta)
  colnames(design_lm) <- gsub("Condition", "", colnames(design_lm))
  
  # Fit model once per split
  fit <- lmFit(DreamAI_clean, design_lm)
  
  split_results <- list()
  
  for (contrast in contrasts) {
    parts <- trimws(strsplit(contrast, " - ")[[1]])
    group1 <- parts[1]
    group2 <- parts[2]
    
    coef_name <- paste0(group1, " - ", group2)
    
    cont <- makeContrasts(
      contrasts = coef_name,
      levels = colnames(design_lm)
    )
    
    fit2 <- contrasts.fit(fit, cont)
    fit2 <- treat(fit2, lfc = log2(1.15), trend = TRUE)
    
    result <- topTreat(fit2, coef = coef_name, adjust.method = "BH",
                       number = Inf, confint = TRUE)
    
    split_results[[contrast]] <- result
  }
  
  return(split_results)
}
# Main loop
DreamAI.results.list <- list()

for (i in seq_along(metadata_split.list)) {
  cat("Running DreamAI+limma split", i, "of", length(metadata_split.list), "\n")
  
  tryCatch({
    DreamAI.results.list[[i]] <- run_DreamAI_split(
      meta = as.data.frame(metadata_split.list[[i]]),
      DreamAI_imputed = DreamAI_imputed,
      contrasts = contrasts
    )
  }, error = function(e) {
    cat("Error in split", i, ":", conditionMessage(e), "\n")
    DreamAI.results.list[[i]] <<- NULL
  })
  
  gc()
  saveRDS(DreamAI.results.list, "DreamAI_results_checkpoint.rds")
}


# 9. mice-pmm+limma ----
run_mice_split <- function(meta, log2norm.df, contrasts) {
  
  # Assign column names from metadata
  colnames(log2norm.df) <- meta$rename
  
  # MICE imputation
  mice_imputed <- mice(log2norm.df, method = "pmm", seed = 456)
  mice_pooled <- Reduce("+", lapply(1:5, function(i) complete(mice_imputed, i))) / 5
  mice_imputed <- mice_pooled
  
  # Build design matrix
  class <- as.factor(meta$Condition)
  design_lm <- model.matrix(~0 + class)
  colnames(design_lm) <- sub("^class", "", colnames(design_lm))
  
  # Fit model once per split
  fit <- lmFit(mice_imputed, design_lm)
  
  split_results <- list()
  
  for (contrast in contrasts) {
    parts <- trimws(strsplit(contrast, " - ")[[1]])
    group1 <- parts[1]
    group2 <- parts[2]
    
    coef_name <- paste0(group1, " - ", group2)
    
    cont <- makeContrasts(
      contrasts = coef_name,
      levels = colnames(design_lm)
    )
    
    fit2 <- contrasts.fit(fit, cont)
    fit2 <- treat(fit2, lfc = log2(1.15), trend = TRUE)
    
    result <- topTreat(fit2, coef = coef_name, adjust.method = "BH",
                       number = Inf, confint = TRUE)
    
    split_results[[contrast]] <- result
  }
  
  return(split_results)
}

# Main loop
mice.results.list <- list()

for (i in seq_along(metadata_split.list)) {
  cat("Running MICE split", i, "of", length(metadata_split.list), "\n")
  
  tryCatch({
    mice.results.list[[i]] <- run_mice_split(
      meta = as.data.frame(metadata_split.list[[i]]),
      log2norm.df = log2norm.df,
      contrasts = contrasts
    )
  }, error = function(e) {
    cat("Error in split", i, ":", conditionMessage(e), "\n")
    mice.results.list[[i]] <<- NULL
  })
  
  gc()
  saveRDS(mice.results.list, "mice_results_checkpoint.rds")
}


# 10. missBayes ----
run_missBayes_split <- function(meta, log2norm.df, contrasts) {

  # Build group factor
  group <- as.factor(meta$Condition)

  # Build comparison matrix
  comparison <- makeContrasts(
    contrasts = contrasts,
    levels = levels(group)
  )

  # Run BayesMissingModel
  set.seed(746)
  result <- BayesMissingModel(
    values = log2norm.df,
    groups = group,
    comparisons = comparison,
    parallel = TRUE,
    n.adapt = 500,
    burn.in = 500,
    n.iter = 10000,
    n.chains = 2,
    mcmcDiag = TRUE
  )

  return(result)
}

# Main loop
missBayes.results.list <- list()

for (i in seq_along(metadata_split.list)) {
  cat("Running BayesMissingModel split", i, "of", length(metadata_split.list), "\n")

  tryCatch({
    missBayes.results.list[[i]] <- run_missBayes_split(
      meta = as.data.frame(metadata_split.list[[i]]),
      log2norm.df = log2norm.df,
      contrasts = contrasts
    )
  }, error = function(e) {
    cat("Error in split", i, ":", conditionMessage(e), "\n")
    missBayes.results.list[[i]] <<- NULL
  })

  gc()
}


# 11. ND+limma ----
run_ND_split <- function(meta, log2norm.df, contrasts) {
  
  # Normal distribution imputation
  ND_imputed <- impute_normal(log2norm.df)
  
  # Build design matrix
  design_lm <- model.matrix(~0 + Condition, data = meta)
  colnames(design_lm) <- gsub("Condition", "", colnames(design_lm))
  
  # Fit model
  fit <- lmFit(ND_imputed, design_lm)
  
  split_results <- list()
  
  for (contrast in contrasts) {
    parts <- trimws(strsplit(contrast, " - ")[[1]])
    group1 <- parts[1]
    group2 <- parts[2]
    
    coef_name <- paste0(group1, " - ", group2)
    
    cont <- makeContrasts(
      contrasts = coef_name,
      levels = colnames(design_lm)
    )
    
    fit2 <- contrasts.fit(fit, cont)
    fit2 <- treat(fit2, lfc = log2(1.15), trend = TRUE)
    
    result <- topTreat(fit2, coef = coef_name, adjust.method = "BH",
                       number = Inf, confint = TRUE)
    
    split_results[[contrast]] <- result
  }
  
  return(split_results)
}

# Main loop
ND.results.list <- list()

for (i in seq_along(metadata_split.list)) {
  cat("Running ND+limma split", i, "of", length(metadata_split.list), "\n")
  
  tryCatch({
    ND.results.list[[i]] <- run_ND_split(
      meta = as.data.frame(metadata_split.list[[i]]),
      log2norm.df = log2norm.df,
      contrasts = contrasts
    )
  }, error = function(e) {
    cat("Error in split", i, ":", conditionMessage(e), "\n")
    ND.results.list[[i]] <<- NULL
  })
  
  gc()
  saveRDS(ND.results.list, "ND_results_checkpoint.rds")
}



# 12. limma ----
run_limma_split <- function(meta, log2norm.df, contrasts) {
  
  # Build design matrix
  design_lm <- model.matrix(~0 + Condition, data = meta)
  colnames(design_lm) <- gsub("Condition", "", colnames(design_lm))
  
  # Fit model once per split
  fit <- lmFit(log2norm.df, design_lm)
  
  split_results <- list()
  
  for (contrast in contrasts) {
    parts <- trimws(strsplit(contrast, " - ")[[1]])
    group1 <- parts[1]
    group2 <- parts[2]
    
    coef_name <- paste0(group1, " - ", group2)
    
    cont <- makeContrasts(
      contrasts = coef_name,
      levels = colnames(design_lm)
    )
    
    fit2 <- contrasts.fit(fit, cont)
    fit2 <- treat(fit2, lfc = log2(1.15), trend = TRUE)
    
    result <- topTreat(fit2, coef = coef_name, adjust.method = "BH",
                       number = Inf, confint = TRUE)
    
    split_results[[contrast]] <- result
  }
  
  return(split_results)
}

# Main loop
limma.results.list <- list()

for (i in seq_along(metadata_split.list)) {
  cat("Running limma split", i, "of", length(metadata_split.list), "\n")
  
  tryCatch({
    limma.results.list[[i]] <- run_limma_split(
      meta = as.data.frame(metadata_split.list[[i]]),
      log2norm.df = log2norm.df,
      contrasts = contrasts
    )
  }, error = function(e) {
    cat("Error in split", i, ":", conditionMessage(e), "\n")
    limma.results.list[[i]] <<- NULL
  })
  
  gc()
  saveRDS(limma.results.list, "limma_results_checkpoint.rds")
}

# 13. distributions of number of discoveries per split for each method ----
# Count discoveries in every dataframe
discovery_counts_proda <- map_dfr(
  seq_along(proDA.results.list),
  function(i) {
    
    map_dfr(
      seq_along(proDA.results.list[[i]]),
      function(j) {
        
        df <- proDA.results.list[[i]][[j]]
        
        # proteins with valid test results
        valid_proteins <- sum(
          !is.na(df$adj_pval) 
        )
        
        # discoveries
        discoveries <- sum(
          df$adj_pval <= 0.05,
          na.rm = TRUE
        )
        
        tibble(
          outer_list = i,
          inner_list = j,
          n_discoveries = discoveries,
          n_tested = valid_proteins,
          prop_discoveries = discoveries / valid_proteins
        )
      }
    )
  }
)

# MSstats is different 
discovery_counts_MSstats <- map_dfr(
  seq_along(MSstats.results.list),
  function(i) {
    
    df <- MSstats.results.list[[i]]
    
    map_dfr(
      unique(df$Label),
      function(label) {
        
        df_sub <- df[df$Label == label, ]
        
        # proteins with valid test results
        valid_proteins <- sum(!is.na(df_sub$adj.pvalue))
        
        # discoveries
        discoveries <- sum(df_sub$adj.pvalue <= 0.05, na.rm = TRUE)
        
        tibble(
          split = i,
          contrast = label,
          n_discoveries = discoveries,
          n_tested = valid_proteins,
          prop_discoveries = discoveries / valid_proteins
        )
      }
    )
  }
)



table(discovery_counts_missbayes$n_discoveries)

ggplot(discovery_counts_MSstats,
       aes(x = factor(n_discoveries))) +
  geom_bar(
    width = 0.8,
    fill = "steelblue",
    color = "black"
  ) +
  labs(
    title = NULL,
    subtitle = "MSstats",
    x = "Number of discoveries",
    y = "Frequency"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      linewidth = 0.8
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    axis.title = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 14,
      face = "italic",
      hjust = 0.5
    )
  )


ggplot(discovery_counts_missbayes,
       aes(x = n_discoveries)) +
  geom_histogram(
    binwidth = 10,
    fill = "steelblue",
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq(0, 1000, by = 100)
  ) +
  labs(
    title = NULL,
    subtitle = "MissBayes (threshold = 0.05)",
    x = "Number of discoveries",
    y = "Frequency"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      linewidth = 0.8
    ),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 14,
      face = "italic",
      hjust = 0.5
    )
  )
# 14. proportion of splits with at least one discovery ----
sum(discovery_counts_MSstats$n_discoveries > 0) / nrow(discovery_counts_MSstats)

# 15. proportion of proteins called significant ----

# 16. distribution of estimated logFCs around zero ----
# extract ALL values from column "Median"
all_FC <- map_dfr(
  seq_along(proDA.results.list),
  function(i) {
    
    map_dfr(
      seq_along(proDA.results.list[[i]]),
      function(j) {
        
        df <- proDA.results.list[[i]][[j]]
        
        tibble(
          outer_list = i,
          inner_list = j,
          Median = df$diff
        )
      }
    )
  }
)

all_FC_MSstats <- map_dfr(
  seq_along(MSstats.results.list),
  function(i) {
    
    df <- MSstats.results.list[[i]]
    
    map_dfr(
      unique(df$Label),
      function(label) {
        
        df_sub <- df[df$Label == label, ]
        
        tibble(
          outer_list = i,
          contrast   = label,
          Median     = df_sub$log2FC
        )
      }
    )
  }
)

# Calculate mean and SD
fc_stats <- all_FC_MSstats %>%
  summarise(
    mean_fc = mean(Median, na.rm = TRUE),
    sd_fc   = sd(Median, na.rm = TRUE)
  )

# Create label string
stats_label <- paste0(
  "Mean = ", round(fc_stats$mean_fc, 3), "\n",
  "SD = ",   round(fc_stats$sd_fc, 3)
)

# Plot
ggplot(all_FC_MSstats, aes(x = Median)) +
  geom_histogram(
    bins = 50,
    fill = "steelblue",
    color = "black"
  ) +
  annotate("text",
           x = 1.5, y = Inf,       # top right corner
           label = stats_label,
           hjust = 1, vjust = 1.5,
           fontface = "bold",
           size = 4.5
  ) +
  xlim(-2, 2) +
  labs(
    title    = NULL,
    subtitle = "MSstats",
    x        = "logFC",
    y        = "Frequency"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.title   = element_text(face = "bold"),
    axis.text    = element_text(color = "black"),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5)
  )

# 17. distribution of discovery rate ----
# Combine all methods
method_list <- list(
  MSstats         = discovery_counts_MSstats,
  MSqRob2         = discovery_counts_msqrob2,
  DEqMS           = discovery_counts_deqms,
  proDA           = discovery_counts_proda,
  limpa           = discovery_counts_limpa,
  DreamAI_limma   = discovery_counts_DreamAI,
  mice_limma      = discovery_counts_mice,
  msImpute_limma  = discovery_counts_msImpute,
  ND_limma        = discovery_counts_ND,
  `missBayes (0.05)` = discovery_counts_missbayes,
  `missBayes (0.01)` = discovery_counts_missbayes_0.01
)

all_results <- bind_rows(
  lapply(names(method_list), function(m) {
    method_list[[m]] %>% mutate(method = m)
  })
)

# Fix factor order
all_results$method <- factor(all_results$method, levels = names(method_list))

# Compute medians
medians <- all_results %>%
  group_by(method) %>%
  summarise(median_fdr = median(prop_discoveries, na.rm = TRUE))

set3_colors <- RColorBrewer::brewer.pal(11, "Set3")
custom_colors <- set3_colors
names(custom_colors) <- levels(all_results$method)

custom_colors["missBayes (0.01)"] <- "#E41A1C"  

ggplot(all_results, aes(x = method, y = prop_discoveries, fill = method)) +
  geom_boxplot(width = 0.9, alpha = 0.8) +
  geom_text(
    data = medians,
    aes(x = method, label = round(median_fdr, 3)),
    y = 0.19,          # just inside the new ylim upper bound
    fontface = "bold",
    size = 4
  ) +
  geom_hline(yintercept = 0.05, color = "red", linewidth = 1, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 0.2)) +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Methods", y = "Proportion of discoveries") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y  = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold",
                                margin = ggplot2::margin(t = 10)),
    axis.title.y = element_text(size = 14, face = "bold",
                                margin = ggplot2::margin(r = 10)),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )
