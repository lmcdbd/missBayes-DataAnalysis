library(data.table)
library(matrixStats)
library(tidyr)
library(ggplot2)
library(dplyr)
library(missBayes)
library(Matrix)

# Import data ----
pg<-fread("CQE/report.pg_matrix.tsv")
pr<-fread("CQE/report.pr_matrix.tsv")
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
# only keep proteins with missing values
log2.df <- log2all.df[rowSums(is.na(log2all.df)) >0,]

group <- as.factor(metadata$Condition)
comparison <- makeContrasts("G70_10_20 - G70_20_10", levels = levels(group))

# MCMC settings to test
mcmc_settings <- data.frame(
  No       = 1:7,
  n_adapt  = c(500, 500, 500, 500,  500, 1000,  500),
  burn_in  = c(500, 500, 500, 1000, 500,  500,  500),
  n_chains = c(2,   3,   3,   3,    2,    2,    4),
  n_iter   = c(5000,3333,5000,5000,10000,10000,10000)
)

n_repeats <- 3

# Results storage
results_list <- list()

for (s in 1:nrow(mcmc_settings)) {
  setting <- mcmc_settings[s, ]

  for (r in 1:n_repeats) {
    cat(sprintf("Setting %d, Repeat %d\n", s, r))

    set.seed(746 + r)  # Different seed per repeat but reproducible

    start_time <- Sys.time()

    result <- BayesMissingModel(
      values     = log2.df,
      groups     = group,
      comparisons = comparison,
      parallel   = TRUE,
      n.adapt    = setting$n_adapt,
      burn.in    = setting$burn_in,
      n.chains   = setting$n_chains,
      n.iter     = setting$n_iter,
      mcmcDiag   = TRUE
    )

    finish_time <- Sys.time()
    runtime_min <- as.numeric(difftime(finish_time, start_time, units = "mins"))

    # --- MCMC diagnostics ---
    median_ESS     <- median(result[[1]]$ESS, na.rm = TRUE)
    pct_ESS_gt1000 <- mean(result[[1]]$ESS > 1000, na.rm = TRUE) * 100
    median_MCSE    <- median(result[[1]]$MCSE, na.rm = TRUE)
    mean_rhat      <- mean(result[[1]]$max_rhat, na.rm = TRUE)
    pct_rhat_le101 <- mean(result[[1]]$max_rhat <= 1.01, na.rm = TRUE) * 100

    # --- Assign true values based on rownames ---
    rn <- rownames(result[[1]])
    true_vals <- ifelse(grepl("HUMAN", rn),  0,
                        ifelse(grepl("ECOLI", rn),  1,
                               ifelse(grepl("YEAST", rn), -1, NA)))

    # --- RMSE on Median column ---
    predicted <- result[[1]]$Median
    rmse <- sqrt(mean((predicted - true_vals)^2, na.rm = TRUE))

    # --- Accuracy and Precision ---
    # TP: ECOLI with pGtROPE >= 95, or YEAST with pLtROPE >= 95
    TP <- sum(
      (grepl("ECOLI", rn) & result[[1]]$pGtROPE >= 95) |
        (grepl("YEAST", rn) & result[[1]]$pLtROPE >= 95),
      na.rm = TRUE
    )

    # FP: HUMAN with either pGtROPE or pLtROPE >= 95
    FP <- sum(
      grepl("HUMAN", rn) & (result[[1]]$pGtROPE >= 95 | result[[1]]$pLtROPE >= 95),
      na.rm = TRUE
    )

    # FN: ECOLI or YEAST that are NOT TP
    FN <- sum(
      (grepl("ECOLI", rn) & result[[1]]$pGtROPE < 95) |
        (grepl("YEAST", rn) & result[[1]]$pLtROPE < 95),
      na.rm = TRUE
    )

    accuracy  <- TP / (TP + FN)
    precision <- TP / (TP + FP)

    # --- Store results ---
    results_list[[length(results_list) + 1]] <- data.frame(
      No              = setting$No,
      Repeat          = r,
      n_adapt         = setting$n_adapt,
      burn_in         = setting$burn_in,
      n_chains        = setting$n_chains,
      n_iter          = setting$n_iter,
      Median_ESS      = median_ESS,
      Pct_ESS_gt1000  = pct_ESS_gt1000,
      Median_MCSE     = median_MCSE,
      Mean_Rhat       = mean_rhat,
      Pct_Rhat_le1.01 = pct_rhat_le101,
      RMSE            = rmse,
      Accuracy        = accuracy,
      Precision       = precision,
      Runtime_min     = runtime_min
    )
  }
}

# Combine all results
results_df <- do.call(rbind, results_list)

# Summarise across repeats (mean over 3 runs per setting)
summary_df <- aggregate(. ~ No + n_adapt + burn_in + n_chains + n_iter,
                        data = results_df[, !names(results_df) %in% "Repeat"],
                        FUN = mean)
summary_df <- summary_df[order(summary_df$No), ]

