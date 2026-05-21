library(data.table)
library(matrixStats)
library(tidyr)
library(ggplot2)
library(dplyr)

pg<-fread("phagoFACS/report.pg_matrix.tsv")
pr<-fread("phagoFACS/report.pr_matrix.tsv")
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
log2norm.df_R <- normalizeMedianValues(log2all.df)

r <- 6

for (x in 11:17) {
  
  # create file name for each plot
  file_name <- paste0(output_dir, "missingpattern_", x , "_2", ".png")
  
  # open high-res PNG device
  png(file_name, width = 4, height = 4, units = "in", res = 300)
  
  row_mean_x <- list()
  
  for (i in 1:(ncol(log2norm.df_R) %/% r)) {
    group_i <- log2norm.df_R[, (r * (i - 1) + 1):(r * i)]
    mean_is_x <- rowMeans(group_i, na.rm = TRUE) <= x  & rowMeans(group_i, na.rm = TRUE) >= x - 0.5
    group_i_x <- group_i[mean_is_x, ]
    nona_counts_x_i <- apply(group_i_x, 1, function(row) sum(!is.na(row)))
    row_mean_x[[i]] <- nona_counts_x_i
  }
  
  # combine
  row_mean_x_vector <- unname(unlist(row_mean_x))
  row_mean_x_vector <- row_mean_x_vector[row_mean_x_vector != 0]
  freq_table_x <- table(row_mean_x_vector)
  
  # choose ylim only for the second plot 
  if (x == 11) {
    ylim_vals <- c(0, 200)   
  } else if (x == 12) {
    ylim_vals <- c(0, 950)
  } else {
    ylim_vals <- NULL       
  }
  
  
  # barplot
  bp <- barplot(freq_table_x,
                main = paste("Group mean in [", x - 0.5, ",", x , "]"),
                xlab = "Observed Number",
                ylab = "Frequency",
                col = "skyblue",
                ylim = ylim_vals)
  
  # ---- FIT p via MLE ----
  n <- r
  data <- as.numeric(row_mean_x_vector)
  
  fit <- optimize(function(p) -loglik_zbinom(p, data, n),
                  interval = c(1e-6, 1-1e-6))
  p_hat <- fit$minimum
  
  # fitted distribution
  k_vals <- as.numeric(names(freq_table_x))
  probs <- dzbinom(k_vals, n, p_hat)
  fitted_counts <- probs * sum(freq_table_x)
  
  # overlay
  points(bp, fitted_counts, pch = 19, col = "red", cex = 1.2)
  lines(bp, fitted_counts, col = "red", lwd = 2)
  
  # show fitted p in subtitle
  mtext(sprintf("p = %.3f", p_hat), side = 3, line = -1, cex = 0.8)
  
  dev.off()
  
  cat("Saved plot:", file_name, "\n")
}

