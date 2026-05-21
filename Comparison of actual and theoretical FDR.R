# 1. import CQE dataset, use 70_10_20 comparison
pg<-fread("CQE/report.pg_matrix.tsv")
pr<-fread("CQE/report.pr_matrix.tsv")
t <- table(unique(pr[,c('Protein.Group','Stripped.Sequence')])$Protein.Group)
pg$Peptide.Count <- t[match(pg$Protein.Group,names(t))]
metadata <- fread("CQE/metadata.txt", header = TRUE)
#Filter by protein groups with <2 peptides
df.prot <- filter(pg, Peptide.Count > 1)
df.prot2 <- df.prot[!grepl("Cont_",df.prot$Protein.Group),]
df.prot.all<- df.prot2[, 5:29]
protein.matrix_all <- log2(as.matrix(df.prot.all))
log2all.df <- as.data.frame(protein.matrix_all)
rownames(log2all.df) <- df.prot2$Protein.Names
# 2. ROPE 0.1, 0.2, 0.3, 0.4, 0.5
groups <- as.factor(metadata$Condition)
comparison <- makeContrasts("G70_10_20 - G70_20_10", levels = levels(groups))
output_01 <- BayesMissingModel(log2all.df, groups, comparison, ROPE = c(-0.1, 0.1))
output_02 <- BayesMissingModel(log2all.df, groups, comparison, ROPE = c(-0.2, 0.2))
output_03 <- BayesMissingModel(log2all.df, groups, comparison, ROPE = c(-0.3, 0.3))
output_04 <- BayesMissingModel(log2all.df, groups, comparison, ROPE = c(-0.4, 0.4))
output_05 <- BayesMissingModel(log2all.df, groups, comparison, ROPE = c(-0.5, 0.5))

# 3. calculate theoretical FDR

TP <- 0
TN <- 0
FP <- 0
FN <- 0
for (i in 1:nrow(output_05[[1]])) {
  row_name <- rownames(output_05[[1]])[i]
  p_val <- output_05[[1]]$pLtCompVal[i]
  fc <- output_05[[1]]$Median[i]
  pLtROPE <- output_05[[1]]$pLtROPE[i]
  pGtROPE <- output_05[[1]]$pGtROPE[i]

  # Skip iteration if p_val or fc is NA
  if (is.na(p_val) || is.na(fc)) next

  if (grepl("_HUMAN", row_name) && pLtROPE < 90 && pGtROPE < 90) {
    TN <- TN + 1
  } else if (grepl("_YEAST", row_name) && pLtROPE > 90 ) {
    TP <- TP + 1
  } else if (grepl("_ECOLI", row_name) && pGtROPE > 90 ) {
    TP <- TP + 1
  } else if ((grepl("_HUMAN", row_name) && pLtROPE >= 90) | (grepl("_HUMAN", row_name) && pGtROPE >= 90) |
             (grepl("_YEAST", row_name) && pGtROPE >= 90) | (grepl("_ECOLI", row_name) && pLtROPE >= 90)) {
    FP <- FP + 1
  } else if (grepl("_ECOLI", row_name)|grepl("_YEAST", row_name)){
    FN <- FN + 1
  }
}

cat("True Positives:", TP, "\n")
cat("True Negatives:", TN, "\n")
cat("False Positives:", FP, "\n")
cat("False Negatives:", FN, "\n")

actual_FDR_05 <- FP/ (FP + TP)
# theoretical FDR

row_name <- rownames(output_05[[1]])

is_positive <- (grepl("_YEAST", row_name) & output_05[[1]]$pLtROPE > 90) |
  (grepl("_ECOLI", row_name) & output_05[[1]]$pGtROPE > 90) |
  (grepl("_HUMAN", row_name) & (output_05[[1]]$pLtROPE >= 90 | output_05[[1]]$pGtROPE >= 90)) |
  (grepl("_YEAST", row_name) & output_05[[1]]$pGtROPE >= 90) |
  (grepl("_ECOLI", row_name) & output_05[[1]]$pLtROPE >= 90)

# Select the positive rows
positive_rows <- output_05[[1]][is_positive, ]
positive_rows <- positive_rows[!is.na(positive_rows$Median), ]

result_sums <- sapply(1:nrow(positive_rows), function(i) {
  if (positive_rows$pLtROPE[i] >= 90) {
    # Sum pInROPE and pGtROPE when pLtROPE > 90
    100 - positive_rows$pLtROPE[i]
  } else if (positive_rows$pGtROPE[i] >= 90) {
    # Sum pInROPE and pLtROPE when pGtROPE > 90
    100 - positive_rows$pGtROPE[i]
  } else {
    # For cases that don't meet either condition (shouldn't happen with your selection)
    NA
  }
})

theoretical_FDR_05 <- mean(result_sums) / 100



actual_FDR <- c(actual_FDR_01, actual_FDR_02, actual_FDR_03, actual_FDR_04, actual_FDR_05) * 100
theoretical_FDR <- c(theoretical_FDR_01, theoretical_FDR_02, theoretical_FDR_03, theoretical_FDR_04, theoretical_FDR_05) * 100


# plot alpha = 0.05
x_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)

par(lwd = 2, cex.axis = 1.2, cex.lab = 1.3, font.lab = 2,
    mar = c(5, 5, 4, 2) + 0.1, mgp = c(3.5, 1, 0),
    family = "sans")

# Plot with thicker lines and improved styling
plot(x_values, actual_FDR,
     type = "b",  # "b" for both points and lines
     pch = 19,    # Solid circle points
     col = "blue3",
     lwd = 3,     # Thicker lines
     cex = 1.5,   # Larger points
     xlab = "Delta",
     ylab = "FDR (%)",
     ylim = c(0, max(c(actual_FDR, theoretical_FDR)) * 1.1),
     xaxt = "n",  # Remove default x-axis
     las = 1,     # Horizontal y-axis labels
     font.lab = 2,
     cex.lab = 1.4,
     cex.axis = 1.2,
     bty = "l",   # L-shaped box
     panel.first = grid(lwd = 1, col = "gray90", lty = 3))  # Add grid first

# Add theoretical FDR with thicker lines
points(x_values, theoretical_FDR,
       type = "b",
       pch = 17,   # Triangle points
       col = "red3",
       lwd = 3,    # Thicker lines
       cex = 1.5)  # Larger points

# Custom x-axis with thicker line
axis(1, at = x_values, labels = x_values,
     lwd = 2,      # Thicker axis line
     lwd.ticks = 2, # Thicker ticks
     cex.axis = 1.2)

# Add y-axis with thicker line (since we removed default x-axis)
axis(2, lwd = 2, lwd.ticks = 2, cex.axis = 1.2, las = 1)

# Add legend with thicker border
legend("topright",
       legend = c("Actual FDR", "Theoretical FDR"),
       col = c("blue3", "red3"),
       pch = c(19, 17),
       pt.cex = 1.5,  # Larger legend points
       lty = 1,
       lwd = 3,       # Thicker legend lines
       cex = 1.2,     # Larger legend text
       bg = "white",  # White background
       box.lwd = 2,   # Thicker legend border
       inset = 0.02)  # Slight inset from edge

# Add a box around the plot with thick line
box(lwd = 3)

# Reset graphical parameters
par(lwd = 1, cex.axis = 1, cex.lab = 1, font.lab = 1,
    mar = c(5, 4, 4, 2) + 0.1, mgp = c(3, 1, 0),
    bty = "o", family = "")

dev.off()




