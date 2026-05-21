pg<-fread("report.pg_matrix.tsv")
pr<-fread("report.pr_matrix.tsv")
t <- table(unique(pr[,c('Protein.Group','Stripped.Sequence')])$Protein.Group)
pg$Peptide.Count <- t[match(pg$Protein.Group,names(t))]
#Filter by protein groups with <2 peptides
df.prot <- filter(pg, Peptide.Count > 1)
df.prot2 <- df.prot[!grepl("Cont_",df.prot$Protein.Group),]
df.prot.all<- df.prot2[, 5:15]
protein.matrix_all <- log2(as.matrix(df.prot.all))
log2all.df <- as.data.frame(protein.matrix_all)
unique_names <- make.names(df.prot2$Protein.Names, unique = TRUE)
rownames(log2all.df) <- unique_names

library(matrixStats)
group_means_all <- c()
group_sd_all <- c()
for (i in 1:(ncol(log2all.df) %/% r)) {
  #na_counts_x_i <- c()
  group_i <- log2all.df[, (r * (i - 1) + 1):(r * i)]
  group_i_mean <- rowMeans(group_i, na.rm = TRUE)
  group_i_var   <- rowVars(as.matrix(group_i), na.rm = TRUE)
  #group_i_RSD  <- group_i_sd / group_i_mean
  group_i_sd <- sqrt(group_i_var)
  group_means_all <- c(group_means_all, group_i_mean)
  group_sd_all  <- c(group_sd_all, group_i_sd)
}

# Set margins (mar) and text size (cex.axis, cex.lab) for publication
par(mar = c(5, 5, 4, 2) + 0.1,  # Adjust margins: bottom, left, top, right
    mgp = c(3, 1, 0),           # Adjust axis title (3), labels (1), and line (0) positions
    cex.axis = 1.2,             # Increase size of axis numbers
    cex.lab = 1.4,              # Increase size of axis labels
    las = 1,                    # Make axis labels always horizontal
    pty = "s")                  # This makes the plot region SQUARE ("s" for square)

# Create the plot
plot(group_means_all, group_sd_all,
     xlab = "Group Mean Intensity",  # x-axis label
     ylab = "Standard Deviation",    # y-axis label
     pch = 16,                       # Solid circle points
     col = "steelblue",              # Point color
     cex = 0.5,                      # Increase point size
     ylim = c(0, 4),                 # y-axis limits
     xlim = c(10,33), # x-axis limits, with 5% padding
     frame.plot = FALSE,             # Remove the box around the plot
     axes = FALSE)                   # Turn off default axes

# Add custom axes with thicker lines and longer tick marks
axis(1, lwd = 2, lwd.ticks = 2)  # Bottom axis (x-axis)
axis(2, lwd = 2, lwd.ticks = 2)  # Left axis (y-axis)
dev.off()
