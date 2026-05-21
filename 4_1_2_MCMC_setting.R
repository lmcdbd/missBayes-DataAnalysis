
library(data.table)
library(matrixStats)
library(ggplot2)
library(dplyr)
library(limma)
library(missBayes)

# import CQE dataset 
pg<-fread("/Users/mengchunli/Documents/work/Paper revision/paper documents/CQE/report.pg_matrix.tsv")
pr<-fread("/Users/mengchunli/Documents/work/Paper revision/paper documents/CQE/report.pr_matrix.tsv")
t <- table(unique(pr[,c('Protein.Group','Stripped.Sequence')])$Protein.Group)
pg$Peptide.Count <- t[match(pg$Protein.Group,names(t))]
metadata <- fread("/Users/mengchunli/Documents/work/Paper revision/paper documents/CQE/metadata.txt", header = TRUE)
#Filter by protein groups with <2 peptides
df.prot <- filter(pg, Peptide.Count > 1)
df.prot2 <- df.prot[!grepl("Cont_",df.prot$Protein.Group),]
df.prot.all<- df.prot2[, 5:29]
protein.matrix_all <- log2(as.matrix(df.prot.all))
log2all.df <- as.data.frame(protein.matrix_all)
rownames(log2all.df) <- df.prot2$Protein.Names

# run model
group <- as.factor(metadata$Condition)
comparison_A <- makeContrasts("G70_20_10 - G70_10_20", levels = levels(group))
comparison_B <- makeContrasts("G35_25_40 - G35_40_25", levels = levels(group))
set.seed(746)
output_A <- BayesMissingModel(values = log2all.df, groups = group, comparisons = comparison_A, parallel = TRUE,
                             n.adapt = 500, burn.in = 500, n.iter = 10000, n.chains = 2, mcmcDiag = TRUE)
output_B <- BayesMissingModel(values = log2all.df, groups = group, comparisons = comparison_B, parallel = TRUE,
                              n.adapt = 500, burn.in = 500, n.iter = 10000, n.chains = 2, mcmcDiag = TRUE)

# 1. missing proportion distribution with different convergence ----

output_B[[1]]$missing_proportion <- log2all.df[,1:10] %>% is.na() %>% rowMeans()
# plotting data
plot_df_B <- output_B[[1]] %>%
  select(Convergence, missing_proportion) %>%
  filter(
    !is.na(Convergence),
    !is.na(missing_proportion)
  ) %>%
  mutate(
    Convergence = factor(
      Convergence,
      levels = c("Strong", "Moderate", "Poor")
    )
  )

plot_df_combined <- bind_rows(plot_df_A, plot_df_B)

# Check counts
table(output_B[[1]]$Convergence, useNA = "ifany")

ggplot(plot_df_combined,
       aes(x = Convergence,
           y = missing_proportion,
           fill = Convergence)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    alpha = 0.9
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,      # diamond
    size = 3,
    fill = "white",
    color = "black"
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Convergence Status",
    y = "Missing Proportion"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )
table(output_B[[1]]$Convergence, useNA = "ifany")

# 2. 95% HDI width distribution with different convergence ----
output_A[[1]]$HDI_width <- output_A[[1]]$HDI_High - output[[1]]$HDI_Low
output_B[[1]]$HDI_width <- output_B[[1]]$HDI_High - output_B[[1]]$HDI_Low

# Prepare plotting data
plot_hdi_A <- output_A[[1]] %>%
  select(Convergence, HDI_width) %>%
  filter(
    !is.na(Convergence),
    !is.na(HDI_width)
  ) %>%
  mutate(
    Dataset = "A",
    Convergence = factor(
      Convergence,
      levels = c("Strong", "Moderate", "Poor")
    )
  )

plot_hdi_B <- output_B[[1]] %>%
  select(Convergence, HDI_width) %>%
  filter(
    !is.na(Convergence),
    !is.na(HDI_width)
  ) %>%
  mutate(
    Dataset = "B",
    Convergence = factor(
      Convergence,
      levels = c("Strong", "Moderate", "Poor")
    )
  )

# Combine datasets
plot_hdi <- bind_rows(plot_hdi_A, plot_hdi_B)

# Plot
ggplot(plot_hdi,
       aes(x = Convergence,
           y = HDI_width,
           fill = Convergence)) +
  
  geom_violin(trim = FALSE,
              alpha = 0.7,
              linewidth = 0.8) +
  
  geom_boxplot(width = 0.15,
               outlier.shape = NA,
               alpha = 0.9) +
  
  stat_summary(fun = median,
               geom = "point",
               shape = 23,
               size = 3,
               fill = "white") +
  
  scale_fill_brewer(palette = "Set2") +
  
  labs(
    x = "Convergence Status",
    y = "95% HDI Width"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )

# 3. reproducibility of results  ----

set.seed(746)
output_B_4_NP <- BayesMissingModel(values = log2all.df, groups = group, comparisons = comparison, parallel = FALSE,
                                n.adapt = 500, burn.in = 500, n.iter = 10000, n.chains = 2, mcmcDiag = TRUE)


runs <- list(output_B[[1]], output_B_2[[1]], output_B_3_NP[[1]], output_B_4_NP[[1]])

# check if any proteins are labeled differently among four runs
add_label <- function(df) {
  df$Label <- ifelse(
    df$pLtROPE >= 95, "Down",
    ifelse(df$pGtROPE >= 95, "Up", "NS")
  )
  df
}

# Apply to each run
runs[[1]] <- add_label(runs[[1]])
runs[[2]] <- add_label(runs[[2]])
runs[[3]] <- add_label(runs[[3]])
runs[[4]] <- add_label(runs[[4]])

# Put all labels together
label_df <- data.frame(
  Run1 = runs[[1]]$Label,
  Run2 = runs[[2]]$Label,
  Run3 = runs[[3]]$Label,
  Run4 = runs[[4]]$Label,
  row.names = rownames(runs[[1]])
)

# Identify proteins with inconsistent labels
inconsistent <- apply(label_df, 1, function(x) length(unique(x)) > 1)

# Assuming each run contains columns:
# Median, missing_proportion

# Common proteins across all runs
common_proteins <- Reduce(intersect, lapply(runs, rownames))

# Extract logFC values
logFC_mat <- sapply(runs, function(x) {
  x[common_proteins, "Median"]
})

colnames(logFC_mat) <- paste0("Run", seq_along(runs))

# Maximum pairwise difference across all runs
max_diff <- apply(logFC_mat, 1, function(x) {
  max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
})

missing_prop <- output_B[[1]][common_proteins, "missing_proportion"]

plot_df <- data.frame(
  Protein = common_proteins,
  MissingProportion = missing_prop,
  MaxLogFCDiff = max_diff
)

ggplot(plot_df, aes(MissingProportion, MaxLogFCDiff)) +
  geom_point(alpha = 0.5) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(linewidth = 1),
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(0.25, "cm")
  ) +
  labs(
    x = "Missing proportion",
    y = expression("|" * Delta * " logFC|")
  )

# 4. comparing two missing models ----
output_A_logi <- BayesMissingModel(values = log2all.df, groups = group, comparisons = comparison_A, parallel = TRUE,
                              n.adapt = 500, burn.in = 500, n.iter = 10000, n.chains = 2, mcmcDiag = TRUE, threshold = 1)
df.ps <- output_A[[1]] %>%
  select(Median, pGtROPE, pLtROPE) %>%
  dplyr::rename(logFC_PS = Median, pGtROPE_PS = pGtROPE, pLtROPE_PS = pLtROPE)

df.gl <- output_A_logi[[1]] %>%
  select(Median, pGtROPE, pLtROPE) %>%
  dplyr::rename(logFC_GL = Median, pGtROPE_GL = pGtROPE, pLtROPE_GL = pLtROPE)

df.ps$Protein <- rownames(df.ps)
df.gl$Protein <- rownames(df.gl)

df.gl$significant_GL <- df.gl$pGtROPE_GL > 95 | df.gl$pLtROPE_GL > 95
df.ps$significant_PS <- df.ps$pGtROPE_PS > 95 | df.ps$pLtROPE_PS > 95

merged <- df.ps %>%
  inner_join(df.gl, by = "Protein", suffix = c("_PS", "_GL"))

# Pearson correlation
cor_test <- cor.test(
  merged$logFC_PS,
  merged$logFC_GL,
  method = "pearson"
)

r_value <- unname(cor_test$estimate)

lims <- range(
  c(merged$logFC_PS, merged$logFC_GL),
  na.rm = TRUE
)

ggplot(merged, aes(logFC_PS, logFC_GL)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(linetype = 2) +
  coord_equal(xlim = lims, ylim = lims) +
  annotate(
    "text",
    x = lims[1] + 0.05 * diff(lims),
    y = lims[2] - 0.05 * diff(lims),
    label = paste0("Pearson r = ", round(r_value, 3)),
    hjust = 0,
    vjust = 1,
    size = 5
  ) +
  labs(
    x = "Protein-specific cutoff logFC",
    y = "Global logistic logFC",
    title = NULL
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 1.0),
    axis.ticks = element_line(linewidth = 1.0),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

table(merged$significant_PS, merged$significant_GL)

set_ps <- merged$Protein[merged$significant_PS]
set_gl <- merged$Protein[merged$significant_GL]

venn_list <- list(
  "Protein-specific cutoff" = set_ps,
  "Global logistic" = set_gl
)

ggVennDiagram(
  venn_list,
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme_void()
