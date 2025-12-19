library(ggplot2)
library(reshape2)
library(readr)


# Load simulated data given by ChatGPT
z_scores <- read_csv("~/Downloads/bubble_heatmap_matrix.csv")
fdr_vals <- read_csv("~/Downloads/bubble_heatmap_fdr.csv")

# Convert to matrices with features as rownames
Z <- as.data.frame(z_scores)
FDR <- as.data.frame(fdr_vals)
rownames(Z) <- Z[[1]]
Z[[1]] <- NULL
rownames(FDR) <- FDR[[1]]
FDR[[1]] <- NULL

# Melt into long format ----
Z_long <- melt(as.matrix(Z), varnames = c("Feature", "Gene"), value.name = "Zscore")
FDR_long <- melt(as.matrix(FDR), varnames = c("Feature", "Gene"), value.name = "FDR")

# Merge into one data frame
df <- merge(Z_long, FDR_long, by = c("Feature", "Gene"))

# Compute bubble size as -log10(FDR), cap for readability
df$neglogFDR <- -log10(pmax(df$FDR, 1e-6))
cap <- quantile(df$neglogFDR, 0.95)
df$neglogFDR <- pmin(df$neglogFDR, cap)

####GRAPH FOR HEATMAP 
p <- ggplot(df, aes(x = Gene, y = Feature)) +
  geom_point(aes(size = neglogFDR, color = Zscore)) +   # no border
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "Feature z-score\n(vs NTC)") +
  scale_size_continuous(name = expression(-log[10]*"(FDR)")) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  ggtitle("Neuron screen: multi-feature phenotype signatures")

print(p)

library(ggplot2)
library(dplyr)
library(ggrepel)

set.seed(7)

# Set some simulated standard deviations, help from Chat
sd_ntc   <- 0.20   #tighter NTC cluster (↓ = tighter)
sd_geneX <- 0.60   #spread of gene effects on x
sd_noise <- 0.35   #noise orthogonal to correlation line (↓ = tighter)
rho      <- 0.60   #correlation strength between x and y (↑ = more diagonal)

n_genes <- 220 #random numbers to simulate my data
n_ntc   <- 60

# Simulate normalized features aka z with NTCs near 0,0
ntc <- tibble(
  x = rnorm(n_ntc, 0, sd_ntc),
  y = rnorm(n_ntc, 0, sd_ntc),
  class = "NTC",
  gene  = paste0("NTC_", seq_len(n_ntc))
)

# Genes: y = rho*x + noise
gx <- rnorm(n_genes, 0, sd_geneX)
gy <- rho * gx + rnorm(n_genes, 0, sd_noise)
genes <- tibble(
  x = gx, y = gy, class = "Gene",
  gene = paste0("GENE_", sprintf("%03d", seq_len(n_genes)))
)

# A few labeled hits
hits <- tribble(
  ~gene,   ~x,    ~y,
  "MAP7",   1.8,   2.2,
  "CDH2",  -2.1,  -1.4,
  "NLGN1", -2.6,  -2.7
) %>% mutate(class = "Gene")

df <- bind_rows(genes, ntc, hits)

# Correlation & regression (genes only) 
r_val <- cor(df %>% filter(class=="Gene") %>% pull(x),
             df %>% filter(class=="Gene") %>% pull(y))
fit <- lm(y ~ x, data = df %>% filter(class=="Gene"))

                       ###### GRAPH FOR COVERAGE #####
plot2 <- ggplot(df, aes(x = x, y = y)) +
  # genes (gray circles, slightly larger, semi-opaque)
  geom_point(data = df %>% filter(class=="Gene"),
             color = "gray50", size = 2.6, alpha = 0.55) +
  # NTCs (tight orange cluster)
  geom_point(data = df %>% filter(class=="NTC"),
             color = "#F28E2B", size = 3.0, alpha = 0.95) +
  # Regression line (genes)
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2],
              color = "black", linewidth = 0.7) +
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray55") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray55") +
  # Labels for example hits
  geom_text_repel(data = hits, aes(label = gene),
                  size = 3, fontface = "bold", color = "black",
                  box.padding = 0.15, point.padding = 0.15, max.overlaps = 100) +
  labs(
    title = "Normalized phenotypes (lower-variance demo)",
    subtitle = paste0("r = ", sprintf("%.2f", r_val)),
    x = "# of foci (normalized z-score)",
    y = "Area of foci (normalized z-score)"
  ) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) +  # tighter frame
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")

print(plot2)

library(ggplot2)
library(dplyr)

set.seed(2025)

# Simulate gene-level effects for two cell types
n_genes <- 1200

df <- tibble(
  Gene = paste0("GENE_", 1:n_genes),
  Neuron_effect   = rnorm(n_genes, 0, 0.6),
  Myotube_effect  = rnorm(n_genes, 0, 0.6)
)

# Add a few true "differential" hits
hits <- c("MAP7", "CDH2", "NLGN1", "SFPQ", "PTEN")
df$Neuron_effect[df$Gene %in% hits] <- c(1.5, 1.2, -1.4, 0.9, -1.2)
df$Myotube_effect[df$Gene %in% hits] <- c(0.1, -0.2, 0.1, -0.4, 0.3)

# Calc differential effect
df <- df %>%
  mutate(
    Diff = Neuron_effect - Myotube_effect,
    pval = 2 * pnorm(-abs(Diff / 0.35)),    # fake p-values (assume SEM=0.35)
    negLog10P = -log10(pval),
    sig = ifelse(abs(Diff) > 1 & pval < 0.05, "Significant", "NS")
  )

                      ###### VOLCANO PLOT #####
ggplot(df, aes(x = Diff, y = negLog10P)) +
  geom_point(aes(color = sig), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("gray70", "#D55E00")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  labs(
    x = expression(Delta*" effect (Neuron – Myotube)"),
    y = expression(-log[10]*"(p-value)"),
    title = "Neuron vs myotube differential effects in # of foci"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  # Label example hits
  geom_text(
    data = subset(df, Gene %in% hits),
    aes(label = Gene),
    size = 3, fontface = "bold", vjust = -0.8
  )

install.packages("pwr")   
library(pwr)

pwr.t.test(
  d = 0.30,
  sig.level = 0.05,
  power = 0.80,
  type = "two.sample",
  alternative = "two.sided"
)

set.seed(123)

# number of guides & genes
n_guides <- 8431
n_genes  <- 2627
guides_per_gene <- 3

# Simulate coverage (cells per guide)
# Negative binomial gives realistic overdispersed coverage
guide_cov <- rnbinom(n_guides, size = 50, mu = 100)

# Put into a data frame
guide_data <- data.frame(
  guide_id = paste0("g", 1:n_guides),
  gene_id  = rep(paste0("gene", 1:n_genes), each = guides_per_gene)[1:n_guides],
  coverage = guide_cov
)

head(guide_data)

library(dplyr)

gene_data <- guide_data %>%
  group_by(gene_id) %>%
  summarize(mean_coverage = mean(coverage))

head(gene_data)


                         ###### HISTOGRAM FOR COVERAGE PLOT #####
ggplot(guide_data, aes(x = coverage)) +
  geom_histogram(binwidth = 10, fill = "steelblue", alpha = 0.8) +
  geom_vline(xintercept = 100, color = "red", linetype = "dashed", size = 1) +
  theme_classic() +
  labs(
    title = "Coverage Distribution per Guide",
    x = "Cells per Guide",
    y = "Number of Guides"
  )

ggplot(gene_data, aes(x = mean_coverage)) +
  geom_histogram(binwidth = 10, fill = "orchid", alpha = 0.8) +
  geom_vline(xintercept = 100, color = "red", linetype = "dashed", size = 1) +
  theme_classic() +
  labs(
    title = "Coverage Distribution per Gene (Mean of 3 Guides)",
    x = "Mean Cells per Gene",
    y = "Number of Genes"
  )


library(MASS)
library(dplyr)
library(ggplot2)

set.seed(123)

Sigma <- matrix(c(
  0.4, 0.2, 0.2, 0.1,
  0.2, 0.4, 0.2, 0.1,
  0.2, 0.2, 0.4, 0.1,
  0.1, 0.1, 0.1, 0.4
), nrow=4, byrow=TRUE)

mean_NT    <- c(5, 50, 100, 300)
mean_DDX5  <- c(5.5, 52, 103, 303)
mean_MBNL1 <- c(6.2, 55, 110, 305)

n_guides_NT    <- 60
n_guides_DDX5  <- 10
n_guides_MBNL1 <- 10
cells_per_guide <- 80

sim_group <- function(n_guides, mean_vec, label) {
  do.call(rbind, lapply(1:n_guides, function(g) {
    X <- mvrnorm(cells_per_guide, mean_vec, Sigma)
    df <- as.data.frame(X)
    df$guide <- paste0(label, "_", g)
    df$group <- label
    df
  }))
}

NT    <- sim_group(n_guides_NT,    mean_NT,    "NT")
DDX5  <- sim_group(n_guides_DDX5,  mean_DDX5,  "DDX5 KO")
MBNL1 <- sim_group(n_guides_MBNL1, mean_MBNL1, "MBNL1 KO")

all_cells <- rbind(NT, DDX5, MBNL1)
colnames(all_cells)[1:4] <- c("foci_count","foci_area","intensity","nuc_area")

# collapse to per-guide means
guide_level <- all_cells %>%
  group_by(group, guide) %>%
  summarize(
    foci_count = mean(foci_count),
    foci_area  = mean(foci_area),
    intensity  = mean(intensity),
    nuc_area   = mean(nuc_area),
    .groups = "drop"
  )

# PCA on guide means
pca <- prcomp(guide_level[,3:6], scale. = TRUE)
scores <- as.data.frame(pca$x)
scores$group <- guide_level$group

ggplot(scores, aes(PC1, PC2, color = group)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values=c("NT"="#1f77b4","DDX5 KO"="#ff7f0e","MBNL1 KO"="#9467bd")) +
  theme_classic(base_size = 16) +
  labs(title = "PCA of Guide-Level Optical Phenotypes",
       x = "PC1", y = "PC2")


# PC1 loadings = contribution of each feature to PC1
loadings_df <- data.frame(
  feature = rownames(pca$rotation),
  loading = pca$rotation[,1]          # PC1 is column 1
)

# Order by absolute loading (biggest contributors on top)
loadings_df <- loadings_df %>%
  arrange(desc(abs(loading)))

loadings_df$feature <- factor(loadings_df$feature,
                              levels = loadings_df$feature)  # keep this order

                             ###### PCA ANALYSIS PLOT #####

ggplot(loadings_df, aes(x = loading, y = feature)) +
  geom_col(fill = "steelblue") +
  theme_classic(base_size = 14) +
  labs(
    x = "PC1 loadings",
    y = "Feature",
    title = "Contribution of Optical Features to PC1"
  )
ggplot(loadings_df, aes(x = loading, y = feature)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(loading, 3)),
            hjust = ifelse(loadings_df$loading > 0, 1.05, -0.05),
            color = "white", size = 3) +
  theme_classic(base_size = 14) +
  labs(
    x = "PC1 loadings",
    y = "Feature",
    title = "Contribution of Optical Features to PC1"
  )


#Create a data frame of your 4 features ----
df <- data.frame(
  Feature = c("Nucleus area",
              "Foci intensity",
              "Foci number",
              "Foci area"),
  Loading = c(0.21, 0.18, 0.15, 0.17)   # <- replace with real values
)

#Reorder features by loading ----
df <- df %>% arrange(Loading) %>%
  mutate(Feature = factor(Feature, levels = Feature))

                       ###### PCA LOADINGS PLOT #####
ggplot(df, aes(x = Loading, y = Feature)) +
  geom_bar(stat = "identity", fill = "#8FB3E8", width = 0.7) +
  
  # Text labels
  geom_text(aes(label = round(Loading, 3)),
            hjust = -0.15,
            size = 4.5) +
  


#Simulate guide-level coverage (around 175)
n_guides <- 8431
guide_cov <- rnbinom(n_guides, size = 80, mu = 175)

guide_df <- data.frame(guide = paste0("g", 1:n_guides),
                       coverage = guide_cov)

#Gene-level coverage (3 guides per gene)
n_genes <- 2627
guides_per_gene <- 3

guide_df$gene <- rep(paste0("gene", 1:n_genes), each = guides_per_gene)[1:n_guides]

gene_df <- guide_df %>%
  group_by(gene) %>%
  summarise(mean_cov = mean(coverage))

#Coverage per guide
p1 <- ggplot(guide_df, aes(x = coverage)) +
  geom_histogram(binwidth = 10, fill = "#6A93C9", alpha = 0.8) +
  geom_vline(xintercept = 175, color = "red", linetype = "dashed", size = 1) +
  theme_classic(base_size = 14) +
  labs(
    title = "Coverage Distribution per Guide",
    x = "Cells per Guide",
    y = "Number of Guides"
  )
p1

install.packages("UpSetR")   # once
library(UpSetR)
set.seed(123)

hits_for_upset <- data.frame(
  gene       = paste0("GENE_", sprintf("%03d", 1:150)),
  Neuron     = rbinom(150, 1, 0.20),  # 20% neuron hits
  Fibroblast = rbinom(150, 1, 0.15),  # 15% fibro hits
  Myotube    = rbinom(150, 1, 0.10)   # 10% myotube hits
)

UpSet(
  hits_for_upset,
  sets        = c("Neuron", "Fibroblast", "Myotube"),
  nintersects = 10,
  order.by    = "freq",
  keep.order  = TRUE,
  mainbar.y.label = "Number of hit genes",
  sets.x.label    = "Hits per cell type"
)

