library(readxl)
library(dplyr)
library(ggplot2)

# 1. Load data
df <- read_excel("/Users/ekaskaur/Desktop/biol4000 data.xlsx")

# 2. Select + clean the 7 morphological traits 
# (columns are coerced to numeric because some contain stray "NA" text)
traits <- df %>%
  transmute(
    species               = `geno.species`,
    `Glume ratio`          = as.numeric(`Glume ratio`),
    `Glume to seed ratio`  = as.numeric(`Glume to seed ratio`),
    `Lower branch length`  = as.numeric(`Lower branch length`),
    `Panicle length`       = as.numeric(`Pinnacle length`),   
    `Num of culms`         = as.numeric(`Num of culms`),
    `Culm height (low)`    = as.numeric(`Culm height (low)`),
    `Culm height (high)`   = as.numeric(`Culm height (high)`)
  ) %>%
  filter(species %in% c("FH", "HY", "FC"))   # drop blank/"NA" rows

# 3. Mean-impute missing trait values 
# Several specimens are missing culm-height / culm-count measurements.
# We fill these with the column mean so every labelled specimen can be
# plotted (this is what reproduces the full point count in the figure).
trait_mat <- as.matrix(traits[, -1])
for (j in seq_len(ncol(trait_mat))) {
  trait_mat[is.na(trait_mat[, j]), j] <- mean(trait_mat[, j], na.rm = TRUE)
}

# 4. Run PCA (scaled/standardised traits) 
pca <- prcomp(trait_mat, scale. = TRUE)
var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)

scores <- as.data.frame(pca$x[, 1:2])
scores$species <- factor(traits$species, levels = c("FH", "HY", "FC"))

# 5. Build scaled loading arrows
loadings <- as.data.frame(pca$rotation[, 1:2])
colnames(loadings) <- c("PC1", "PC2")
loadings$trait <- rownames(loadings)

arrow_scale <- 0.8 * min(
  (max(scores$PC1) - min(scores$PC1)) / (max(loadings$PC1) - min(loadings$PC1)),
  (max(scores$PC2) - min(scores$PC2)) / (max(loadings$PC2) - min(loadings$PC2))
)
loadings$xend <- loadings$PC1 * arrow_scale
loadings$yend <- loadings$PC2 * arrow_scale
loadings$xlab <- loadings$xend * 1.18
loadings$ylab <- loadings$yend * 1.18

# nudge a couple of labels that sit too close to each other / their arrow
nudge <- list(
  `Panicle length`      = c(dx = 0.05, dy = 0.35),
  `Lower branch length` = c(dx = 0.15, dy = -0.35)
)
for (nm in names(nudge)) {
  i <- which(loadings$trait == nm)
  loadings$xlab[i] <- loadings$xlab[i] + nudge[[nm]]["dx"]
  loadings$ylab[i] <- loadings$ylab[i] + nudge[[nm]]["dy"]
}

# 6. Colours
cols <- c(FH = "#4CAF50", HY = "#4169E1", FC = "#F08080")

# 7. Plot 
p <- ggplot(scores, aes(PC1, PC2, color = species, fill = species)) +
  stat_ellipse(geom = "polygon", aes(color = species), alpha = 0.15, level = 0.95, linewidth = 0.3) +
  stat_ellipse(geom = "path", linewidth = 0.8, level = 0.95) +
  geom_point(size = 3) +
  geom_segment(
    data = loadings, inherit.aes = FALSE,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.018, "npc"), type = "open"),
    linewidth = 0.6, color = "black"
  ) +
  geom_text(
    data = loadings, inherit.aes = FALSE,
    aes(x = xlab, y = ylab, label = trait),
    fontface = "bold", size = 4.2, color = "black"
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(
    title = "PCA of Morphological Traits",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    color = "geno.species", fill = "geno.species"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    axis.line = element_line(linewidth = 0.5)
  )

ggsave("pca_biplot.png", p, width = 9, height = 7, dpi = 300, bg = "white")
cat("Variance explained PC1/PC2:", var_exp, "\n")
cat("n per species:\n"); print(table(scores$species))
