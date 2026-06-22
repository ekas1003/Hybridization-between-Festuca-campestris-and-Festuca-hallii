library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork) 

# 1. Load data
file_path <- "~/Desktop/biol4000 data.xlsx"

raw <- read_excel(file_path, sheet = "Sheet1")
raw <- raw[-1, ]   

df <- raw %>%
  dplyr::select(
    geno.species,
    `Glume ratio`,
    `Glume to seed ratio`,
    `Lower branch length`,
    `Pinnacle length`            
  ) %>%                          
  dplyr::filter(geno.species %in% c("FH", "HY", "FC")) %>%
  mutate(
    geno.species = factor(geno.species, levels = c("FH", "HY", "FC")),
    `Glume ratio`            = as.numeric(`Glume ratio`),
    `Glume to seed ratio`    = as.numeric(`Glume to seed ratio`),
    `Lower branch length`    = as.numeric(`Lower branch length`),
    `Pinnacle length`        = as.numeric(`Pinnacle length`)
  )

# 2. Colors
fill_colors <- c(FH = "#4CAF50", HY = "#4169E1", FC = "#F08080")

# 3. build one panel with stats + significance brackets
make_panel <- function(data, yvar, ylab, title, digits = 2, p_adjust = "bonferroni") {

  d <- data %>% dplyr::filter(!is.na(.data[[yvar]]))

  kw <- kruskal.test(d[[yvar]], d$geno.species)
  cat("\n====", title, "====\n"); print(kw)

  pw <- pairwise.wilcox.test(d[[yvar]], d$geno.species, p.adjust.method = p_adjust)
  print(pw)

  p_FH_HY <- round(pw$p.value["HY", "FH"], digits)
  p_FH_FC <- round(pw$p.value["FC", "FH"], digits)
  p_HY_FC <- round(pw$p.value["FC", "HY"], digits)

  ymax <- max(d[[yvar]], na.rm = TRUE)
  rng  <- ymax - min(d[[yvar]], na.rm = TRUE)
  step <- rng * 0.10
  tick <- step * 0.25

  y1 <- ymax + step * 1.0   # FH vs HY  (lowest bracket)
  y2 <- ymax + step * 2.4   # HY vs FC  (middle bracket)
  y3 <- ymax + step * 3.8   # FH vs FC  (full-span, highest bracket)

  ggplot(d, aes(x = geno.species, y = .data[[yvar]], fill = geno.species)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(aes(color = geno.species), width = 0.1, size = 2.2, alpha = 0.9) +
    scale_fill_manual(values = fill_colors, name = "geno.species") +
    scale_color_manual(values = fill_colors, guide = "none") +

    annotate("segment", x = 1, xend = 2, y = y1, yend = y1) +
    annotate("segment", x = 1, xend = 1, y = y1 - tick, yend = y1) +
    annotate("segment", x = 2, xend = 2, y = y1 - tick, yend = y1) +
    annotate("text", x = 1.5, y = y1 + step * 0.35, label = paste0(p_FH_HY), size = 4) +

    annotate("segment", x = 2, xend = 3, y = y2, yend = y2) +
    annotate("segment", x = 2, xend = 2, y = y2 - tick, yend = y2) +
    annotate("segment", x = 3, xend = 3, y = y2 - tick, yend = y2) +
    annotate("text", x = 2.5, y = y2 + step * 0.35, label = paste0(p_HY_FC), size = 4) +

    annotate("segment", x = 1, xend = 3, y = y3, yend = y3) +
    annotate("segment", x = 1, xend = 1, y = y3 - tick, yend = y3) +
    annotate("segment", x = 3, xend = 3, y = y3 - tick, yend = y3) +
    annotate("text", x = 2, y = y3 + step * 0.35, label = paste0(p_FH_FC), size = 4) +

    labs(title = title, x = "Species", y = ylab) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10)
    )
}

# 4. Build the 4 panels 
p1 <- make_panel(df, "Glume ratio",         "Glume Ratio",              "Glume ratio vs Species",          digits = 2, p_adjust = "bonferroni")
p2 <- make_panel(df, "Glume to seed ratio", "Glume to Seed Ratio",      "Glume to seed ratio vs Species",  digits = 3, p_adjust = "bonferroni")
p3 <- make_panel(df, "Lower branch length", "Lower Branch Length (mm)", "Lower branch length vs Species",  digits = 2, p_adjust = "bonferroni")
p4 <- make_panel(df, "Pinnacle length",     "Panicle Length (mm)",      "Panicle length vs Species",       digits = 2, p_adjust = "none")

# 5. Combine grid with one shared legend
combined <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined <- combined +
  plot_annotation(
    title = "Morphological Traits vs Species",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  )

combined

ggsave("morphological_traits_vs_species.png", combined,
       width = 14, height = 10, dpi = 300, bg = "white")

