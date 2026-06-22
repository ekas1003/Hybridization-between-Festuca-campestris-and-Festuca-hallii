library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

# 1. Load data
file_path <- "~/Desktop/biol4000 data.xlsx"   

raw <- read_excel(file_path, sheet = "Sheet1")
raw <- raw[-1, ]   # drop the blank row sitting directly under the header

df <- raw %>%
  dplyr::select(
    geno.species,
    mean.qscore.FC,
    `Glume ratio`,
    `Glume to seed ratio`,
    `Lower branch length`,
    `Pinnacle length`
  ) %>%
  dplyr::filter(geno.species %in% c("FH", "HY", "FC")) %>%
  mutate(
    geno.species          = factor(geno.species, levels = c("FH", "HY", "FC")),
    mean.qscore.FC         = as.numeric(mean.qscore.FC),
    `Glume ratio`           = as.numeric(`Glume ratio`),
    `Glume to seed ratio`   = as.numeric(`Glume to seed ratio`),
    `Lower branch length`   = as.numeric(`Lower branch length`),
    `Pinnacle length`       = as.numeric(`Pinnacle length`)
  )

# 2. Colors
fill_colors <- c(FH = "#4CAF50", HY = "#6BAED6", FC = "#F08080")

# 3. build one scatter panel with fitted regression line
make_scatter_panel <- function(data, yvar, ylab, title, xvar = "mean.qscore.FC") {

  d <- data %>% dplyr::filter(!is.na(.data[[yvar]]), !is.na(.data[[xvar]]))

  
  fit_data <- data.frame(x = d[[xvar]], y = d[[yvar]])
  model <- lm(y ~ x, data = fit_data)

  cat("\n====", title, "====\n")
  print(summary(model))

  intercept <- coef(model)[1]
  slope     <- coef(model)[2]

  ggplot(d, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(aes(color = geno.species), size = 3) +
    geom_abline(intercept = intercept, slope = slope,
                color = "black", linewidth = 1.2) +
    scale_color_manual(values = fill_colors, name = "geno.species") +
    labs(title = title, x = "Mean Q-score (FC)", y = ylab) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10)
    )
}

# 4. Build the 4 panels 
p1 <- make_scatter_panel(df, "Pinnacle length",     "Panicle length (mm)",      "Panicle length vs Mean Q-score (FC)")
p2 <- make_scatter_panel(df, "Glume ratio",         "Glume ratio",              "Glume ratio vs Mean Q-score (FC)")
p3 <- make_scatter_panel(df, "Glume to seed ratio", "Glume to seed ratio",      "Glume to seed ratio vs Mean Q-score (FC)")
p4 <- make_scatter_panel(df, "Lower branch length", "Lower branch length (mm)", "Lower branch length vs Mean Q-score (FC)")

# 5. Combine with one shared legend
combined <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined <- combined +
  plot_annotation(
    title = "Morphological Traits vs Mean Q-score (FC)",
    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  )

combined

ggsave("morphological_traits_vs_qscore.png", combined,
       width = 14, height = 10, dpi = 300, bg = "white")

