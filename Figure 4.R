# ===== Clear environment =====
rm(list = ls())

# ===== Load required libraries =====
library(ggplot2)
library(ggnewscale)

# ===== Set file path and load data =====
ve_df <- read.csv("D:/results/Jiangsu_VE_predictions.csv")

# Convert months column to date format
ve_df$months <- as.Date(ve_df$months, format = "%Y/%m/%d")

# Check conversion
print(head(ve_df$months))

# Remove missing values from numeric columns
ve_df <- ve_df[complete.cases(ve_df[, c("fit", "lwr", "upr")]), ]

# Check platforms column
unique_platforms <- unique(ve_df$platforms)

# ===== Variant key dates =====
id_date_delta     <- as.Date("2021-07-01", format = "%Y-%m-%d")
id_date_BA.5_1    <- as.Date("2022-02-01", format = "%Y-%m-%d")
id_date_BF.7      <- as.Date("2022-11-01", format = "%Y-%m-%d")
id_date_BA.5_2    <- as.Date("2022-12-01", format = "%Y-%m-%d")
id_date_XBB.1.16  <- as.Date("2023-02-01", format = "%Y-%m-%d")
id_date_EG.5.1    <- as.Date("2023-05-01", format = "%Y-%m-%d")

# ======================== Figure 2-4: One plot per platform ========================
platforms_list <- unique(ve_df$platforms)
platform_colors <- c('#FF6666', '#66CC99', "#6666FF", "#FFCC33")
platform_fills  <- alpha(c("#FF9999", "#33CC33", "mediumpurple2", "#FFFF00"), 0.3)

for (i in seq_along(platforms_list)) {
  plt_name <- platforms_list[i]
  plt_color <- platform_colors[i]
  plt_fill <- platform_fills[i]
  
  ve_sub <- ve_df[ve_df$platforms == plt_name, ]
  
  ve_plot_sub <- ggplot() +
    geom_line(data = ve_sub, aes(x = months, y = fit), color = plt_color, size = 1.0) +
    geom_ribbon(data = ve_sub, aes(x = months, ymin = lwr, ymax = upr), fill = plt_fill, alpha = 0.3) +
    
    theme_classic() +
    annotate("text", x = id_date_delta + 85, y = 118, label = "  Delta", hjust = 0, vjust = 1.5, size = 3.9, colour = "black") +
    annotate("text", x = id_date_BA.5_1 + 110, y = 118, label = "  BA.5", hjust = 0, vjust = 1.5, size = 3.9, colour = "black") +
    geom_vline(xintercept = id_date_BA.5_1, size = 0.7, linetype = "dotted", color = "grey36") +
    annotate("text", x = id_date_BF.7 + 1, y = 118, label = "BF.7", hjust = 0, vjust = 1.5, size = 3.9, colour = "black") +
    geom_vline(xintercept = id_date_BF.7, size = 0.7, linetype = "dotted", color = "grey36") +
    annotate("text", x = id_date_BA.5_2 + 10, y = 118, label = "  BA.5", hjust = 0, vjust = 1.5, size = 3.9, colour = "black") +
    geom_vline(xintercept = id_date_BA.5_2, size = 0.7, linetype = "dotted", color = "grey36") +
    annotate("text", x = id_date_XBB.1.16 + 10, y = 118, label = " XBB.1.16", hjust = 0, vjust = 1.5, size = 3.9, colour = "black") +
    geom_vline(xintercept = id_date_XBB.1.16, size = 0.7, linetype = "dotted", color = "grey36") +
    annotate("text", x = id_date_EG.5.1 + 18, y = 118, label = "  EG.5.1", hjust = 0, vjust = 1.5, size = 3.9, colour = "black") +
    geom_vline(xintercept = id_date_EG.5.1, size = 0.7, linetype = "dotted", color = "grey36") +
    geom_hline(yintercept = 50, size = 0.7, linetype = "dotted", color = "grey36") +
    geom_hline(yintercept = 0, size = 0.7, linetype = "dotted", color = "grey36") +
    scale_x_date(limits = as.Date(c("2021-07-01", "2023-07-31")),
                 breaks = seq(as.Date("2021-07-01"), as.Date("2023-07-31"), by = "1 month"),
                 date_labels = "%Y-%m", expand = c(0.02, 0), name = "") +
    scale_y_continuous(name = "Predicted efficacy (%)", limits = c(-20, 119), breaks = seq(-20, 100, 20)) +
    
    theme(
      axis.text.x = element_text(vjust = 0.8, hjust = 0.9, angle = 30, colour = "black", size = 13),
      axis.text.y = element_text(colour = "black", size = 14),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold")
    ) +
    ggtitle("")
  
  ggsave(
    filename = paste0("D:/results/Jiangsu_VE_dynamic_predictions/jiangsu_", i, "_", plt_name, ".png"),
    plot = ve_plot_sub,
    width = 10,
    height = 4,
    dpi = 1200
  )
}
