# ===== Clear environment variables =====
rm(list = ls())

# ===== Set working directory to the folder containing the CSV file =====
setwd("D:\\Results\\GeneticDistance_VE_ExternalValidation_SensitivityAnalysis")

# ===== Load required R packages =====
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# ===== Read CSV file =====
data <- read_csv("Predicted_vs_ObservedEfficacy_SensitivityAnalysis.csv", show_col_types = FALSE)

# ===== Adjust factor levels in VE column (ensure order: lower, Mean, upper) =====
data$VE <- factor(data$VE, levels = c("95%CI_lower", "Mean", "95%CI_upper"))

# ===== Replace "+" with "+\n" in TechnicalName for better visualization =====
data$TechnicalName <- str_replace(data$TechnicalName, "\\+", "+\n")

# ===== Set position_dodge width (slightly closer) =====
position_adjust <- position_dodge(width = 0.2)

# ===== Plot =====
p <- ggplot(data, aes(x = TechnicalName, y = PredictedEfficacy, group = VE)) +
  # Plot predicted efficacy for each vaccine
  geom_point(aes(color = VE), size = 3, position = position_adjust) +
  # Plot predicted efficacy confidence intervals
  geom_errorbar(aes(ymin = Lower_Pred, ymax = Upper_Pred, color = VE), 
                width = 0.1, position = position_adjust) +
  # Rotate vaccine names 90 degrees
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14)) +
  # Customize grid and background
  theme(panel.grid.major.x = element_blank(),  
        panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted", size = 0.3),
        panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.3),
        plot.margin = unit(c(2, 0.5, 0.5, 0.5), "cm")) +  
  # Set y-axis range
  scale_y_continuous(limits = c(-20, 112), breaks = seq(0, 100, 20)) +
  # Customize legend labels, colors, and title
  scale_color_manual(values = c("95%CI_lower" = "#9B870C", "Mean" = "#FF69B4", "95%CI_upper" = "#006400"),
                     labels = c("95%CI_lower" = "95% CI Lower", 
                                "Mean" = "Predicted efficacy", 
                                "95%CI_upper" = "95% CI Upper")) +
  # Axis labels and legend title
  labs(y = "Predicted efficacy (%)", x = NULL, color = "Model-derived efficacy") +
  # Customize font size and legend style
  theme(axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA))

# ===== Print plot =====
print(p)

# ===== Save high-resolution PNG =====
ggsave("C:/Users/Desktop/Efficacy_Comparison_HighRes.png", plot = p, width = 8, height = 6, dpi = 1200)
