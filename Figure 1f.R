# Clear environment variables
rm(list = ls())

# Set working directory to the path containing CSV file
setwd("D:\\Results\\NeutralizingAntibody_Protection_InternalExternalValidation_Calibration")

# Load required R packages
library(ggplot2)
library(dplyr)
library(readr)

# Read CSV file
data <- read_csv("Predicted_vs_Observed_Efficacy_ExternalValidation_Sensitivity.csv", show_col_types = FALSE)

# Adjust factor levels in NAb column to ensure order: lower, Mean, upper
data$NAb <- factor(data$NAb, levels = c("95%CI_lower", "Mean", "95%CI_upper"))

# Set position_dodge width to make points closer
position_adjust <- position_dodge(width = 0.2)

# Plot
p <- ggplot(data, aes(x = TechnicalName, y = PredictedEfficacy, group = NAb)) +
  # Plot predicted efficacy for each vaccine
  geom_point(aes(color = NAb), size = 3, position = position_adjust) +
  # Plot efficacy confidence intervals
  geom_errorbar(aes(ymin = Lower_Pred, ymax = Upper_Pred, color = NAb), 
                width = 0.1, position = position_adjust) +
  # Rotate vaccine names by 90 degrees
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18)) +
  # Style settings
  theme(panel.grid.major.x = element_blank(),  
        panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted", size = 0.3),
        panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.3),
        plot.margin = unit(c(2, 0.5, 0.5, 0.5), "cm")) +  # Margins: top, right, bottom, left
  # Set y-axis range
  scale_y_continuous(limits = c(0, 102), breaks = seq(0, 100, 20)) +
  # Adjust legend labels, modify colors and legend title
  scale_color_manual(values = c("95%CI_lower" = "#9B870C", "Mean" = "#FF69B4", "95%CI_upper" = "#006400"),
                     labels = c("95%CI_lower" = "95% CI lower", "Mean" = "Mean", "95%CI_upper" = "95% CI upper")) +
  # Add y-axis label, remove x-axis label, edit legend title
  labs(y = "Protective efficacy (%)", x = NULL, color = "Neutralization level") +
  # Adjust font and legend background
  theme(axis.title.y = element_text(size = 20),  
        axis.text = element_text(size = 18),  
        legend.position = "bottom",
        legend.text = element_text(size = 18),  
        legend.title = element_text(size = 18),  
        legend.key = element_rect(fill = "white", color = NA),  
        legend.background = element_rect(fill = "white", color = NA))

# Print figure
print(p)

# Save as high-resolution PNG file
ggsave("C:/Users/Desktop/Efficacy_Comparison_HighRes.png", plot = p, width = 8, height = 8, dpi = 1200)
