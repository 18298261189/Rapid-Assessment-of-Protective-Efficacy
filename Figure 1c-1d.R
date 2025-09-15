# Clear environment
rm(list = ls())

# Set working directory to the path containing CSV file
setwd("D:\\Results\\NeutralizingAntibody_Protection_InternalExternalValidation_Calibration")

# Load required packages
library(ggplot2)
library(dplyr)
library(readr)

# Read CSV file and suppress column type messages
data <- read_csv("Predicted_vs_Observed_Efficacy_ExternalValidation.csv", show_col_types = FALSE)

# Ensure Strains column is text
data <- data %>% mutate(Strains = as.character(Strains))

# Adjust vaccine names, handle line break for 2CoronaVac+1ChAdOx1nCoV-19
data$TechnicalName <- gsub("\\+", "\\+\n", data$TechnicalName)

# Define grid style with dotted horizontal lines, remove vertical grid lines
grid_dotted <- element_line(color = "lightgray", linetype = "dotted")

# Plot figure
p <- ggplot() + 
  # Observed efficacy with 95% CI
  geom_point(data = data, aes(x = TechnicalName, y = ObservedEfficacy, color = "Reported efficacy"), size = 3) +
  geom_errorbar(data = data, aes(x = TechnicalName, ymin = Lower, ymax = Upper, color = "Reported efficacy"), width = 0.1) +
  # Predicted efficacy with 95% CI
  geom_point(data = data, aes(x = TechnicalName, y = PredictedEfficacy, color = "Predicted efficacy"), size = 3, position = position_nudge(x = 0.2)) +
  geom_errorbar(data = data, aes(x = TechnicalName, ymin = Lower_Pred, ymax = Upper_Pred, color = "Predicted efficacy"), width = 0.1, position = position_nudge(x = 0.2)) +
  # Rotate vaccine names 90 degrees and adjust position
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14)) +  
  # Set horizontal grid line style as dotted, disable minor grid lines
  theme(panel.grid.major.x = element_blank(),  
        panel.grid.minor = element_blank(),     
        panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted", size = 0.3),  
        panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.3)) +  
  # Set y-axis range 0–102
  scale_y_continuous(limits = c(0, 102), breaks = seq(0, 100, 20)) +
  # Add strain names at the top
  annotate("text", x = 1:length(data$TechnicalName), y = 110, label = data$Strains, color = "black", angle = 0, hjust = 0.5, size = 14) +  
  # Adjust margins to make space for strain labels
  theme(plot.margin = margin(t = 80, r = 20, b = 20, l = 20)) +  
  # Customize legend labels and set colors
  scale_color_manual(values = c("Reported efficacy" = "#00BFFF", "Predicted efficacy" = "#FF69B4")) +
  # Customize legend style, white background, remove title, increase font size
  theme(legend.position = "right",
        legend.background = element_rect(fill = "white", color = NA),  
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 14)) +  
  # Set y-axis label, remove x-axis label, enlarge font
  labs(y = "Protective efficacy (%)", x = NULL, color = NULL) +
  # Fully customize legend at bottom, horizontal orientation
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.background = element_blank(),
    legend.key = element_rect(fill = "white", color = "white"),
    legend.text = element_text(size = 14),
    axis.text.y = element_text(size = 14)  # set y-axis tick font size
  )

p

# Save as high-resolution PNG to desktop (dpi = 1200)
ggsave("C:/Users/Desktop/Efficacy_Comparison_HighRes.png", plot = p, width = 8, height = 6, dpi = 1200)
