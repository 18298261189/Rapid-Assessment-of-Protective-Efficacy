# ===== Clear environment =====
rm(list = ls())

# ===== Set working directory to the path of the CSV file =====
setwd("D:\\Results\\GeneticDistance_Protection_Model_InternalValidation")

# ===== Load required packages =====
library(ggplot2)
library(dplyr)
library(readr)

# ===== Read CSV file (suppress column type messages) =====
data <- read_csv("Predicted_vs_Observed_Efficacy_InternalValidation_ImprovedGeneticDistance_Coeff2.0.csv", show_col_types = FALSE)

# Ensure Strains column is treated as character
data <- data %>% mutate(Strains = as.character(Strains))

# Adjust vaccine names: add line breaks for long names like 2CoronaVac+1ChAdOx1nCoV-19
data$TechnicalName <- gsub("\\+", "\\+\n", data$TechnicalName)

# Background grid style: dotted horizontal lines, remove vertical grid lines
grid_dotted <- element_line(color = "lightgray", linetype = "dotted")

# ===== Plot figure =====
p <- ggplot() + 
  # Observed efficacy with 95% CI
  geom_point(data = data, aes(x = TechnicalName, y = ObservedEfficacy, color = "Observed efficacy"), size = 3) +
  geom_errorbar(data = data, aes(x = TechnicalName, ymin = Lower, ymax = Upper, color = "Observed efficacy"), width = 0.1) +
  # Predicted efficacy with 95% CI
  geom_point(data = data, aes(x = TechnicalName, y = PredictedEfficacy, color = "Predicted efficacy"), size = 3, position = position_nudge(x = 0.2)) +
  geom_errorbar(data = data, aes(x = TechnicalName, ymin = Lower_Pred, ymax = Upper_Pred, color = "Predicted efficacy"), width = 0.1, position = position_nudge(x = 0.2)) +
  # Rotate vaccine names 90 degrees
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14)) +  
  # Horizontal grid style: dotted, disable minor grids
  theme(panel.grid.major.x = element_blank(),  
        panel.grid.minor = element_blank(),     
        panel.grid.major.y = element_line(color = "lightgray", linetype = "dotted", size = 0.3),  
        panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.3)) +  
  # Set y-axis range
  scale_y_continuous(limits = c(-22, 102), breaks = seq(0, 100, 20)) +
  # Add strain names at top
  annotate("text", x = 1:length(data$TechnicalName), y = 110, label = data$Strains, color = "black", angle = 0, hjust = 0.5, size = 6) +
  # Adjust margins for space at top
  theme(plot.margin = margin(t = 80, r = 20, b = 20, l = 20)) +
  # Modify legend labels and colors
  scale_color_manual(values = c("Observed efficacy" = "#00BFFF", "Predicted efficacy" = "#FF69B4")) +
  # Customize legend style
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "white", color = NA),  
        legend.key = element_rect(fill = "white", color = "transparent"),
        legend.text = element_text(size = 14)) +  
  # Axis labels and font size
  labs(y = "Efficacy (%)", x = NULL, color = NULL) +
  theme(axis.title.y = element_text(size = 16),  
        axis.text = element_text(size = 14))  

p

# ===== Save as high-resolution PNG to desktop =====
ggsave("C:/Users/Desktop/Efficacy_Comparison_HighRes.png", plot = p, width = 8, height = 6, dpi = 1200)
