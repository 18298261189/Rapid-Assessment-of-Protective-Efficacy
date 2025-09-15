# Clear all variables from current environment
rm(list = ls())

# Load required package
library(epiR)  # For calculating CCC (concordance correlation coefficient)

# Read CSV file
file_path <- "D:\\Results\\NeutralizingAntibody_Protection_InternalExternalValidation_Calibration\\Predicted_vs_Observed_Efficacy_InternalExternal.csv"
data <- read.csv(file_path)

# # Multiply ObservedEfficacy and PredictedEfficacy columns by 100 (if needed)
# data$ObservedEfficacy <- data$ObservedEfficacy * 100
# data$PredictedEfficacy <- data$PredictedEfficacy * 100

# Calculate concordance correlation coefficient (CCC)
ccc_result <- epi.ccc(data$ObservedEfficacy, data$PredictedEfficacy)
ccc_value <- round(ccc_result$rho.c[1], 2)      # Extract and round CCC value
ccc_ci_lower <- round(ccc_result$rho.c[2], 2)   # Extract and round 95% CI lower bound
ccc_ci_upper <- round(ccc_result$rho.c[3], 2)   # Extract and round 95% CI upper bound

# Format CCC and 95% CI as legend string
ccc_legend <- paste("CCC =", ccc_value, "(95% CI:", ccc_ci_lower, "-", ccc_ci_upper, ")")

# Draw calibration plot with higher resolution and larger point size
png("D:\\JinLairun\\PhD\\Supervisor_LiJingxin\\Thesis\\Main_Manuscript\\Figures_and_Rcode\\NeutralizingAntibody_Protection_InternalExternalValidation_Calibration\\CalibrationPlot_Internal_External.png", 
    width = 8, height = 8, units = "in", res = 1200)

# Set margins, increase left margin
par(mar = c(5, 6, 2, 2))  # bottom, left, top, right

# Set plotting parameters, increase axis and label font size
plot(
  data$ObservedEfficacy, data$PredictedEfficacy, 
  xlab = "Reported protective efficacy (%)", 
  ylab = "Predicted protective efficacy (%)", 
  xlim = c(0, 100), 
  ylim = c(0, 100), 
  pch = 16,       # Solid circle points
  col = "skyblue", 
  main = " ",     # No main title
  cex = 4.0,      # Point size
  cex.lab = 1.8,  # Axis label font size
  cex.axis = 1.5  # Axis tick font size
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "black")

# Add CCC value with 95% CI as legend in top-left corner
legend("topleft", legend = ccc_legend, bty = "n", cex = 1.5)

# Save the figure
dev.off()
