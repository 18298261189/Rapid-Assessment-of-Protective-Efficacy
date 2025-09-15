# ===== Clear all variables in the current environment =====
rm(list = ls())

# ===== Load required package =====
library(epiR)  # For calculating CCC

# ===== Read CSV file =====
file_path <- "D:\\Results\\Predicted_vs_Observed_Efficacy_Internal+ExternalValidation_Coeff2.0.csv"
data <- read.csv(file_path)

# ===== Calculate concordance correlation coefficient (CCC) =====
ccc_result <- epi.ccc(data$ObservedEfficacy, data$PredictedEfficacy)

ccc_value <- formatC(ccc_result$rho.c[[1]], format = "f", digits = 2)
ccc_ci_lower <- formatC(ccc_result$rho.c[[2]], format = "f", digits = 2)
ccc_ci_upper <- formatC(ccc_result$rho.c[[3]], format = "f", digits = 2)

cat("Concordance correlation coefficient (CCC) =", ccc_value, "\n")
cat("95% confidence interval: [", ccc_ci_lower, ",", ccc_ci_upper, "]\n")

# ===== Format CCC and CI string for legend =====
ccc_legend <- paste("CCC =", ccc_value, "(95% CI:", ccc_ci_lower, "-", ccc_ci_upper, ")")

# ===== Plot calibration curve (high resolution, larger point size) =====
png("D:\\Results\\CalibrationPlot_Internal+ExternalValidation_Coeff2.0.png",
    width = 8, height = 8, units = "in", res = 1200)

# Set margins (bottom, left, top, right)
par(mar = c(5, 6, 2, 2))

# Plot observed vs predicted efficacy
plot(
  data$ObservedEfficacy, data$PredictedEfficacy,
  xlab = "Observed efficacy (%)",
  ylab = "Predicted efficacy (%)",
  xlim = c(0, 100),
  ylim = c(0, 100),
  pch = 16,   # Solid circles
  col = "skyblue",
  main = " ", # No main title
  cex = 4.0,      # Point size
  cex.lab = 1.8,  # Axis label size
  cex.axis = 1.5  # Axis tick label size
)

# Add diagonal line
abline(0, 1, lty = 2, col = "black")

# Add CCC value with CI as legend (top left)
legend("topleft", legend = ccc_legend, bty = "n", cex = 1.5)

# Save figure
dev.off()

cat("Concordance correlation coefficient (CCC) =", ccc_value, "\n")
cat("95% confidence interval: [", ccc_ci_lower, ",", ccc_ci_upper, "]\n")
