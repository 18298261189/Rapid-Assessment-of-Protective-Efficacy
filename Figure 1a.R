rm(list=ls())

library(ggplot2)

# Define ProbRemainUninfected function
ProbRemainUninfected <- function(logTitre, logk, C50) {
  return(1 / (1 + exp(-exp(logk) * (logTitre - C50))))
}

# Define LogisticModel_PercentUninfected function
LogisticModel_PercentUninfected <- function(mu_titre, sig_titre, logk, C50) {
  Step <- sig_titre * 0.001
  IntegralVector <- seq(mu_titre - 5 * sig_titre, mu_titre + 5 * sig_titre, by=Step)
  return(sum(ProbRemainUninfected(IntegralVector, logk, C50) * dnorm(IntegralVector, mu_titre, sig_titre)) * Step)
}

# Set parameters
PooledSD <- 0.4550
logk <- log(2.2924)
C50 <- log10(0.1784)
Cov <- matrix(c(0.03106460, 0.010755914, 0.01075591, 0.005727749), ncol=2)

# Define three groups of NeutMean and NeutConv
NeutMean1 <- 1772
NeutConv1 <- 320
NeutMean2 <- 234
NeutConv2 <- 320
NeutMean3 <- 42
NeutConv3 <- 320

# Calculate three groups of NeutRatio_Reported
NeutRatio_Reported1 <- log10(NeutMean1 / NeutConv1)
NeutRatio_Reported2 <- log10(NeutMean2 / NeutConv2)
NeutRatio_Reported3 <- log10(NeutMean3 / NeutConv3)

# Calculate three groups of Covaxin_PointEstimate
Covaxin_PointEstimate1 <- LogisticModel_PercentUninfected(NeutRatio_Reported1, PooledSD, logk, C50)
Covaxin_PointEstimate2 <- LogisticModel_PercentUninfected(NeutRatio_Reported2, PooledSD, logk, C50)
Covaxin_PointEstimate3 <- LogisticModel_PercentUninfected(NeutRatio_Reported3, PooledSD, logk, C50)

# Compute gradients
numericGradient <- function(func, params){
  epsilon <- 1e-5
  grads <- numeric(length(params))
  for(i in 1:length(params)){
    params1 <- params
    params2 <- params
    params1[i] <- params1[i] + epsilon
    params2[i] <- params2[i] - epsilon
    grads[i] <- (func(params1) - func(params2)) / (2 * epsilon)
  }
  return(grads)
}

covaxin_temp <- function(p_temp, neut_ratio_reported) {
  LogisticModel_PercentUninfected(neut_ratio_reported, PooledSD, p_temp[1], p_temp[2])
}

# Calculate gradients and confidence intervals for three groups
grad_covaxin1 <- numericGradient(function(p_temp) covaxin_temp(p_temp, NeutRatio_Reported1), c(logk, C50))
grad_covaxin2 <- numericGradient(function(p_temp) covaxin_temp(p_temp, NeutRatio_Reported2), c(logk, C50))
grad_covaxin3 <- numericGradient(function(p_temp) covaxin_temp(p_temp, NeutRatio_Reported3), c(logk, C50))

Covaxin_LowerB1 <- Covaxin_PointEstimate1 - 1.96 * sqrt(sum(grad_covaxin1 %*% Cov %*% grad_covaxin1))
Covaxin_UpperB1 <- Covaxin_PointEstimate1 + 1.96 * sqrt(sum(grad_covaxin1 %*% Cov %*% grad_covaxin1))

Covaxin_LowerB2 <- Covaxin_PointEstimate2 - 1.96 * sqrt(sum(grad_covaxin2 %*% Cov %*% grad_covaxin2))
Covaxin_UpperB2 <- Covaxin_PointEstimate2 + 1.96 * sqrt(sum(grad_covaxin2 %*% Cov %*% grad_covaxin2))

Covaxin_LowerB3 <- Covaxin_PointEstimate3 - 1.96 * sqrt(sum(grad_covaxin3 %*% Cov %*% grad_covaxin3))
Covaxin_UpperB3 <- Covaxin_PointEstimate3 + 1.96 * sqrt(sum(grad_covaxin3 %*% Cov %*% grad_covaxin3))

# Print predicted efficacy and 95% CI
print(paste("2CoronaVac+1 aerosolized Ad5-nCoV: ", round(100 * Covaxin_PointEstimate1, 2), "% [", round(100 * Covaxin_LowerB1, 2), "%, ", round(100 * Covaxin_UpperB1, 2), "%]"))
print(paste("2CoronaVac+1 intramuscular Ad5-nCoV: ", round(100 * Covaxin_PointEstimate2, 2), "% [", round(100 * Covaxin_LowerB2, 2), "%, ", round(100 * Covaxin_UpperB2, 2), "%]"))
print(paste("3CoronaVac: ", round(100 * Covaxin_PointEstimate3, 2), "% [", round(100 * Covaxin_LowerB3, 2), "%, ", round(100 * Covaxin_UpperB3, 2), "%]"))

# Create Covaxin table including three vaccines
CovaxinTable <- data.frame(
  NeutRatio_Reported = c(NeutRatio_Reported1, NeutRatio_Reported2, NeutRatio_Reported3),
  Efficacy = c(100 * Covaxin_PointEstimate1, 100 * Covaxin_PointEstimate2, 100 * Covaxin_PointEstimate3),
  Lower = c(100 * Covaxin_LowerB1, 100 * Covaxin_LowerB2, 100 * Covaxin_LowerB3),
  Upper = c(100 * Covaxin_UpperB1, 100 * Covaxin_UpperB2, 100 * Covaxin_UpperB3),
  TechnicalName = c("2CoronaVac+1 aerosolized Ad5-nCoV", "2CoronaVac+1 intramuscular Ad5-nCoV", "3CoronaVac"),
  VaccineColor = c("darkcyan", "darkorange", "darkgreen")
)

# Define neutralization values
NeutValue <- seq(0.1, 11, by = 0.001)

# Calculate efficacy under Logistic model
Efficacy_Logistic_Raw <- sapply(log10(NeutValue), function(x) {
  LogisticModel_PercentUninfected(x, PooledSD, logk, C50)
})

# Calculate CI for Logistic model
grad1 <- numeric(length(NeutValue))
grad2 <- numeric(length(NeutValue))
Lower_Pred <- numeric(length(NeutValue))
Upper_Pred <- numeric(length(NeutValue))

for (i in 1:length(NeutValue)) {
  f_temp <- function(p_temp) {
    LogisticModel_PercentUninfected(log10(NeutValue[i]), PooledSD, p_temp[1], p_temp[2])
  }
  grad1[i] <- numericGradient(f_temp, c(logk, C50))[1]
  grad2[i] <- numericGradient(f_temp, c(logk, C50))[2]
  G <- cbind(grad1[i], grad2[i])
  Lower_Pred[i] <- Efficacy_Logistic_Raw[i] - 1.96 * sqrt(G %*% Cov %*% t(G))
  Upper_Pred[i] <- Efficacy_Logistic_Raw[i] + 1.96 * sqrt(G %*% Cov %*% t(G))
}

# Create data frame for plotting
LogisticModel_withPoolSD <- data.frame(
  "NeutRatio_Reported" = log10(NeutValue),
  "Efficacy" = 100 * Efficacy_Logistic_Raw,
  "Lower" = 100 * Lower_Pred,
  "Upper" = 100 * Upper_Pred
)

CovaxinTable$TechnicalName <- factor(CovaxinTable$TechnicalName, 
                                     levels = c("2CoronaVac+1 aerosolized Ad5-nCoV", 
                                                "2CoronaVac+1 intramuscular Ad5-nCoV", 
                                                "3CoronaVac"))

### Add previously published baseline immunity efficacy prediction curve
# Define group B parameters
PooledSD_B <- 0.4647
logk_B <- log(3.0977)
C50_B <- log10(0.2011)

# Calculate efficacy for group B
Efficacy_Logistic_Raw_B <- sapply(log10(NeutValue), function(x) {
  LogisticModel_PercentUninfected(x, PooledSD_B, logk_B, C50_B)
})

# Calculate CI for group B
grad1_B <- numeric(length(NeutValue))
grad2_B <- numeric(length(NeutValue))
Lower_Pred_B <- numeric(length(NeutValue))
Upper_Pred_B <- numeric(length(NeutValue))

for (i in 1:length(NeutValue)) {
  f_temp_B <- function(p_temp) {
    LogisticModel_PercentUninfected(log10(NeutValue[i]), PooledSD_B, p_temp[1], p_temp[2])
  }
  grad1_B[i] <- numericGradient(f_temp_B, c(logk_B, C50_B))[1]
  grad2_B[i] <- numericGradient(f_temp_B, c(logk_B, C50_B))[2]
  G_B <- cbind(grad1_B[i], grad2_B[i])
  Lower_Pred_B[i] <- Efficacy_Logistic_Raw_B[i] - 1.96 * sqrt(G_B %*% Cov %*% t(G_B))
  Upper_Pred_B[i] <- Efficacy_Logistic_Raw_B[i] + 1.96 * sqrt(G_B %*% Cov %*% t(G_B))
}

# Create data frame for group B curve
LogisticModel_withPoolSD_B <- data.frame(
  "NeutRatio_Reported" = log10(NeutValue),
  "Efficacy" = 100 * Efficacy_Logistic_Raw_B,
  "Lower" = 100 * Lower_Pred_B,
  "Upper" = 100 * Upper_Pred_B
)

Figure11_with_B <- ggplot() +
  # Group A 95% CI band (optimized)
  geom_ribbon(data = LogisticModel_withPoolSD, aes(x = 10^NeutRatio_Reported, ymin = Lower, ymax = Upper, fill = "Optimized"), alpha = 0.3) +
  
  # Group B 95% CI band (before optimization)
  geom_ribbon(data = LogisticModel_withPoolSD_B, aes(x = 10^NeutRatio_Reported, ymin = Lower, ymax = Upper, fill = "Before optimization"), alpha = 0.3) +
  
  # Group A fitted curve (optimized)
  geom_line(data = LogisticModel_withPoolSD, aes(x = 10^NeutRatio_Reported, y = Efficacy, color = "Optimized"), size = 0.5) +
  
  # Group B fitted curve (before optimization)
  geom_line(data = LogisticModel_withPoolSD_B, aes(x = 10^NeutRatio_Reported, y = Efficacy, color = "Before optimization"), size = 0.5) +
  
  # Set x axis to log scale and breaks
  scale_x_log10(lim = c(0.1, 11), breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8), labels = c(0.125, 0.25, 0.5, 1, 2, 4, 8)) +
  scale_y_continuous(lim = c(0, 100), breaks = seq(0, 100, by = 10)) +
  
  # Axis labels
  xlab("Mean neutralization level (normalized to convalescent sera)") +
  ylab("Predicted protective efficacy") +
  
  # Theme adjustments
  theme_linedraw() +
  theme(
    panel.grid.major = element_line(linetype = "99", color = "grey30"),
    panel.grid.minor = element_line(linetype = "99", color = "grey30"),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linetype = "solid"),
    
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  
  # Manual colors and merge legends
  scale_fill_manual(values = c("Optimized" = "lightcoral", "Before optimization" = "darkgray"), 
                    labels = c("Optimized", "Before optimization")) +
  scale_color_manual(values = c("Optimized" = "red", "Before optimization" = "darkgray"), 
                     labels = c("Optimized", "Before optimization")) +
  
  guides(fill = guide_legend(order = 1, title = NULL),
         color = guide_legend(order = 1, title = NULL))

# Save file path
file_path_B <- "C:/Users/Desktop/Figure11_with_AB_Optimized.png"

# Save figure
ggsave(file_path_B, plot = Figure11_with_B, width = 7, height = 5, dpi = 1200)

# Show figure
print(Figure11_with_B)
