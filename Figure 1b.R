rm(list=ls())
library(ggplot2)
library(numDeriv)  # Load numDeriv package to use grad function

# Set working directory and create directory for figures
setwd("D:\\Results\\Predicted_Efficacy_vs_Variant_FoldDecrease")

# Define logistic model
ProbRemainUninfected <- function(logTitre, logk, C50) {
  1 / (1 + exp(-exp(logk) * (logTitre - C50)))
}

# Logistic model function
LogisticModel_PercentUninfected <- function(mu_titre, sig_titre, logk, C50) {
  Step <- sig_titre * 0.001
  IntegralVector <- seq(mu_titre - 5 * sig_titre, mu_titre + 5 * sig_titre, by = Step)
  sum(ProbRemainUninfected(IntegralVector, logk, C50) * dnorm(IntegralVector, mu_titre, sig_titre)) * Step
}

# Set parameters for the three vaccines
std10 <- 0.4550
std <- log(10^std10)
logk <- log(2.2924 / log(10))
C50 <- log10(0.1784)  # C50 for logistic model

# Define initial efficacies for three vaccines
initial_efficacies <- c(95.40, 87.32, 61.40)
vaccine_names <- c("2Coronavac + \n1 aerosolized Ad5-nCoV", "2Coronavac + \n1 intramuscular Ad5-nCoV", "3Coronavac")

# Define fold changes (from 1 to 55 with 0.1 interval)
fold_changes <- seq(1, 55, by = 0.1)

# Create a dataframe for storing results
Results <- data.frame()

# Create a covariance matrix (based on your earlier example)
Cov <- matrix(c(0.03106460, 0.010755914, 0.01075591, 0.005727749), ncol=2)

# Loop over each vaccine
for (i in 1:length(initial_efficacies)) {
  
  # Set initial efficacy and corresponding threshold
  ef <- initial_efficacies[i]
  threshold <- qnorm(1 - ef / 100, mean = 0, sd = std)
  
  # Estimate the logistic threshold
  threshold_logistic <- -nlm(function(mu) { abs(LogisticModel_PercentUninfected(mu, std, logk, 0) - ef / 100) }, -threshold)$estimate
  
  # Calculate efficacy for each fold change
  for (fold_change in fold_changes) {
    mt_threshold_logistic <- threshold_logistic + log(fold_change)
    efficacy_variant <- LogisticModel_PercentUninfected(0, std, logk, mt_threshold_logistic) * 100
    
    # Apply consistent 95% CI calculations for the same efficacy values
    f_temp <- function(p_temp) {
      LogisticModel_PercentUninfected(0, std, p_temp[1], p_temp[2])
    }
    
    grad <- grad(f_temp, c(logk, mt_threshold_logistic))
    G <- cbind(grad[1], grad[2])
    lower_ci <- efficacy_variant - 1.96 * sqrt(G %*% Cov %*% t(G)) * 100
    upper_ci <- efficacy_variant + 1.96 * sqrt(G %*% Cov %*% t(G)) * 100
    
    # Ensure lower CI is not less than 0
    lower_ci <- max(lower_ci, 0)
    
    # Store results in a dataframe
    Results <- rbind(Results, data.frame(FoldChange = fold_change, Efficacy = efficacy_variant, Lower = lower_ci, Upper = upper_ci, Vaccine = vaccine_names[i]))
  }
}

# Set the desired order for the legend
Results$Vaccine <- factor(Results$Vaccine, levels = c("2Coronavac + \n1 aerosolized Ad5-nCoV", "2Coronavac + \n1 intramuscular Ad5-nCoV", "3Coronavac"))

# Plotting the results with adjusted 95% CI and desired legend order
ggplot(Results, aes(x = FoldChange, y = Efficacy, color = Vaccine, group = Vaccine)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Vaccine), alpha = 0.2, color = NA) +  # Remove ribbon border
  theme_classic() +
  xlab("Fold reduction in neutralizing antibody titers (relative to prototype strain)") +
  ylab("Predicted protective efficacy against variants (%)") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(1, 57, 4)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100), breaks = seq(0, 100, 10)) +
  # scale_color_manual(values = c("cyan", "orange", "green")) +  # Adjust the colors if needed
  # scale_fill_manual(values = c("cyan", "orange", "green")) +  # Adjust the fill colors if needed
  labs(color = NULL, fill = NULL) +  # Remove legend titles
  theme(
    text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),
    legend.key.height = unit(1.5, "cm"),
    plot.title = element_text(size = 20, hjust = 0.5)
  )

# Save the figure as PNG
ggsave("Vaccine_Efficacy_with_95CI.png", width = 10, height = 6)
