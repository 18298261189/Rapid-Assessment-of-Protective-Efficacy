# ===== Clear environment =====
rm(list=ls())

# ===== Load required packages =====
library(ggplot2)
library(numDeriv)  

# ===== Set working directory and create folder =====
setwd("D:/results/VE_persistence_by_strains")
dir.create("figures", showWarnings = FALSE)

# Function to tag figure files with date
myfun_date <- function() format(Sys.time(), "%Y %m %d %H%M%S")
fig_height <- 6
fig_width <- 12

# ===== Initial efficacy and half-lives for vaccines =====
ef <- c(95.40, 87.32, 61.40)          # Initial efficacy
hl <- c(117.0195, 65.1616, 73.5492)   # Antibody half-lives
start_day <- c(14, 14, 14)            # Starting day for each vaccine

# ===== Calculate decay rates =====
dr <- -log(2)/hl  
std10 <- 0.4550                       # Antibody SD (log10 scale)
std <- log(10^std10)                  # SD in natural log
k <- 2.2924/log(10)                   # Logistic k value
logk <- log(k)
max_days <- 360
full_days <- 0:max_days               # Day 0 to 360

# ===== Logistic functions =====
ProbRemainUninfected <- function(logTitre, logk, C50) {
  1/(1 + exp(-exp(logk) * (logTitre - C50)))
}

LogisticModel_PercentUninfected <- function(mu_titre, sig_titre, logk, C50) {
  Output <- NULL
  if (length(C50) == 1) C50 <- rep(C50, length(mu_titre))
  if (length(logk) == 1) logk <- rep(logk, length(mu_titre))
  for (i in 1:length(mu_titre)) {
    Step <- sig_titre[i] * 0.001
    IntegralVector <- seq(mu_titre[i] - 5 * sig_titre[i], mu_titre[i] + 5 * sig_titre[i], by = Step)
    Output[i] <- sum(ProbRemainUninfected(IntegralVector, logk[i], C50[i]) *
                       dnorm(IntegralVector, mu_titre[i], sig_titre[i])) * Step
  }
  Output
}

# ===== Initial thresholds =====
threshold <- qnorm(1 - ef / 100, mean = 0, sd = std)

# ===== Back-calculate logistic thresholds =====
threshold_logistic <- NULL
for (i in 1:length(threshold)) {
  threshold_logistic[i] <- -nlm(
    function(mu) { abs(LogisticModel_PercentUninfected(mu, std, logk, 0) - ef[i] / 100) },
    -threshold[i]
  )$estimate
}

# ===== Vaccine labels =====
labels <- c("2CoronaVac+\n1 Aerosolized Ad5-nCoV",
            "2CoronaVac+\n1 Intramuscular Ad5-nCoV",
            "3CoronaVac")

# ===== Compute efficacy over time (with 95% CI) =====
Results <- NULL
Lower_Pred <- NULL
Upper_Pred <- NULL

for (i in 1:length(threshold)) {
  mn <- rep(NA, length(full_days))  
  days <- full_days[full_days >= start_day[i]]  
  mn[full_days >= start_day[i]] <- 0 + dr[i] * (days - start_day[i])  
  
  efficacy2 <- rep(NA, length(full_days))
  lower <- rep(NA, length(full_days))
  upper <- rep(NA, length(full_days))
  
  Cov <- matrix(c(0.04, 0, 0, 0.01), nrow = 2)  # Example covariance matrix
  for (j in which(full_days >= start_day[i])) {
    efficacy2[j] <- LogisticModel_PercentUninfected(mn[j], std, logk, threshold_logistic[i]) * 100
    
    grad1 <- grad(function(p) LogisticModel_PercentUninfected(mn[j], std, p[1], p[2]),
                  c(logk, threshold_logistic[i]))[1]
    grad2 <- grad(function(p) LogisticModel_PercentUninfected(mn[j], std, p[1], p[2]),
                  c(logk, threshold_logistic[i]))[2]
    G <- cbind(grad1, grad2)
    
    Lower_Pred[j] <- max(efficacy2[j] - 1.96 * sqrt(G %*% Cov %*% t(G)) * 100, 0)
    Upper_Pred[j] <- efficacy2[j] + 1.96 * sqrt(G %*% Cov %*% t(G)) * 100
  }
  
  temp_results <- data.frame(Efficacy_logistic = efficacy2, 
                             Lower_CI = Lower_Pred, 
                             Upper_CI = Upper_Pred,
                             initialEfficacy = ef[i], 
                             day = full_days, 
                             Vaccine = labels[i])
  
  temp_results <- temp_results[temp_results$day >= start_day[i], ]
  Results <- rbind(Results, temp_results)
}

# ===== Factor order for plotting =====
Results$Vaccine <- factor(Results$Vaccine, 
                          levels = c("2CoronaVac+\n1 Aerosolized Ad5-nCoV",
                                     "2CoronaVac+\n1 Intramuscular Ad5-nCoV",
                                     "3CoronaVac"))

# ===== Plot and save figure =====
ggplot(Results, aes(x = day, y = Efficacy_logistic, color = Vaccine, group = Vaccine)) + 
  geom_line() + 
  theme_classic() +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI, fill = Vaccine), alpha = 0.2, linetype = 0) +
  geom_hline(yintercept = ef, linetype = "dotted", size = 0.1) +
  xlab("Days after booster") +
  ylab("Predicted efficacy (%)") +
  scale_x_continuous(expand = c(0, 0), breaks = c(14, 90, 180, 270, 360), limits = c(0, 360))+
  scale_y_continuous(expand = c(0, 1), breaks = seq(0, 100, 10), limits = c(0, 105)) +
  labs(color = "", fill = "") +
  theme(
    legend.key.height = unit(1.0, "cm"),
    axis.line = element_line(size = 0.2),
    axis.ticks = element_line(size = 0.2),
    axis.text.x = element_text(size = 6),
    axis.ticks.length = unit(0.05, "cm"),
    axis.text.y = element_text(size = 6)
  )

ggsave(
  filename = paste0("figures/Figure2A_logistic_with_CI_", myfun_date(), ".png"),
  width = fig_width,
  height = fig_height,
  units = "cm",
  dpi = 1200
)
