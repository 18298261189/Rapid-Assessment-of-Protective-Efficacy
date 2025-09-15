# ===== Clean environment =====
rm(list = ls())

# ===== Load required packages =====
library(lmerTest)
library(dplyr)
library(numDeriv)
library(ggplot2)
library(plotly)
library(reshape2)

# ===== Read data and fit model =====
df <- read.csv("D:/results/3Dplot/Rcode/VE-and-whole-genome-mismatch.csv",
               header = TRUE, fill = TRUE)
df$platforms <- as.factor(df$platforms)
df$VE <- df$VE2  # Use VE2 as the final VE for modeling
model <- lmer(VE ~ RBD + (1 | platforms), data = df)

# ===== Set RBD range (0 to 16) and predict efficacy at day 14 =====
RBD_seq <- seq(0, 16, length.out = 17)
platform_name <- "2CoronaVac+1AerosolizedAd5-nCoV"
pred_data <- data.frame(RBD = RBD_seq,
                        platforms = factor(platform_name, levels = levels(df$platforms)),
                        TimeAfterBoosterDose = 14)

pred_data$efficacy <- predict(model, newdata = pred_data, re.form = NULL)

# ===== Export csv1.csv (predicted VE on the line) =====
write.csv(pred_data[, c("RBD", "efficacy")],
          file = "D:/results/3Dplot/Rcode/csv1.csv", row.names = FALSE)

# ===== Read csv1 to get VE and RBD =====
csv1 <- read.csv("D:/results/3Dplot/Rcode/csv1.csv")
ef <- csv1$efficacy                         
RBD_seq <- csv1$RBD                         
hl <- rep(117.0195, length(ef))             
start_day <- rep(14, length(ef))            
labels <- rep("2CoronaVac+1AerosolizedAd5-nCoV", length(ef))

# ===== Model parameters =====
dr <- -log(2)/hl
std10 <- 0.4550
std <- log(10^std10)
k <- 2.2924/log(10)
logk <- log(k)
full_days <- seq(14, 360, by = 10)  

# ===== Define logistic functions =====
ProbRemainUninfected <- function(logTitre, logk, C50) {
  1/(1 + exp(-exp(logk) * (logTitre - C50)))
}

LogisticModel_PercentUninfected <- function(mu_titre, sig_titre, logk, C50) {
  Output <- NULL
  if (length(C50) == 1) C50 <- rep(C50, length(mu_titre))
  if (length(logk) == 1) logk <- rep(logk, length(mu_titre))
  for (i in 1:length(mu_titre)) {
    Step <- sig_titre * 0.001
    IntegralVector <- seq(mu_titre[i] - 5 * sig_titre,
                          mu_titre[i] + 5 * sig_titre, by = Step)
    Output[i] <- sum(ProbRemainUninfected(IntegralVector, logk[i], C50[i]) *
                       dnorm(IntegralVector, mu_titre[i], sig_titre)) * Step
  }
  Output
}

# ===== Back-calculate C50 =====
threshold <- qnorm(1 - ef / 100, mean = 0, sd = std)
threshold_logistic <- sapply(1:length(threshold), function(i) {
  -nlm(function(mu) abs(LogisticModel_PercentUninfected(mu, std, logk, 0) - ef[i] / 100),
       -threshold[i])$estimate
})

# ===== Main loop: calculate VE for each RBD trajectory =====
Results <- data.frame()

for (i in 1:length(ef)) {
  days <- full_days[full_days >= start_day[i]]
  mn <- dr[i] * (days - start_day[i])  
  ve <- LogisticModel_PercentUninfected(mn, std, logk, threshold_logistic[i]) * 100
  
  temp_result <- data.frame(
    day = days,
    RBD = RBD_seq[i],
    VE = ve
  )
  Results <- rbind(Results, temp_result)
}

# ===== Export 3D dataset as csv2.csv =====
write.csv(Results, file = "D:/results/3Dplot/Rcode/csv2.csv", row.names = FALSE)

# ===== Filter day = 90 and plot 2D VE vs RBD =====
subset_day <- 90
df_plot <- Results[Results$day == subset_day, ]

p <- ggplot(df_plot, aes(x = RBD, y = VE)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = paste("VE vs RBD at day", subset_day),
       x = "RBD genetic distance",
       y = "Vaccine efficacy (%)") +
  theme_minimal()

print(p)

ggsave(filename = "D:/results/3Dplot/Rcode/VE_vs_RBD.png",
       plot = p,
       width = 8, height = 6, dpi = 300)

# ===== 3D surface plot =====
df <- read.csv("D:/results/3Dplot/Rcode/csv2.csv")

df_matrix <- acast(df, day ~ RBD, value.var = "VE")

fig <- plot_ly(
  z = ~df_matrix,
  x = as.numeric(colnames(df_matrix)),   
  y = as.numeric(rownames(df_matrix)),   
  type = "surface"
)

fig <- fig %>%
  layout(
    scene = list(
      xaxis = list(title = "RBD genetic distance"),
      yaxis = list(title = "Days after booster"),
      zaxis = list(title = "Predicted efficacy (%)")
    ),
    title = "3D surface plot of predicted vaccine efficacy"
  )

# ===== Show plot =====
fig

# ===== Save interactive 3D plot as HTML (optional) =====
htmlwidgets::saveWidget(fig, "D:/results/3Dplot/Rcode/VE_surface_plot.html")
