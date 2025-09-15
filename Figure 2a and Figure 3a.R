# ===== Clear environment =====
rm(list = ls())

# ===== Load required packages =====
library(ggplot2)
library(lmerTest)
library(performance)
library(dplyr)

# ===== Read data =====
df <- read.csv("D:\\Results\\VE-and-whole-genome-mismatch.csv",
               header = TRUE, fill = TRUE)

df$platforms <- as.factor(df$platforms)
df$VE <- df$VE2  # Use VE2 column as final VE for modeling

# ===== Mixed-effects model fitting =====
Model_combined <- lmer(VE ~ RBD + (1 | platforms), data = df)
summary(Model_combined)

# Extract random effect standard deviations
re_sd <- as.data.frame(VarCorr(Model_combined))

# Extract platform and residual standard deviation
platform_sd <- re_sd$sdcor[re_sd$grp == "platforms"]
resid_sd <- re_sd$sdcor[re_sd$grp == "Residual"]

# Display with 4 decimal places
cat("Standard deviation of platform intercepts:", format(platform_sd, digits = 4, nsmall = 4), "\n")
cat("Standard deviation of residuals          :", format(resid_sd, digits = 4, nsmall = 4), "\n")

# ===== Show fixed-effect coefficients, keep t values with 4 decimals =====
coefs <- summary(Model_combined)$coefficients
coefs_df <- as.data.frame(coefs)
print(coefs_df)

# ===== Show R² (keep 4 decimals) =====
r2_values <- performance::r2(Model_combined)
cat("Marginal R² (fixed effects):", format(round(r2_values$R2_marginal, 4), nsmall = 4), "\n")
cat("Conditional R² (fixed + random effects):", format(round(r2_values$R2_conditional, 4), nsmall = 4), "\n")

# ===== Standardized regression coefficient =====
# Method: standardize numeric predictor RBD, fit model with same random structure
df_z <- df
df_z$VE_z  <- as.numeric((df$VE))
df_z$RBD_z <- as.numeric(scale(df$RBD))

Model_std <- lmer(VE_z ~ RBD_z + (1 | platforms), data = df_z)

cat("\n=== Standardized model (dependent variable and RBD standardized) ===\n")
std_coefs <- summary(Model_std)$coefficients
std_coefs_df <- as.data.frame(std_coefs)

print(std_coefs_df)

# ===== Create prediction dataset (TimeAfterBoosterDose fixed at 90) =====
pred_data <- expand.grid(
  RBD = seq(0, 30, length.out = 301),
  platforms = factor(levels(df$platforms), levels = levels(df$platforms)),
  TimeAfterBoosterDose = 90
)

# ===== Model predictions and SE calculation =====
pred_data$VE <- predict(Model_combined, newdata = pred_data, re.form = NULL)
X <- model.matrix(formula(Model_combined, fixed.only = TRUE), data = pred_data)
beta_cov <- vcov(Model_combined)
pred_data$se_unweighted <- sqrt(rowSums((X %*% beta_cov) * X))

# ===== Sample size info (not used for SE) =====
samplesize_sum <- df %>%
  group_by(platforms) %>%
  summarise(total_n = sum(samplesize, na.rm = TRUE)) %>%
  ungroup()
pred_data <- pred_data %>% left_join(samplesize_sum, by = "platforms")

# ===== Build confidence intervals =====
pred_data <- pred_data %>%
  mutate(
    se = se_unweighted,
    lower_CI = VE - 1.96 * se,
    upper_CI = VE + 1.96 * se
  )

# ===== Platform color settings =====
platform_colors <- c(
  "2CoronaVac+1AerosolizedAd5-nCoV" = "#E69F00",
  "2CoronaVac+1Ad5-nCoV"            = "#56B4E9",
  "3CoronaVac"                      = "#009E73",
  "2CoronaVac+1BNT162b2"            = "#F0E442",
  "2CoronaVac+1ChAdOx1"             = "#0072B2"
)

# ===== Select three regimens of interest =====
platforms_to_keep <- c("2CoronaVac+1AerosolizedAd5-nCoV", 
                       "2CoronaVac+1Ad5-nCoV", 
                       "3CoronaVac")

# ===== Plot overall prediction =====
pp_total <- ggplot() +
  geom_ribbon(data = pred_data, aes(x = RBD, ymin = lower_CI, ymax = upper_CI, fill = platforms), alpha = 0.2) +
  geom_line(data = pred_data, aes(x = RBD, y = VE, color = platforms), size = 1.2) +
  geom_point(data = df, aes(x = RBD, y = VE, color = platforms, shape = platforms),
             position = position_jitter(width = 0, height = 0),
             size = 10, alpha = 0.7) +
  labs(
    x = "RBD genetic distance",
    y = "Predicted protective efficacy (%)",
    title = "Predicted protective efficacy across regimens",
    color = "Vaccination regimen",
    fill = "Vaccination regimen",
    shape = "Vaccination regimen"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 30),
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30)
  ) +
  scale_color_manual(values = platform_colors) +
  scale_fill_manual(values = platform_colors) +
  scale_shape_manual(values = rep(16, 5)) +
  scale_y_continuous(limits = c(-22, 120), breaks = seq(0, 100, by = 20)) +
  scale_x_continuous(limits = c(0, 16), breaks = seq(0, 16, by = 2))

print(pp_total)

# ===== Subplots (TimeAfterBoosterDose = 90) =====
for (platform_name in levels(df$platforms)) {
  df_sub <- df %>% filter(platforms == platform_name)
  
  pred_data_sub <- expand.grid(
    RBD = seq(0, 30, length.out = 301),
    platforms = factor(platform_name, levels = levels(df$platforms)),
    TimeAfterBoosterDose = 90
  )
  
  pred_data_sub$VE <- predict(Model_combined, newdata = pred_data_sub, re.form = NULL)
  X_sub <- model.matrix(formula(Model_combined, fixed.only = TRUE), data = pred_data_sub)
  pred_data_sub$se <- sqrt(rowSums((X_sub %*% beta_cov) * X_sub))
  
  pred_data_sub <- pred_data_sub %>%
    mutate(
      lower_CI = VE - 1.96 * se,
      upper_CI = VE + 1.96 * se
    )
  
  color_this <- platform_colors[as.character(platform_name)]
  
  pp_sub <- ggplot() +
    geom_ribbon(data = pred_data_sub,
                aes(x = RBD, ymin = lower_CI, ymax = upper_CI),
                fill = color_this, alpha = 0.2) +
    geom_line(data = pred_data_sub,
              aes(x = RBD, y = VE),
              color = color_this, size = 1.3) +
    geom_point(data = df_sub,
               aes(x = RBD, y = VE),
               size = 10, alpha = 0.7, shape = 16,
               color = color_this) +
    labs(
      x = "RBD genetic distance",
      y = "Predicted protective efficacy (%)",
      title = platform_name
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      text = element_text(size = 30),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 30)
    ) +
    scale_y_continuous(limits = c(-22, 120), breaks = seq(0, 100, by = 20)) +
    scale_x_continuous(limits = c(0, 18), breaks = seq(0, 18, by = 3))
  
  print(pp_sub)
  
  ggsave(filename = paste0("D:\\Results\\Subplot_", platform_name, "_95CI_5Regimens_Coeff2.0.png"),
         plot = pp_sub, width = 10, height = 8, dpi = 1200)
}






###Figure 3a

## Prediction
# ===== Prediction at a specified RBD genetic distance (time fixed at 90 days) =====
RBD_input <- 17

pred_input <- data.frame(
  RBD = RBD_input,
  platforms = factor(levels(df$platforms), levels = levels(df$platforms)),
  TimeAfterBoosterDose = 90
)

# Predicted VE
pred_input$VE <- predict(Model_combined, newdata = pred_input, re.form = NULL)

# Standard error of prediction
X_input <- model.matrix(formula(Model_combined, fixed.only = TRUE), data = pred_input)
se_input <- sqrt(rowSums((X_input %*% beta_cov) * X_input))

# 95% CI
pred_input$lower_CI <- pred_input$VE - 1.96 * se_input
pred_input$upper_CI <- pred_input$VE + 1.96 * se_input

# Format output
pred_input_fmt <- pred_input %>%
  mutate(across(c(VE, lower_CI, upper_CI), ~ formatC(.x, format = "f", digits = 2))) %>%
  rename(
    Platform = platforms,
    GeneticDistance = RBD,
    TimePoint = TimeAfterBoosterDose,
    PredictedVE = VE,
    Lower95CI = lower_CI,
    Upper95CI = upper_CI
  )

cat("At genetic distance RBD =", RBD_input, "and observation time = 90 days, the predicted vaccine efficacy for each regimen is:\n")
print(pred_input_fmt, row.names = FALSE)

##### Note: TimeAfterBoosterDose here represents the observation time point.
