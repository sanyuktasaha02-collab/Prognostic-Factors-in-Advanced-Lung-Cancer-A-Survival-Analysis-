library(survival)


# Load dataset
lung <- survival::lung

str(lung)
summary(lung)
dim(lung)
head(lung)

cat("Missing values in original data:\n")
print(colSums(is.na(lung)))

cat("Status counts:\n")
print(table(lung$status))


# Prepare the data
# status = 2 means death; status = 1 means censored
lung$event <- as.integer(lung$status == 2)

# Original coding: 1 = Male, 2 = Female
lung$sex <- factor(
  lung$sex,
  levels = c(1, 2),
  labels = c("Male", "Female")
)

# Data for Kaplan-Meier analysis
km_data <- lung[
  complete.cases(lung[, c("time", "event", "sex")]),
]

# Complete cases for regression and parametric models
d <- na.omit(
  lung[, c("time", "event", "age", "sex", "ph.ecog", "wt.loss")]
)

cat("Observations for KM analysis:", nrow(km_data), "\n")
cat("Complete observations for regression:", nrow(d), "\n")
cat("Deaths in regression dataset:", sum(d$event), "\n")

stopifnot(nrow(km_data) > 0, nrow(d) > 0)


# KAPLAN-MEIER SURVIVAL CURVES

km <- survfit(
  Surv(time, event) ~ sex,
  data = km_data
)

print(km)

plot(
  km,
  col = c("red", "blue"),
  lwd = 2,
  xlab = "Time (days)",
  ylab = "Survival probability",
  main = "Kaplan-Meier Survival by Sex",
  mark.time = TRUE,
  conf.int = FALSE
)

legend(
  "topright",
  legend = c("Male", "Female"),
  col = c("red", "blue"),
  lwd = 2,
  bty = "n"
)


# LOG-RANK TEST
lr <- survdiff(
  Surv(time, event) ~ sex,
  data = km_data
)

p_logrank <- pchisq(
  lr$chisq,
  df = length(lr$n) - 1,
  lower.tail = FALSE
)

print(lr)
cat("Log-rank chi-square:", lr$chisq, "\n")
cat("Log-rank p-value:", p_logrank, "\n")

#  COX PROPORTIONAL HAZARDS MODEL

cox_model <- coxph(
  Surv(time, event) ~ age + sex + ph.ecog + wt.loss,
  data = d,
  x = TRUE
)

print(summary(cox_model))


# Hazard ratios and 95% confidence intervals
cox_results <- exp(
  cbind(
    HR = coef(cox_model),
    confint(cox_model)
  )
)

cat("Hazard ratios and 95% confidence intervals:\n")
print(cox_results)

#  CHECK PROPORTIONAL HAZARDS ASSUMPTION

ph_test <- cox.zph(cox_model)

print(ph_test)
plot(ph_test)


#  HARRELL'S C-INDEX
c_index <- summary(cox_model)$concordance

cat("C-index:", c_index[1], "\n")
cat("Standard error:", c_index[2], "\n")


#  IDENTIFY IMPORTANT PROGNOSTIC FACTORS
factor_test <- drop1(cox_model, test = "Chisq")

cat("Likelihood-ratio tests for individual predictors:\n")
print(factor_test)


# PARAMETRIC SURVIVAL MODELS
# Weibull model
weibull_model <- survreg(
  Surv(time, event) ~ age + sex + ph.ecog + wt.loss,
  data = d,
  dist = "weibull"
)

cat("Weibull model results:\n")
print(summary(weibull_model))


# Exponential model
exponential_model <- survreg(
  Surv(time, event) ~ age + sex + ph.ecog + wt.loss,
  data = d,
  dist = "exponential"
)


# Lognormal model
lognormal_model <- survreg(
  Surv(time, event) ~ age + sex + ph.ecog + wt.loss,
  data = d,
  dist = "lognormal"
)


# Compare models using AIC
aic_results <- AIC(
  exponential_model,
  weibull_model,
  lognormal_model
)

cat("Parametric model comparison using AIC:\n")
print(aic_results)


# COMPARE COX AND WEIBULL PREDICTED SURVIVAL
# Example patient profile
patient <- data.frame(
  age = 65,
  sex = factor("Female", levels = levels(d$sex)),
  ph.ecog = 1,
  wt.loss = median(d$wt.loss)
)

# Predicted survival from the Cox model
cox_curve <- survfit(
  cox_model,
  newdata = patient
)

# Time grid for Weibull prediction
times <- seq(
  0,
  max(d$time),
  by = 1
)

# Weibull accelerated failure time prediction
lp <- as.numeric(
  predict(weibull_model, newdata = patient, type = "lp")
)

scale_param <- weibull_model$scale

# Weibull survival function
weibull_survival <- exp(
  -(times / exp(lp))^(1 / scale_param)
)

# Plot both predicted curves
plot(
  cox_curve,
  col = "blue",
  lwd = 2,
  conf.int = FALSE,
  xlab = "Time (days)",
  ylab = "Predicted survival probability",
  main = "Cox vs Weibull Predicted Survival",
  xlim = c(0, max(d$time)),
  ylim = c(0, 1)
)

lines(
  times,
  weibull_survival,
  col = "red",
  lwd = 2,
  lty = 2
)

legend(
  "topright",
  legend = c("Cox model", "Weibull model"),
  col = c("blue", "red"),
  lwd = 2,
  lty = c(1, 2),
  bty = "n"
)


#  SUMMARY OF MAIN RESULTS
cat("\nPROJECT SUMMARY \n")
cat("Log-rank p-value:", signif(p_logrank, 4), "\n")
cat("C-index:", round(c_index[1], 3), "\n")
cat("Cox model observations:", nrow(d), "\n")
cat("Cox model deaths:", sum(d$event), "\n")
cat("Global PH test p-value:", 
    signif(ph_test$table["GLOBAL", "p"], 4), "\n")

cat("\nAIC comparison:\n")
print(aic_results)

