PROGNOSTIC FACTORS IN ADVANCED LUNG CANCER: A SURVIVAL ANALYSIS
================================================================
 
Author  : Sanyukta Saha
Tools   : R (survival package)
Dataset : survival::lung (NCCTG advanced lung cancer study)
 
 
1. OVERVIEW
-----------
This project analyses how long patients with advanced lung cancer survive
and which patient characteristics are linked to shorter or longer survival.
Because the outcome is time to death and some patients are still alive at the
end of follow-up (censored), survival analysis methods are used instead of
ordinary regression.
 
 
2. OBJECTIVES
-------------
- Describe the data and handle missing values
- Compare survival curves for men and women
- Estimate the independent effect of age, sex, performance status and
  weight loss on the risk of death
- Check the assumptions and predictive ability of the Cox model
- Fit parametric survival models and compare them with the Cox model
 
 
3. DATA
-------
228 patients, 10 variables; survival time in days; 165 deaths (72.4%) and
63 censored (27.6%).
 
Variables used in the models:
  time     Survival time in days
  status   1 = censored, 2 = dead (recoded to event: 1 = death, 0 = censored)
  age      Age in years
  sex      1 = male, 2 = female (recoded as a factor)
  ph.ecog  ECOG performance score rated by the physician (0 = best, higher = worse)
  wt.loss  Weight loss in the last six months (pounds)
 
Not used: inst, ph.karno, pat.karno (overlap with ph.ecog) and meal.cal
(47 missing values).
 
Missing data: 15 patients (6.6%) had missing model variables and were removed,
leaving 213 patients and 151 deaths for the regression and parametric models.
 
 
4. METHODS
----------
1. Data inspection and missing-value check
2. Kaplan-Meier survival curves by sex
3. Log-rank test to compare the two curves
4. Cox proportional hazards model (age + sex + ph.ecog + wt.loss)
5. Proportional hazards check (cox.zph and Schoenfeld residual plots)
6. Concordance index (C-index) to measure discrimination
7. Likelihood-ratio tests for each predictor (drop1)
8. Parametric models: Weibull, exponential and log-normal, compared by AIC
9. Predicted survival for an example patient: Cox vs Weibull
 
 
5. KEY RESULTS
--------------
Kaplan-Meier and log-rank
  - Median survival: 426 days (female) vs 270 days (male)
  - Log-rank chi-square = 10.3, 1 df, p = 0.0013
 
Cox model (n = 213, 151 deaths)
  - Female sex : HR 0.554 (95% CI 0.393 to 0.781), p = 0.0008
  - ph.ecog    : HR 1.674 per point (95% CI 1.308 to 2.143), p < 0.0001
  - Age        : HR 1.013 (95% CI 0.995 to 1.033), p = 0.165
  - Weight loss: HR 0.991 (95% CI 0.978 to 1.004), p = 0.176
 
Model checks
  - Proportional hazards holds for all variables (global p = 0.314)
  - C-index = 0.647 (SE 0.026), modest discrimination
  - Removing sex or ph.ecog worsens the fit; removing age or weight loss does not
 
Parametric models
  - AIC: Weibull 2106.98, exponential 2128.13, log-normal 2136.39
  - Weibull is best; scale 0.716 (shape about 1.40, hazard rises over time)
  - Weibull time ratios: female 1.52, ph.ecog 0.69 per point
  - Converted Weibull hazard ratios (0.56 and 1.67) match the Cox results
 
 
6. CONCLUSIONS
--------------
- Women survived longer than men, and the difference remains after adjusting
  for age, performance status and weight loss.
- Poorer performance status (ph.ecog) is the strongest adverse factor.
- Age and weight loss added little once sex and ph.ecog were known.
- The Cox and Weibull models agree, giving consistent evidence.
 
 
7. HOW TO RUN
-------------
Install the package once:  install.packages("survival")
Then run the analysis script from top to bottom in R. It loads survival::lung,
recodes the variables, fits all models and produces the plots.
 
 
8. LIMITATIONS
--------------
- Complete-case analysis (15 patients removed); multiple imputation would be
  more rigorous.
- C-index was not corrected for optimism and no external validation was done.
- Observational data cannot prove causation; unmeasured factors such as stage
  and treatment may confound the associations.
 
