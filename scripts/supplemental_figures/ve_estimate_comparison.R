library(here)
library(tidyverse)
library(survival)

source(here("src", "full_VE_model.R"))
source(here("src", "data_generation.R"))

### Final ve sweeping across mean pre-vaccination risk and true vaccine protection (NO WANING)

set.seed(0)

# Scenario parameters
pars = list(
    start_time = 0,
    end_time = 200,
    dt = 14,
    lambda = 0.005,
    lambda_negative = 0.015,
    theta_0 = 1 - seq(0.1, 0.9, 0.1),
    eta = -1,
    epsilon_v = seq(0.1, 1.0, 0.1),
    epsilon_u = seq(0.1, 1.0, 0.1),
    alpha_v = 20,
    alpha_u = 20,
    mu_v = 1,
    mu_u = 1,
    pi_pos = 1,
    pi_neg = 1,
    c = 0.5,
    N = 5e3
)

pars$time <- pars$end_time

# Scenario options (no waning, continuous pre-vaccination risk, VE from cumulative attack rates)
opts = list(
    waning = FALSE,
    heterogeneity = TRUE,
    instantaneous = FALSE
)
cumul_ve_opts <- opts

# Create copy of options for estimating VE from instantaneous incidence rates
insnt_ve_opts <- opts
insnt_ve_opts$instantaneous <- TRUE

# Generate parameter sets
pars_dt <- crossing(!!!pars) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup()

print(paste0("generating ", nrow(pars_dt)," no waning scenario linelists"))

# Generate VE data (analytical and logistic regression models)
final_ve_dt <- pars_dt %>%
    rowwise_mutate(counts = list(generate_counts(pars, opts))) %>%
    mutate(linelist = map(counts, function(x) generate_linelist(x))) %>%
    mutate(
        insnt_ve = map_vec(pars, function(x) estimate_math_ve(x, insnt_ve_opts)),
        cumul_ve = map_vec(pars, function(x) estimate_math_ve(x, cumul_ve_opts)),
        cohrt_ve = map_vec(pars, function(x) estimate_ve_cohort(x, cumul_ve_opts)),
        ulreg_ve = map_vec(linelist, function(x) estimate_ve_uncond_logreg(x)),
        clreg_ve = map_vec(linelist, function(x) estimate_ve_cond_logreg(x))
    )

# Calculate absolute average difference between VE estimates
no_waning_ve_compare <- final_ve_dt %>%
    select(ends_with("_ve")) %>%
    pivot_longer(!cumul_ve) %>%
    mutate(ve_diff = value - cumul_ve) %>%
    group_by(name) %>%
    arrange(ve_diff) %>%
    summarize(
        mean_diff = mean(ve_diff, na.rm = TRUE),
        lower_diff = quantile(ve_diff, probs = c(0.025), na.rm = TRUE),
        upper_diff = quantile(ve_diff, probs = c(0.975), na.rm = TRUE)
    )

print("Mean and 95% interval of difference between VE estimates and VE^cumulative")
print("Scenarios: varying true vax protection and pre-vax risk means")
print(no_waning_ve_compare)

rm(pars_dt, final_ve_dt)

### Final ve sweeping across waning rate, mean pre-vaccination risk, and
### true vaccine protection (WITH WANING)

set.seed(0)

# Scenario parameters
pars = list(
    start_time = 0,
    end_time = 200,
    dt = 14,
    lambda = 0.005,
    lambda_negative = 0.015,
    theta_0 = 1 - c(0.3, 0.6, 0.9),
    eta = c(30, 180, 360, 1440),
    epsilon_v = seq(0.1, 1.0, 0.1),
    epsilon_u = seq(0.1, 1.0, 0.1),
    alpha_v = Inf,
    alpha_u = Inf,
    mu_v = 1,
    mu_u = 1,
    pi_pos = 1,
    pi_neg = 1,
    c = 0.5,
    N = 5e3
)

pars$time <- pars$end_time

# Scenario options (waning, point-density pre-vaccination risk, VE from cumulative attack rates)
opts = list(
    waning = TRUE,
    heterogeneity = FALSE,
    instantaneous = FALSE
)
cumul_ve_opts <- opts

# Create copy of options for estimating VE from instantaneous incidence rates
insnt_ve_opts <- opts
insnt_ve_opts$instantaneous <- TRUE

# Generate parameter sets
pars_dt <- crossing(!!!pars) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup()

print(paste0("generating ", nrow(pars_dt)," waning scenario linelists"))

# Generate VE data (analytical and logistic regression models)
final_ve_dt <- pars_dt %>%
    rowwise_mutate(counts = list(generate_counts(pars, opts))) %>%
    mutate(linelist = map(counts, function(x) generate_linelist(x))) %>%
    mutate(
        insnt_ve = map_vec(pars, function(x) estimate_math_ve(x, insnt_ve_opts)),
        cumul_ve = map_vec(pars, function(x) estimate_math_ve(x, cumul_ve_opts)),
        cohrt_ve = map_vec(pars, function(x) estimate_ve_cohort(x, cumul_ve_opts)),
        ulreg_ve = map_vec(linelist, function(x) estimate_ve_uncond_logreg(x)),
        clreg_ve = map_vec(linelist, function(x) estimate_ve_cond_logreg(x))
    ) %>%
    unnest(pars)

# Calculate absolute average difference between VE estimates
waning_ve_compare <- final_ve_dt %>%
    select(ends_with("_ve")) %>%
    pivot_longer(!cumul_ve) %>%
    mutate(ve_diff = value - cumul_ve) %>%
    group_by(name) %>%
    arrange(ve_diff) %>%
    summarize(
        mean_diff = mean(ve_diff, na.rm = TRUE),
        lower_diff = quantile(ve_diff, probs = c(0.025), na.rm = TRUE),
        upper_diff = quantile(ve_diff, probs = c(0.975), na.rm = TRUE)
    )

print("Mean and 95% interval of difference between VE estimates and VE^cumulative")
print("Scenarios: varying true vax protection, pre-vax risk means and vax waning")
print(waning_ve_compare)