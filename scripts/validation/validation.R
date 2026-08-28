library(here)
library(tidyverse)
library(survival)

source(here("src", "full_VE_model.R"))
source(here("src", "data_generation.R"))

vax_term <- function(t, p) {
    p$alpha_v / (p$alpha_v + (p$epsilon_v * p$theta_0 * p$lambda * t))
}

unvax_term <- function(t, p) {
    p$alpha_u / (p$alpha_u + (p$epsilon_u * p$lambda * t))
}

# main text equation 1a
ve_instantaneous <- function(p) {
    numer <- vax_term(p$time, p)^(p$alpha_v + 1)
    denom <- unvax_term(p$time, p)^(p$alpha_u + 1)

    ret <- p$theta_0 * (p$epsilon_v / p$epsilon_u) * (numer / denom)
    return((1 - ret) * 100)
}

# main text equation 1b
ve_cumulative <- function(p) {
    numer <- 1 - (vax_term(p$time, p)^(p$alpha_v))
    denom <- 1 - (unvax_term(p$time, p)^(p$alpha_u))

    return((1 - (numer / denom)) * 100)
}

# main text equation 1c
ve_cohort <- function(p) {
    suscep_v <- function(t, p) {
        return(vax_term(t, p)^(p$alpha_v))
    }

    suscep_u <- function(t, p) {
        return(unvax_term(t, p)^(p$alpha_u))
    }

    car_v <- 1 - suscep_v(p$time, p)
    car_u <- 1 - suscep_u(p$time, p)

    e_pyar_v <- integrate(
        f = suscep_v,
        lower = 0,
        upper = p$time,
        p = p
    )$value

    e_pyar_u <- integrate(
        f = suscep_u,
        lower = 0,
        upper = p$time,
        p = p
    )$value

    numer <- car_v / e_pyar_v
    denom <- car_u / e_pyar_u

    return((1 - (numer / denom)) * 100)
}

# supplemental equation S.12 with waning
theta <- function(t, p) {
    exponent <- -t * (log(2) / p$eta)

    return(1 - ((1 - p$theta_0) * exp(exponent)))
}

# supplemental equation S.13 with waning
Theta <- function(t, p) {
    numer <- (2^(-t / p$eta)) * (2^(t / p$eta) - 1) * p$eta * (p$theta_0 - 1)

    return(t + (numer / log(2)))
}

# supplemental equation S.48
ve_waning_instantaneous <- function(p) {
    Theta_t <- Theta(p$time, p)
    theta_t <- theta(p$time, p)

    exponent <- -p$lambda * ((p$epsilon_v * Theta_t) - (p$epsilon_u * p$time))
    ret <- theta_t * (p$epsilon_v / p$epsilon_u) * exp(exponent)

    return((1 - ret) * 100)
}

# supplemental equation S.49
ve_waning_cumulative <- function(p) {
    Theta_t <- Theta(p$time, p)
    theta_t <- theta(p$time, p)

    numer <- 1 - exp(-p$epsilon_v * p$lambda * Theta_t)
    denom <- 1 - exp(-p$epsilon_u * p$lambda * p$time)

    return((1 - (numer / denom)) * 100)
}

ve_waning_cohort <- function(p) {
    # supplemental equation S.45
    suscep_v <- function(t, p) {
        return(exp(-p$lambda * p$epsilon_v * Theta(t, p)))
    }

    # supplemental equation S.44
    suscep_u <- function(t, p) {
        return(exp(-p$lambda * p$epsilon_u * t))
    }

    car_v <- 1 - suscep_v(p$time, p)
    car_u <- 1 - suscep_u(p$time, p)

    e_pyar_v <- integrate(
        f = suscep_v,
        lower = 0,
        upper = p$time,
        p = p
    )$value

    e_pyar_u <- integrate(
        f = suscep_u,
        lower = 0,
        upper = p$time,
        p = p
    )$value

    numer <- car_v / e_pyar_v
    denom <- car_u / e_pyar_u

    return((1 - (numer / denom)) * 100)
}

count_diffs_above_threshold <- function(.df, threshold) {
    .df %>%
        filter(ve_diff > threshold) %>%
        nrow()
}

count_na_diffs <- function(.df) {
    .df %>%
        filter(is.na(ve_diff)) %>%
        nrow()
}

print_validation_summary <- function(insnt_dt, cumul_dt, cohrt_dt, ve_diff_threshold) {
    n_large_insnt_ve_diffs <- count_diffs_above_threshold(insnt_dt, ve_diff_threshold)
    n_large_cumul_ve_diffs <- count_diffs_above_threshold(cumul_dt, ve_diff_threshold)
    n_large_cohrt_ve_diffs <- count_diffs_above_threshold(cohrt_dt, ve_diff_threshold)

    n_na_insnt_ve_diffs <- count_na_diffs(insnt_dt)
    n_na_cumul_ve_diffs <- count_na_diffs(cumul_dt)
    n_na_cohrt_ve_diffs <- count_na_diffs(cohrt_dt)

    print(paste0("Number of times the difference between VE from code implementation and analytical equations is greater than ", ve_diff_threshold))
    print(paste0("(number of VE differences that are NA)"))
    print(paste0("VE^instantaneous: ", n_large_insnt_ve_diffs, " (", n_na_insnt_ve_diffs, " NAs)"))
    print(paste0("VE^cumulative: ", n_large_cumul_ve_diffs, " (", n_na_cumul_ve_diffs, " NAs)"))
    print(paste0("VE^cohort: ", n_large_cohrt_ve_diffs, " (", n_na_cohrt_ve_diffs, " NAs)"))
}

### NO WANING
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

# Parameters to estimate final VE 
# start time is set near 0, not at zero to avoid irrelavent NAs at time = 0
pars = list(
    start_time = 1e-3,
    end_time = 200,
    dt = 10,
    lambda = c(0.001, 0.005, 0.01),
    theta_0 = 1 - seq(0.1, 0.9, 0.1),
    epsilon_v = seq(0.1, 1.0, 0.1),
    epsilon_u = seq(0.1, 1.0, 0.1),
    alpha_v = c(0.2, 2, 20),
    alpha_u = c(0.2, 2, 20)
)

cumul_ve_comp <- generate_par_sets(pars, include_early = FALSE) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        cumul_ve1 = map_vec(pars, function(x) estimate_math_ve(x, cumul_ve_opts)),
        cumul_ve2 = map_vec(pars, function(x) ve_cumulative(x)),
        ve_diff = abs(cumul_ve1 - cumul_ve2)
    ) %>%
    arrange(desc(ve_diff))

insnt_ve_comp <- generate_par_sets(pars, include_early = FALSE) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        insnt_ve1 = map_vec(pars, function(x) estimate_math_ve(x, insnt_ve_opts)),
        insnt_ve2 = map_vec(pars, function(x) ve_instantaneous(x)),
        ve_diff = abs(insnt_ve1 - insnt_ve2)
    ) %>%
    arrange(desc(ve_diff))

cohrt_ve_comp <- generate_par_sets(pars, include_early = FALSE) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        cohrt_ve1 = map_vec(pars, function(x) estimate_ve_cohort(x, cumul_ve_opts)),
        cohrt_ve2 = map_vec(pars, function(x) ve_cohort(x)),
        ve_diff = abs(cohrt_ve1 - cohrt_ve2)
    ) %>%
    arrange(desc(ve_diff))

print("Summary of validation for non-waning models (main text equations 1a,b,c)")
print("------------------------------------------------------------------------")
print_validation_summary(insnt_ve_comp, cumul_ve_comp, cohrt_ve_comp, 1e-10)

# clean up memory before next validation
rm(insnt_ve_comp, cumul_ve_comp, cohrt_ve_comp)

### WITH WANING
# Scenario options (no waning, continuous pre-vaccination risk, VE from cumulative attack rates)
opts = list(
    waning = TRUE,
    heterogeneity = FALSE,
    instantaneous = FALSE
)
cumul_ve_opts <- opts

# Create copy of options for estimating VE from instantaneous incidence rates
insnt_ve_opts <- opts
insnt_ve_opts$instantaneous <- TRUE

# Parameters to estimate final VE 
# start time is set near 0, not at zero to avoid irrelavent NAs at time = 0
pars = list(
    start_time = 1e-3,
    end_time = 200,
    dt = 10,
    lambda = c(0.001, 0.005, 0.01),
    theta_0 = 1 - c(0.3, 0.6, 0.9),
    eta = c(30, 360, 1440),
    epsilon_v = seq(0.1, 1.0, 0.1),
    epsilon_u = seq(0.1, 1.0, 0.1),
    alpha_v = c(0.2, 2, 20),
    alpha_u = c(0.2, 2, 20)
)

cumul_ve_comp <- generate_par_sets(pars, include_early = FALSE) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        cumul_ve1 = map_vec(pars, function(x) estimate_math_ve(x, cumul_ve_opts)),
        cumul_ve2 = map_vec(pars, function(x) ve_waning_cumulative(x)),
        ve_diff = abs(cumul_ve1 - cumul_ve2)
    ) %>%
    arrange(desc(ve_diff))

insnt_ve_comp <- generate_par_sets(pars, include_early = FALSE) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        insnt_ve1 = map_vec(pars, function(x) estimate_math_ve(x, insnt_ve_opts)),
        insnt_ve2 = map_vec(pars, function(x) ve_waning_instantaneous(x)),
        ve_diff = abs(insnt_ve1 - insnt_ve2)
    ) %>%
    arrange(desc(ve_diff))

cohrt_ve_comp <- generate_par_sets(pars, include_early = FALSE) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        cohrt_ve1 = map_vec(pars, function(x) estimate_ve_cohort(x, cumul_ve_opts)),
        cohrt_ve2 = map_vec(pars, function(x) ve_waning_cohort(x)),
        ve_diff = abs(cohrt_ve1 - cohrt_ve2)
    ) %>%
    arrange(desc(ve_diff))

print("Summary of validation for waning models (supp equations S.44,45,48,49)")
print("----------------------------------------------------------------------")
print_validation_summary(insnt_ve_comp, cumul_ve_comp, cohrt_ve_comp, 1e-10)
