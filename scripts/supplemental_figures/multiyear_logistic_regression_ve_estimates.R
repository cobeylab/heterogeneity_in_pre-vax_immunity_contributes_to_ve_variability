library(here)
library(tidyverse)
library(optparse)
library(RcppTOML)
library(survival)
library(cowplot)
library(paletteer)

# maintain consistent RNG seed
set.seed(0)

# parse command line arguments
if (interactive()) {
    parser <- OptionParser()
    parser <- add_option(parser, "--config", help = "path to config TOML file (rel path from this script)")
    parsed_args <- parse_args(parser, args = c("--config=fig3.toml"))
} else {
    parser <- OptionParser()
    parser <- add_option(parser, "--config", help = "path to config TOML file (rel path from this script)")
    parsed_args <- parse_args(parser)
}

if (is.null(parsed_args[["config"]])) {
  cat("ERROR: --config is required\n")
  quit("no", status = 1, runLast = FALSE)
}

# parse TOML parameter file
script_dir <- here("scripts", "fig3")
config_path <- here(script_dir, parsed_args[["config"]])

config <- parseTOML(config_path)
pars <- config$parameters
exp_name <- pars$exp_name

dat_dir <- here("data", exp_name)
res_filepath <- here(dat_dir, "results.csv")

# calculates the test-negative infection hazard following the same seasonal forcing
# used in the simulation model (i.e., shift = 0.5)
test_neg_infection_hazard <- function(time, beta, amplifier, shift = 0.5) {
    return(beta * (1 + (amplifier * cos(2 * pi * (time + shift)))))
}

# estimate vaccinated test-negative infections at a given time
# (taking into account the time step width)
vax_test_neg_infections <- function(time, pars, shift = 0.5) {
    beta_test_neg <- pars$beta * 0.75
    haz <- test_neg_infection_hazard(time, beta_test_neg, pars$seasonality_amplifier, shift)
    return(pars$pop_size * pars$vax_coverage * haz * pars$dt)
}

# estimate unvaccinated test-negative infections at a given time
# (taking into account the time step width)
unvax_test_neg_infections <- function(time, pars, shift = 0.5) {
    beta_test_neg <- pars$beta * 0.75
    haz <- test_neg_infection_hazard(time, beta_test_neg, pars$seasonality_amplifier, shift)
    return(pars$pop_size * (1 - pars$vax_coverage) * haz * pars$dt)
}

per_10k <- 1e3 / pars$pop_size

# gather test-positive and test-negative infections and adjust per 10k people
dt <- read_csv(res_filepath, show_col_types = FALSE) %>%
    filter(t > 0.005) %>%
    mutate(
        step = t / pars$dt,
        year = ceiling(t),
        time = t - (year - 1)
    ) %>%
    group_by(exp, t) %>%
    mutate(
        v_tp_inf = mean(vax_inf) * per_10k,
        u_tp_inf = mean(unvax_inf) * per_10k,
        v_tn_inf = vax_test_neg_infections(t, pars) * per_10k,
        u_tn_inf = unvax_test_neg_infections(t, pars) * per_10k
    ) %>%
    ungroup() %>%
    select(exp, t, year, ends_with("_tp_inf"), ends_with("_tn_inf")) %>%
    pivot_longer(!c(exp, t, year), values_to = "per_10k_inc") %>%
    mutate(count = as.integer(per_10k_inc / per_10k))

# generate linelists from the counts of cases and controls
vax_tp_linelist <- dt %>%
    filter(name == "v_tp_inf") %>%
    mutate(
        count = as.integer(per_10k_inc / per_10k),
        vax = 1,
        inf = 1
    ) %>%
    select(-c(per_10k_inc, exp, name)) %>%
    uncount(count)

vax_tn_linelist <- dt %>%
    filter(name == "v_tn_inf") %>%
    mutate(
        count = as.integer(per_10k_inc / per_10k),
        vax = 1,
        inf = 0
    ) %>%
    select(-c(per_10k_inc, exp, name)) %>%
    uncount(count)

unvax_tp_linelist <- dt %>%
    filter(name == "u_tp_inf") %>%
    mutate(
        count = as.integer(per_10k_inc / per_10k),
        vax = 0,
        inf = 1
    ) %>%
    select(-c(per_10k_inc, exp, name)) %>%
    uncount(count)

unvax_tn_linelist <- dt %>%
    filter(name == "u_tn_inf") %>%
    mutate(
        count = as.integer(per_10k_inc / per_10k),
        vax = 0,
        inf = 0
    ) %>%
    select(-c(per_10k_inc, exp, name)) %>%
    uncount(count)

# combine the linelists together
# a random sample of the total linelists for each year is used to minimize memory usage
# and generate sample sizes small enough for conditional logistic regression to succeeed
linelist <- bind_rows(vax_tp_linelist, vax_tn_linelist,
                      unvax_tp_linelist, unvax_tn_linelist) %>%
    mutate(block = as.integer((t %/% 0.04) - ((year - 1) * 25))) %>%
    group_by(year) %>%
    slice_sample(prop = 0.0015) %>%
    nest() %>%
    ungroup()

# remove to reduce memory usage
rm(vax_tp_linelist, vax_tn_linelist, unvax_tp_linelist, unvax_tn_linelist)

# helper function to estimate VE using cumulative attack rates
estimate_ve_from_cARs <- function(x) {
    x %>%
        filter(inf == 1) %>%
        group_by(vax) %>%
        summarize(tot_infs = n()) %>%
        ungroup() %>%
        mutate(cAR = tot_infs / (pars$pop_size * if_else(vax == 1, pars$vax_coverage, 1 - pars$vax_coverage))) %>%
        select(-tot_infs) %>%
        pivot_wider(names_from = vax, names_glue = "{vax}_{.value}", values_from = cAR) %>%
        summarize((1 - (`1_cAR` / `0_cAR`)) * 100) %>%
        as.numeric()
}

# helper function to estimate VE using conditional logistic regression
estimate_clreg_ve <- function(x) {
    reg <- clogit(inf ~ vax + strata(block), data = x)
    or <- exp(reg$coefficients["vax"])
    
    return((1 - or) * 100)
}

# helper function to estimate VE using unconditional logistic regression
estimate_ulreg_ve <- function(x) {  
    reg <- glm(inf ~ vax + block, data = x, family = binomial)
    or <- exp(reg$coefficients["vax"])
    
    return((1 - or) * 100)
}

# estimate annual VE estimates and calculate average absolute difference between
# regression-based VE and VE from cumulative attack rates
ret <- linelist %>%
    mutate(
        ve_clreg = map_vec(data, function(x) estimate_clreg_ve(x)),
        ve_ulreg = map_vec(data, function(x) estimate_ulreg_ve(x)),
        ve_cAR = map_vec(data, function(x) estimate_ve_from_cARs(x))
    ) %>%
    select(year, data, starts_with("ve_")) %>%
    pivot_longer(ends_with("reg")) %>%
    mutate(ve_diff = abs(value - ve_cAR)) %>%
    group_by(name) %>%
    summarize(mean_diff = mean(ve_diff))

print("Average absolute difference between regression-based VE and VE from cumulative attack rates")
print(ret)