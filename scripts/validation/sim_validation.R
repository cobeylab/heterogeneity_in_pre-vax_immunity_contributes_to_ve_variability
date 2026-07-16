library(here)
library(tidyverse)
library(optparse)
library(RcppTOML)
library(cowplot)
library(paletteer)

source(here("src", "full_VE_model.R"))
source(here("src", "data_generation.R"))

# parse command line arguments
if (interactive()) {
    parser <- OptionParser()
    parser <- add_option(parser, "--config", help = "path to config TOML file (rel path from this script)")
    parsed_args <- parse_args(parser, args = c("--config=config.toml"))
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
script_dir <- here("scripts", "validation")
config_path <- here(script_dir, parsed_args[["config"]])

config <- parseTOML(config_path)
pars <- config$parameters
exp_name <- pars$exp_name

# read in simulation data
print("READING DATA")

dat_dir <- here("data", exp_name)
exp_filepath <- here(dat_dir, "experiments.csv")
res_filepath <- here(dat_dir, "results.csv")
suscep_filepath <- here(dat_dir, "suscep", "suscep_1-1.csv")

experiments <- read_csv(exp_filepath, show_col_types = FALSE)
results <- read_csv(res_filepath, show_col_types = FALSE)
suscep <- read_csv(suscep_filepath, show_col_types = FALSE)

sim_dt <- results %>%
    mutate(
        step = t / pars$dt
    ) %>%
    group_by(exp, t) %>%
    mutate(
        v_inc = mean(vax_inf),
        u_inc = mean(unvax_inf)
    ) %>%
    ungroup() %>%
    mutate(
        v_cumul_inc = cumsum(v_inc),
        u_cumul_inc = cumsum(u_inc),
        v_car = v_cumul_inc / vaxd,
        u_car = u_cumul_inc / (pars$pop_size - vaxd),
        ve_cumul = (1 - (v_car / u_car)) * 100
    ) %>%
    select(time = t, ends_with("car"), ve_cumul) %>%
    pivot_longer(!time) %>%
    mutate(data = "sim")

# Scenario options (no waning, continuous risk distributions, VE from cumulative attack rates)
opts = list(
    waning = FALSE,
    heterogeneity = TRUE,
    instantaneous = FALSE
)

# Parameters to draw cumulative attack rate curves
pars <- list(
    start_time = 0,
    end_time = pars$tmax,
    dt = pars$dt,
    lambda = pars$beta,
    theta_0 = 1 - pars$vax_efficacy,
    epsilon_v = pars$vax_mean_pre_vax_suscep,
    epsilon_u = pars$unvax_mean_pre_vax_suscep,
    alpha_v = pars$vax_pre_vax_suscep_shape,
    alpha_u = pars$unvax_pre_vax_suscep_shape
)

math_dt <- generate_par_sets(pars, include_early = FALSE) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        v_car = map_vec(pars, function(x) cumulative_attack_rate(x$time, x, opts, TRUE)),
        u_car = map_vec(pars, function(x) cumulative_attack_rate(x$time, x, opts, FALSE)),
        ve_cumul = map_vec(pars, function(x) estimate_math_ve(x, opts))
    ) %>%
    unnest(pars) %>%
    select(time, ends_with("car"), ve_cumul) %>%
    pivot_longer(!time) %>%
    mutate(data = "math")

bind_rows(sim_dt, math_dt) %>%
    ggplot() +
        aes(x = time, y = value, color = data) +
        geom_line() +
        facet_wrap(vars(name), scales = "free") +
        theme_cowplot() +
        background_grid()
