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
        v_car = (v_cumul_inc / vaxd) * 100,
        u_car = (u_cumul_inc / (pars$pop_size - vaxd)) * 100,
        ve_cumul = (1 - (v_car / u_car)) * 100
    ) %>%
    select(time = t, ends_with("car"), ve_cumul) %>%
    pivot_longer(!time) %>%
    mutate(data = "sim")

# read in simulation data
print("ESTIMATING VE FROM ANALYTICAL MODELS")

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
    theta_0 = 1 - pars$true_vax_protection,
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
        v_car = map_vec(pars, function(x) cumulative_attack_rate(x$time, x, opts, TRUE) * 100),
        u_car = map_vec(pars, function(x) cumulative_attack_rate(x$time, x, opts, FALSE) * 100),
        ve_cumul = map_vec(pars, function(x) estimate_math_ve(x, opts))
    ) %>%
    unnest(pars) %>%
    select(time, ends_with("car"), ve_cumul) %>%
    pivot_longer(!time) %>%
    mutate(data = "math")

# read in simulation data
print("SAVING PLOT")

name_labs <- c(
    "u_car" = "Unvaccinated cumulative\nattack rate",
    "v_car" = "Vaccinated cumulative\nattack rate",
    "ve_cumul" = "Cumulative-attack-rate\nVE estimate"
)

plt <- bind_rows(sim_dt, math_dt) %>%
    ggplot() +
        aes(x = time, y = value, color = data, linetype = data) +
        geom_line(linewidth = 2) +
        facet_wrap(
            vars(name),
            labeller = labeller(
                name = name_labs
            )
        ) +
        scale_linetype_manual(
            name = "Data source",
            breaks = c("math", "sim"),
            values = c("solid", "22"),
            labels = c("Analytical model", "Simulation model")
        ) +
        scale_color_manual(
            name = "Data source",
            breaks = c("math", "sim"),
            values = c("dodgerblue", "darkorange"),
            labels = c("Analytical model", "Simulation model")
        ) +
        coord_cartesian(ylim = c(NA, 100)) +
        theme_cowplot(20) +
        theme(legend.position = "top") +
        background_grid() +
        labs(x = "Time (days)", y = "Value (%)")

fig_dir <- here("plots", "supplemental_figs")
dir.create(fig_dir)

ggsave(
    here(fig_dir, "sim_validation.png"),
    plt,
    width = 15,
    height = 5,
    units = "in",
    bg = "white"
)

ggsave(
    here(fig_dir, "sim_validation.pdf"),
    plt,
    width = 15,
    height = 5,
    units = "in",
    bg = "white"
)