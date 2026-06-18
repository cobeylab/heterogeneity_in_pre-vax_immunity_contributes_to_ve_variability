library(here)
library(tidyverse)
library(cowplot)
library(paletteer)

source(here("src", "full_VE_model.R"))
source(here("src", "data_generation.R"))

# Scenario options (no waning, no continuous risk distributions, VE from cumulative attack rates)
opts = list(
    waning = FALSE,
    heterogeneity = FALSE,
    instantaneous = FALSE
)

# Heatmap scenario parameters
pars <- list(
    start_time = 200,
    end_time = 200,
    dt = 0,
    lambda = c(0.001, 0.005, 0.01),
    theta_0 = 1 - seq(0, 1, 0.01),
    epsilon_v = 1,
    epsilon_u = 1,
    alpha_v = Inf,
    alpha_u = Inf
)

# Generate VE data
dt <- generate_par_sets(pars, include_early = FALSE) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        estd_ve = map_vec(pars, function(x) estimate_math_ve(x, opts))
    ) %>%
    unnest(pars) %>%
    mutate(
        true_vax_protect = (1 - theta_0) * 100
    )

# Plot VE vs true vaccine protection across lambda values
plt <- ggplot(dt) +
    aes(x = true_vax_protect, y = estd_ve, color = factor(lambda)) +
    geom_line(aes(y = true_vax_protect), color = "black", linetype = "44") +
    geom_line(linewidth = 1) +
    scale_color_manual(
        name = "Exogeneous infection hazard",
        breaks = c(0.001, 0.005, 0.01),
        values = c("dodgerblue", "darkorange", "darkorchid")
    ) +
    theme_cowplot(14) +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.05, 0.845),
        legend.box.background = element_rect(fill = "#ffffffaa", color = "#ffffff00")
    ) +
    background_grid() +
    labs(x = "True vaccine protection (%)", y = "VE estimate (%)")

fig_path <- here("plots", "supplemental_figs")
dir.create(fig_path)

ggsave(
    here(fig_path, "lewnard_2018_fig_1c_recreation.png"),
    plt,
    bg = "white",
    height = 5,
    width = 5,
    units = "in"
)

ggsave(
    here(fig_path, "lewnard_2018_fig_1c_recreation.pdf"),
    plt,
    bg = "white",
    height = 5,
    width = 5,
    units = "in"
)