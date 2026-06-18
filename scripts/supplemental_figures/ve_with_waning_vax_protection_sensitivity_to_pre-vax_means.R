library(here)
library(tidyverse)
library(cowplot)
library(paletteer)

source(here("src", "full_VE_model.R"))
source(here("src", "data_generation.R"))

# Scenario options (waning, no continuous risk distributions, VE from cumulative attack rates)
opts = list(
    waning = TRUE,
    heterogeneity = FALSE,
    instantaneous = FALSE
)

# Scenario parameters
pars <- list(
    start_time = 0,
    end_time = 200,
    dt = 1,
    lambda = 0.0015,
    theta_0 = 1 - 0.3,
    epsilon_v = 1,
    epsilon_u = c(1, 0.9, 0.8),
    alpha_v = Inf,
    alpha_u = Inf,
    eta = c(Inf, 90)
)

# Generate VE data
dt <- generate_par_sets(pars) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        tot_vax_effect = map_vec(pars, function(x) total_vaccine_effect(x$time, x, opts)),
        estd_ve = map_vec(pars, function(x) estimate_math_ve(x, opts))
    ) %>%
    unnest(pars) %>%
    mutate(avg_vax_protect = (1 - (tot_vax_effect / time)) * 100) %>%
    select(time, eta, epsilon_u, estd_ve, avg_vax_protect) %>%
    pivot_longer(c(estd_ve, avg_vax_protect))

# Consistent plot aesthetics
lt <- c(
    "estd_ve" = "solid",
    "avg_vax_protect" = "22"
)

cl <- c(
    "estd_ve.1" = "#d95f02",
    "estd_ve.0.9" = "#1b9e77",
    "estd_ve.0.8" = "#7570b3",
    "avg_vax_protect.1" = "#000000",
    "avg_vax_protect.0.9" = "#000000",
    "avg_vax_protect.0.8" = "#000000"
)

labs <- c(
    "estd_ve" = "Vaccine effectiveness",
    "avg_vax_protect" = "True average\nvaccine protection"
)

eta_labs <- c(
    "Inf" = "Constant vaccine protection",
    "90" = "Vaccine protection wanes with 90-day half life"
)

# Plot VE vs average true vaccine protection over time with and without waning
ve_plt <- ggplot(dt) +
    aes(
        x = time,
        y = value,
        linetype = factor(name),
        color = interaction(name, epsilon_u)
    ) +
    geom_hline(yintercept = 0, linetype = "44", color = "gray25") +
    geom_line(linewidth = 1) +
    facet_wrap(
        vars(factor(eta, levels = c(Inf, 90), labels = eta_labs)),
        ncol = 1
    ) +
    coord_cartesian(ylim = c(-10, 31)) +
    scale_linetype_manual(
        name = NULL,
        values = lt,
        labels = labs
    ) +
    scale_color_manual(
        name = "Unvaccinated mean\nbaseline susceptibility",
        values = cl,
        labels = c("1.0", "0.9", "0.8"),
        breaks = c("estd_ve.1", "estd_ve.0.9", "estd_ve.0.8")
    ) +
    theme_cowplot() +
    background_grid() +
    theme(
        legend.position = "top",
        legend.position.inside = c(0.5, 0.5),
        legend.background = element_rect(fill = "#ffffff99"),
        legend.key.spacing.y = unit(0.001, "npc"),
        legend.box = "vertical"
    ) +
    labs(x = "Time (days)", y = "Vaccine effectiveness (%)")

fig_dir <- here("plots", "supplemental_figs")
dir.create(fig_dir)

ggsave(
    here(fig_dir, "ve_with_waning_vs_unvax_pre-vax_mean.png"),
    ve_plt,
    width = 5,
    height = 5,
    units = "in",
    bg = "white"
)

ggsave(
    here(fig_dir, "ve_with_waning_vs_unvax_pre-vax_mean.pdf"),
    ve_plt,
    width = 5,
    height = 5,
    units = "in",
    bg = "white"
)