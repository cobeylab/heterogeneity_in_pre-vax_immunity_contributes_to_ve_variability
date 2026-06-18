library(here)
library(tidyverse)
library(cowplot)
library(paletteer)

source(here("src", "full_VE_model.R"))
source(here("src", "data_generation.R"))

# Helper function to clean up cell label values that are close to zero.
clean_label <- function(x, sigdig) {
  return(ifelse(abs(x) < 1e-7, 0, signif(x, digits = sigdig)))
}

# Scenario options (no waning, no continuous risk distributions, VE from cumulative attack rates)
opts = list(
    waning = TRUE,
    heterogeneity = FALSE,
    instantaneous = FALSE
)

# Heatmap scenario parameters
pars <- list(
    start_time = 200,
    end_time = 200,
    dt = 0,
    lambda = seq(0.001, 0.01, 0.001),
    theta_0 = 1 - c(0.3, 0.6, 0.9),
    epsilon_v = 1,
    epsilon_u = 1,
    alpha_v = Inf,
    alpha_u = Inf,
    eta = c(30, 60, 90, 180, 360, 720, 1440, Inf)
)

# Generate VE data
dt <- generate_par_sets(pars, include_early = FALSE) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        tot_vax_effect = map_vec(pars, function(x) total_vaccine_effect(x$time, x, opts)),
        estd_ve = map_vec(pars, function(x) estimate_math_ve(x, opts))
    ) %>%
    unnest(pars) %>%
    mutate(
        eta_f = factor(eta, levels = c(30, 60, 90, 180, 360, 720, 1440, Inf),
                            labels = c(30, 60, 90, 180, 360, 720, 1440, "No\nwaning")),
        starting_vax_protect = (1 - theta_0) * 100,
        avg_vax_protect = (1 - (tot_vax_effect / time)) * 100,
        bias = estd_ve - avg_vax_protect,
        cell_label = clean_label(bias, sigdig = 2)
    )

# Heatmap plot of final VE vs lambda, vax protection waning rates
plt <- ggplot(dt) +
    aes(
        x = eta_f,
        y = lambda,
        fill = bias,
        label = cell_label
    ) +
    geom_tile(color = NA) +
    geom_text() +
    geom_vline(xintercept = 7.5, color = "white", linewidth = 2) +
    facet_wrap(
        vars(starting_vax_protect),
        labeller = labeller(
            starting_vax_protect = function(x) paste0("Starting true vaccine protection = ", x, "%")
        )
    ) +
    scale_fill_paletteer_c(
        "grDevices::Reds",
        name = "Difference between final VE and true average vaccine protection (% pts)",
        n.breaks = 3,
        limits = c(NA, 0)
    ) +
    coord_cartesian(ylim = c(0, NA), expand = FALSE) +
    scale_y_continuous(breaks = c(0, 0.005, 0.01)) +
    theme_cowplot(14) +
    theme(
        legend.position = "top",
        legend.key.width = unit(0.03, "npc"),
        strip.background = element_rect(fill = "gray90"),
    ) +
    labs(
        x = "True vaccine protection half life (days)",
        y = "Exogenous infection hazard"
    )

fig_dir <- here("plots", "supplemental_figs")
dir.create(fig_dir)

ggsave(
    here(fig_dir, "final_ve_vs_vax_waning_rate.png"),
    plt,
    width = 15,
    height = 5,
    units = "in",
    bg = "white"
)

ggsave(
    here(fig_dir, "final_ve_vs_vax_waning_rate.pdf"),
    plt,
    width = 15,
    height = 5,
    units = "in",
    bg = "white"
)