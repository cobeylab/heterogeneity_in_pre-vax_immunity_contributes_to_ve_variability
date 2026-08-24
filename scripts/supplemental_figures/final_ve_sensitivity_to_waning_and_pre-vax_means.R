library(here)
library(tidyverse)
library(cowplot)
library(paletteer)

source(here("src", "full_VE_model.R"))
source(here("src", "data_generation.R"))

# Helper function to label scenarios when VE crosses zero
crosses_zero <- function(start_eff, VE) {
    starts_pos <- start_eff > 0
    ends_neg <- VE < 0
    if (starts_pos && !ends_neg) {
        return("pos2pos")
    } else if (!starts_pos && ends_neg) {
       return("neg2neg")
    } else if (starts_pos && ends_neg) {
        return("pos2neg")
    } else {
        return("neg2pos")
    }
}

# Helper function to label scenarios when VE crosses average true vaccine protection
crosses_truth <- function(start, end, avg) {
    starts_pos <- start > 0
    ends_neg <- end < 0
    start_over <- start > avg
    end_under <- end < avg

    if (starts_pos && !ends_neg) {
        # positive and...
        if (start_over && !end_under) {
            # always overestimates
            return("ovr2ovr")
        } else if (!start_over && end_under) {
            # always underestimates
            return("und2und")
        } else if (start_over && end_under) {
            # over- then underestimates
            return("ovr2und")
        } else {
            return("pos.other")
        }
    } else if (starts_pos && ends_neg) {
        # crosses zero
        return("pos2neg")
    } else if (!starts_pos && ends_neg) {
        # negative
        return("neg2neg")
    }
}

# Scenario options (waning, no continuous risk distributions, VE from cumulative attack rates)
opts = list(
    waning = TRUE,
    heterogeneity = FALSE,
    instantaneous = FALSE
)

# Scenario parameters
pars <- list(
    start_time = 200,
    end_time = 200,
    dt = 0,
    lambda = 0.005,
    theta_0 = 1 - c(0.3, 0.6, 0.9),
    epsilon_v = seq(0.01, 1.0, 0.01),
    epsilon_u = seq(0.01, 1.0, 0.01),
    alpha_v = Inf,
    alpha_u = Inf,
    eta = c(30, 90, 360)
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
        starting_vax_protect = (1 - theta_0) * 100,
        epsilon_ratio = epsilon_v / epsilon_u,
        starting_vax_protect_adj = (1 - (epsilon_ratio * theta_0)) * 100,
        avg_vax_protect = (1 - (tot_vax_effect / time)) * 100
    ) %>%
    rowwise_mutate(
        neg_ve = crosses_zero(starting_vax_protect_adj, estd_ve),
        overunder = crosses_truth(starting_vax_protect_adj, estd_ve, avg_vax_protect)
    )

# Consistent plot aesthetics
eta_labs <- c(
    "30" = "Vaccine waning half life of 30 days",
    "90" = "Vaccine waning half life of 90 days",
    "360" = "Vaccine waning half life of 1 year"
)

cl <- c(
    "pos2pos" = "#008837",
    "ovr2ovr" = "#d7191c",
    "ovr2und" = "#fdae61",
    "und2und" = "#ffffbf",
    "pos2neg" = "#abd9e9",
    "neg2neg" = "#2c7bb6",
    "neg2pos" = "#c2a5cf"
)

labs <- c(
    "ovr2ovr" = "Overestimates average protection",
    "ovr2und" = "Over- then underestimates average protection",
    "und2und" = "Underestimates average protection",
    "neg2neg" = "Starts and ends negative",
    "pos2neg" = "Starts postive and ends negative"
)

# Heatmap showing final VE outcomes across waning and mean pre-vax scenarios
plt <- ggplot(dt) +
    aes(
        x = epsilon_u,
        y = epsilon_v,
        fill = factor(
            overunder,
            levels = c("ovr2ovr", "ovr2und", "und2und", "pos2neg", "neg2neg")
        )
    ) +
    geom_tile(color = NA) +
    facet_grid(
        rows = vars(factor(eta, levels = c(30, 90, 360), labels = eta_labs)),
        cols = vars(starting_vax_protect),
        labeller = labeller(
            starting_vax_protect = function(x) paste0("Starting true vaccine protection = ", x, "%")
        )
    ) +
    scale_fill_manual(
        name = expression("Final" ~ widehat(VE) ~ "after 200 days..."),
        labels = labs,
        values = cl,
    ) +
    coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05), expand = FALSE) +
    theme_cowplot() +
    background_grid() +
    theme(legend.position = "top") +
    guides(fill = guide_legend(nrow = 3)) +
    labs(
        x = "Average unvaccinated pre-vaccination susceptibility",
        y = "Average vaccinated pre-vaccination susceptibility"
    )

fig_dir <- here("plots", "supplemental_figs")
dir.create(fig_dir)

ggsave(
    here(fig_dir, "final_ve_outcomes_vs_waning_and_pre-vax_means.png"),
    plt,
    width = 10,
    height = 10,
    units = "in",
    bg = "white"
)

ggsave(
    here(fig_dir, "final_ve_outcomes_vs_waning_and_pre-vax_means.pdf"),
    plt,
    width = 10,
    height = 10,
    units = "in",
    bg = "white"
)
