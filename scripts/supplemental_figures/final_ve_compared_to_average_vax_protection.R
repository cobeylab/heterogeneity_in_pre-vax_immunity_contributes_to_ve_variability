library(here)
library(tidyverse)
library(cowplot)
library(paletteer)

source(here("src", "full_VE_model.R"))
source(here("src", "data_generation.R"))

# # Helper function to clean up cell label values that are close to zero.
# clean_label <- function(x, sigdig) {
#   return(ifelse(abs(x) < 1e-7, 0, signif(x, digits = sigdig)))
# }

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
    lambda = 0.01,
    theta_0 = 1 - seq(0.1, 1, 0.1),
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
        bias = estd_ve - avg_vax_protect
    )

# Find true vaccine protection at which bias is highest when there is no waning
intermediate_VE <- round(dt$avg_vax_protect[abs(dt$bias) == max(abs(dt$bias))])

# Panel A: heatmap of final VE bias across starting vaccine protection and waning rate
bias_plt <- ggplot(dt) +
    aes(
        x = eta_f,
        y = factor(starting_vax_protect),
        fill = bias
    ) +
    geom_tile() +
    geom_vline(xintercept = 7.5, color = "white", linewidth = 2) +
    scale_fill_paletteer_c(
        "grDevices::Reds",
        name = "Difference between final VE and true\naverage vaccine protection (% pts)",
        n.breaks = 3,
        limits = c(NA, 0)
    ) +
    theme_cowplot(14) +
    theme(legend.position = "top") +
    labs(
        x = "True vaccine protection half life (days)",
        y = "Starting true vaccine protection (%)"
    )

# Panel B: heatmap of average vaccine protection across starting protection and waning rate
avg_eff_plt <- ggplot(dt) +
    aes(
        x = eta_f,
        y = factor(starting_vax_protect),
        fill = avg_vax_protect
    ) +
    geom_tile() +
    geom_vline(xintercept = 7.5, color = "white", linewidth = 2) +
    scale_fill_gradient2(
        name = "True average vaccine\nprotection after 200 days (%)",
        low = "white",
        mid = "dodgerblue",
        high = "white",
        midpoint = intermediate_VE,
        limits = c(0, 100),
        breaks = c(0, intermediate_VE, 100)
    ) +
    theme_cowplot() +
    theme(legend.position = "top") +
    labs(
        x = "True vaccine protection half life (days)",
        y = "Starting true vaccine protection (%)"
    )

# Combine panels and save
plt <- plot_grid(
  bias_plt, avg_eff_plt,
  nrow = 1,
  labels = "AUTO"
)

fig_dir <- here("plots", "supplemental_figs")
dir.create(fig_dir)

ggsave(
    here(fig_dir, "final_ve_compared_to_avg_vax_protect.png"),
    plt,
    width = 10,
    height = 5,
    units = "in",
    bg = "white"
)

ggsave(
    here(fig_dir, "final_ve_compared_to_avg_vax_protect.pdf"),
    plt,
    width = 10,
    height = 5,
    units = "in",
    bg = "white"
)
