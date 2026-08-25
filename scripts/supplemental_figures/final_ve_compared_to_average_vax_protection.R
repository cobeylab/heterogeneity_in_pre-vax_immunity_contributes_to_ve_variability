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

# Scenario options (waning vax protection, no continuous risk distributions, VE from cumulative attack rates)
opts = list(
    waning = TRUE,
    heterogeneity = FALSE,
    instantaneous = FALSE
)

# Lewnard plot scenario parameters
pars <- list(
    start_time = 200,
    end_time = 200,
    dt = 0,
    lambda = c(0.001, 0.005, 0.01),
    theta_0 = 1 - seq(0, 1, 0.01),
    epsilon_v = 1,
    epsilon_u = 1,
    alpha_v = Inf,
    alpha_u = Inf,
    eta = Inf
)

# Generate lewnard plot data
lewnard_dt <- generate_par_sets(pars, include_early = FALSE) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        estd_ve = map_vec(pars, function(x) estimate_math_ve(x, opts))
    ) %>%
    unnest(pars) %>%
    mutate(
        true_vax_protect = (1 - theta_0) * 100,
        bias = estd_ve - true_vax_protect
    )

# Find true vaccine protection at which bias is highest when there is no waning
intermediate_VE <- round(lewnard_dt$true_vax_protect[abs(lewnard_dt$bias) == max(abs(lewnard_dt$bias))])

lewnard_plt <- ggplot(lewnard_dt) +
    aes(x = true_vax_protect, y = estd_ve, color = factor(lambda)) +
    geom_line(aes(y = true_vax_protect), color = "black", linetype = "44") +
    geom_segment(
        x = 58,
        xend = 58,
        y = 0,
        yend = 58,
        color = "maroon",
        linetype = "33",
        linewidth = 1
    ) +
    annotate(
        "label",
        x = 82.5,
        y = 10,
        label = "Largest distance between VE\nand true protection occurs at\n58% true vaccine protection",
        color = "maroon",
        size = 3
    ) +
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
    labs(x = "True vaccine protection (%)", y = expression("Final " ~ widehat(VE) ~ " (%)"))

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
heatmap_dt <- dt <- generate_par_sets(pars, include_early = FALSE) %>%
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
        true_vax_protect = (1 - theta_0) * 100
    )

# heatmap of final VE bias across starting vaccine protection and waning rate
bias_plt <- ggplot(heatmap_dt) +
    aes(
        x = eta_f,
        y = factor(starting_vax_protect),
        fill = bias
    ) +
    geom_tile() +
    geom_vline(xintercept = 7.5, color = "white", linewidth = 2) +
    scale_fill_paletteer_c(
        "grDevices::Reds",
        name = expression(atop(
            "Difference between final" ~ widehat(VE) ~ "and true",
            "average vaccine protection (% pts)"
        )),
        n.breaks = 3,
        limits = c(NA, 0)
    ) +
    theme_cowplot(14) +
    theme(legend.position = "top") +
    labs(
        x = "True vaccine protection half life (days)",
        y = "Starting true vaccine protection (%)"
    )

# heatmap of average vaccine protection across starting protection and waning rate
avg_eff_plt <- ggplot(heatmap_dt) +
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
  lewnard_plt, bias_plt, avg_eff_plt,
  nrow = 1,
  labels = "AUTO"
)

fig_dir <- here("plots", "supplemental_figs")
dir.create(fig_dir)

ggsave(
    here(fig_dir, "final_ve_compared_to_avg_vax_protect.png"),
    plt,
    width = 15,
    height = 5,
    units = "in",
    bg = "white"
)

ggsave(
    here(fig_dir, "final_ve_compared_to_avg_vax_protect.pdf"),
    plt,
    width = 15,
    height = 5,
    units = "in",
    bg = "white"
)
