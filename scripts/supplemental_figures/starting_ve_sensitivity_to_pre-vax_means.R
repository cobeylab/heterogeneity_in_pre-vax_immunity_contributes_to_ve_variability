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

# VE is constrained to (-infinity, 100]. As a result, negative bias can be much
# larger than positive bias for similar levels of mean pre-vaccination risk. This
# helper function balances the fill color scale to be more equal from 0 to the
# maximum and minimum VE observed in the dataset.
clean_fill <- function(x, xmin = NA, xmax = NA) {
    min <- ifelse(is.na(xmin), abs(min(x)), xmin)
    max <- ifelse(is.na(xmax), abs(max(x)), xmax)

    norm_x <- rep(0, length(x))
    for (i in 1:length(x)) {
        val <- x[i]
        if (val < 0) {
            norm_x[i] <- val / min
        } else if (val > 0) {
            norm_x[i] <- val / max
        }
    }

    return(norm_x)
}

# Scenario options (no waning, continuous risk distributions, VE from cumulative attack rates)
opts = list(
    waning = FALSE,
    heterogeneity = TRUE,
    instantaneous = FALSE
)

# Heatmap scenario parameters (only mean pre-vaccination risk varies)
pars <- list(
    start_time = 200,
    end_time = 200,
    dt = 0,
    lambda = 0.005,
    theta_0 = 1 - seq(0.1, 0.9, 0.1),
    epsilon_v = seq(0.1, 1.0, 0.1),
    epsilon_u = seq(0.1, 1.0, 0.1),
    alpha_v = 20,
    alpha_u = 20
)

# Calculate bias between starting VE and true vaccine protection for each scenario
heatmap_dt <- generate_par_sets(pars, include_early = FALSE) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        true_vax_protection = map_vec(pars, function(x) (1 - x$theta_0) * 100),
        starting_ve = map_vec(pars, function(x) estimate_ve_starting(x)),
        starting_ve_bias = starting_ve - true_vax_protection,
        cell_label = clean_label(starting_ve_bias, sigdig = 2)
    ) %>%
    unnest(pars)

# Cleaning up color fill values for better heatmap presentation
heatmap_dt$fill_val <- clean_fill(heatmap_dt$starting_ve_bias, xmin = NA, xmax = 100)
fill_lim <- round(max(
    abs(min(heatmap_dt$fill_val)),
    abs(max(heatmap_dt$fill_val))
)) * 1.01

# Heatmap of starting VE bias across mean pre-vaccination risk values and true vaccine protection
heatmap_plot <- ggplot(heatmap_dt) +
    aes(
        x = epsilon_u,
        y = epsilon_v,
        fill = fill_val,
        label = cell_label
    ) +
    geom_tile(color = NA) +
    geom_text() +
    facet_wrap(
        vars(true_vax_protection),
        labeller = labeller(
            true_vax_protection = function(x) paste0("True vaccine protection = ", x, "%")
        )
    ) +
    scale_fill_paletteer_c(
        "ggthemes::Orange-Blue Diverging",
        name = expression("Difference between starting" ~ widehat(VE) ~ "and true vaccine protection (% pts)"),
        limits = c(-fill_lim, fill_lim),
        breaks = c(0)
    ) +
    coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05), expand = FALSE) +
    theme_cowplot(14) +
    theme(
        legend.position = "top",
        legend.background = element_rect(fill = "white")
    ) +
    labs(
        x = "Mean unvaccinated pre-vaccination susceptibility",
        y = "Mean vaccinated pre-vaccination susceptibility"
    )

fig_path <- here("plots", "supplemental_figs")
dir.create(fig_path)

ggsave(
  here(fig_path, "starting_ve_vs_mean_pre-vax_suscep.png"),
  heatmap_plot,
  bg = "white",
  height = 13,
  width = 13,
  units = "in"
)

ggsave(
  here(fig_path, "starting_ve_vs_mean_pre-vax_suscep.pdf"),
  heatmap_plot,
  bg = "white",
  height = 13,
  width = 13,
  units = "in"
)