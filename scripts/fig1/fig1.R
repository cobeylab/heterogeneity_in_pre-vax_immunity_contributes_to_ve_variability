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

### Panel A
# Heatmap scenario parameters (only mean pre-vaccination risk varies)
pars <- list(
    start_time = 200,
    end_time = 200,
    dt = 0,
    lambda = 0.005,
    theta_0 = 1 - 0.5,
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
        starting_ve = map_vec(pars, function(x) estimate_ve_starting(x) * 100),
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

# Cells to highlight for VE curves.
cell_highlight <- tibble(
    epsilon_u = seq(0.3, 0.9, 0.1),
    epsilon_v = seq(0.9, 0.3, -0.1),
)

# Panel A: Heatmap of starting VE bias across mean pre-vaccination risk values
heatmap_plot <- ggplot(heatmap_dt) +
    aes(
        x = epsilon_u,
        y = epsilon_v,
        fill = fill_val,
        label = cell_label
    ) +
    geom_tile(color = NA) +
    geom_tile(
        data = cell_highlight,
        aes(x = epsilon_u, y = epsilon_v),
        fill = NA,
        color = "gold",
        linewidth = 1,
        inherit.aes = FALSE
    ) +
    geom_text() +
    facet_grid(
        cols = vars(true_vax_protection),
        labeller = labeller(
            true_vax_protection = function(x) paste0("True vaccine protection = ", x, "%")
        )
    ) +
    scale_fill_paletteer_c(
        "ggthemes::Orange-Blue Diverging",
        name = expression(atop(
            'Difference between starting' ~ widehat(VE),
            'and true vaccine protection (% pts)'
        )),
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

### Panel B
# VE curve scenario parameters for each highlighted cell in the heatmap
pars <- list(
    start_time = 0,
    end_time = 400,
    dt = 1,
    lambda = 0.05,
    theta_0 = 1 - 0.5,
    epsilon_v = 0.3,
    epsilon_u = 0.9,
    alpha_v = 20,
    alpha_u = 20
)

curves_dt <- generate_par_sets(pars)

for (i in 1:6) {
    mean_modifier <- 0.1 * i
    pars <- list(
        start_time = 0,
        end_time = 400,
        dt = 1,
        lambda = 0.05,
        theta_0 = 1 - 0.5,
        epsilon_v = 0.3 + mean_modifier,
        epsilon_u = 0.9 - mean_modifier,
        alpha_v = 20,
        alpha_u = 20
    )

    curves_dt <- curves_dt %>%
        bind_rows(generate_par_sets(pars))
}

# Calculate cumulative VE over time for each scenario
curves_dt <- curves_dt %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        true_vax_protection = map_vec(pars, function(x) (1 - x$theta_0) * 100),
        starting_ve = map_vec(pars, function(x) estimate_ve_starting(x) * 100),
        estd_ve = map_vec(pars, function(x) estimate_math_ve(x, opts)),
        starting_ve_bias = starting_ve - true_vax_protection,
        cell_label = clean_label(starting_ve_bias, sigdig = 2)
    ) %>%
    unnest(pars)

# Match the line colors to the highlighted cells' fill colors (i.e., starting VE bias)
curves_dt$line_col <- clean_fill(
    x = curves_dt$starting_ve_bias,
    xmin = abs(min(heatmap_dt$starting_ve_bias)),
    xmax = 100
)

# Panel B: VE over time for each highlighted scenario in the heatmap
curves_plot <- ggplot(curves_dt) +
    aes(
        x = time,
        y = estd_ve,
        color = line_col,
        group = interaction(epsilon_u, epsilon_v)
    ) +
    geom_hline(
        linewidth = 1,
        yintercept = 50,
        linetype = "22",
        color = "gray50"
    ) +
    geom_line(linewidth = 1) +
    annotate(
        "text",
        x = 200,
        y = 55,
        label = "True vaccine protection"
    ) +
    scale_color_paletteer_c(
        "ggthemes::Orange-Blue Diverging",
        name = "Difference between starting\nVE and true vaccine protection (% pts)",
        limits = c(-fill_lim, fill_lim),
        breaks = c(0)
    ) +
    coord_cartesian(xlim = c(0, 405), ylim = c(-50, 100), expand = FALSE) +
    theme_cowplot() +
    background_grid(minor = "y") +
    theme(legend.position = "none") +
    labs(x = "Time (days)", y = "Estimated vaccine effectiveness (%)")

# Combine panels into one plot and save
combo <- plot_grid(
    heatmap_plot, curves_plot,
    labels = "AUTO"
    )

fig_path <- here("plots")
dir.create(fig_path)

ggsave(
  here(fig_path, "fig1.png"),
  combo,
  bg = "white",
  height = 5,
  width = 10,
  units = "in"
)

ggsave(
  here(fig_path, "fig1.pdf"),
  combo,
  bg = "white",
  height = 5,
  width = 10,
  units = "in"
)