library(here)
library(tidyverse)
library(cowplot)
library(paletteer)

source(here("src", "full_VE_model.R"))
source(here("src", "data_generation.R"))

# Helper function to draw a gamma distribution density curve for a given shape and mean
draw_gamma <- function(shape, mean, label) {
  return(ggplot2::geom_function(
        fun = dgamma,
        args = list(shape = shape, rate = shape / mean),
        xlim = c(0, 2),
        ggplot2::aes(color = factor(shape), linetype = label),
        linewidth = 1,
        n = 500
    ))
}

# Helper function to draw true vaccine protection label on plots
true_eff_label <- function(x, y, label, xoff = 55, yoff = 0.01) {
    box <- annotate(
        "rect",
        xmin = x - xoff,
        xmax = x + xoff,
        ymin = y - yoff,
        ymax = y + yoff,
        fill = "#ffffffaa"
    )

    text <- annotate(
        "text",
        x = x,
        y = y,
        label = label
    )

    return(list(box, text))
}

# Scenario options (no waning, continuous risk distributions, VE from cumulative attack rates)
opts = list(
    waning = FALSE,
    heterogeneity = TRUE,
    instantaneous = FALSE
)

# Parameters to draw cumulative attack rate curves
pars <- list(
    start_time = 0,
    end_time = 200,
    dt = 1,
    lambda = 0.0015,
    theta_0 = 1 - 0.3,
    epsilon_v = 1,
    epsilon_u = 1,
    alpha_v = c(0.2, 2, 20),
    alpha_u = c(0.2, 2, 20)
)

# Generate data for all cumulative attack rate curves
cAR_dt <- generate_par_sets(pars) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        cAR_v = map_vec(pars, function(x) cumulative_attack_rate(x$time, x, opts, TRUE)),
        cAR_u = map_vec(pars, function(x) cumulative_attack_rate(x$time, x, opts, FALSE))
    ) %>%
    unnest(pars)

# Parameters to draw final VE curves
pars <- list(
    start_time = 200,
    end_time = 200,
    dt = 0,
    lambda = 0.0015,
    theta_0 = 1 - seq(0.01, 0.99, 0.01),
    epsilon_v = 1,
    epsilon_u = 1,
    alpha_v = c(0.2, 2, 20),
    alpha_u = c(0.2, 2, 20)
)

# Generate data for final VE estimates
final_ve_dt <- generate_par_sets(pars, include_early = FALSE) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        estd_ve = map_vec(pars, function(x) estimate_math_ve(x, opts)),
        true_vax_prot = map_vec(pars, function(x) (1 - x$theta_0) * 100)
    ) %>%
    unnest(pars)

# Maintain consistent linetypes
lt <- c("Unvaccinated" = "solid", "Vaccinated" = "dashed",
        "gamma_unvax_inf" = "solid", "gamma_vax_inf" = "dashed")

# Maintain consistent line colors
cl <- c("0.2" = "#d95f02", "2" = "#1b9e77", "20" = "#7570b3")

### TOP ROW: vaccinated and unvaccinated populations have indentical distributions
# Panel A: distributions with differnet shape parameters
het_plot_a <- ggplot() +
    geom_vline(
        xintercept = 1,
        color = "gray50",
        linetype = "12",
        linewidth = 1
    ) +
    draw_gamma(shape = 0.2, mean = 1, label = "Unvaccinated") +
    draw_gamma(shape = 0.2, mean = 1, label = "Vaccinated") +
    draw_gamma(shape = 2, mean = 1, label = "Unvaccinated") +
    draw_gamma(shape = 2, mean = 1, label = "Vaccinated") +
    draw_gamma(shape = 20, mean = 1, label = "Unvaccinated") +
    draw_gamma(shape = 20, mean = 1, label = "Vaccinated") +
    scale_color_manual(name = "Shape parameter", values = cl) +
    scale_linetype_manual(name = element_blank(), values = lt) +
    ylim(0, 2) +
    theme_cowplot(14) +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.625, 0.675),
        legend.key.width = unit(35, "pt"),
        legend.spacing.y = unit(1, "pt")
    ) +
    labs(x = "Pre-vaccination susceptibility", y = "Density")

# Panel B: cumulative attack rates with shared distributions
het_plot_b <- cAR_dt %>%
    filter(alpha_v == alpha_u) %>%
    ggplot() +
        aes(x = time, group = interaction(alpha_u, alpha_v)) +
        geom_line(
            aes(y = cAR_u, color = factor(alpha_u), linetype = "Unvaccinated"),
            linewidth = 1
        ) +
        geom_line(
            aes(y = cAR_v, color = factor(alpha_v), linetype = "Vaccinated"),
            linewidth = 1
        ) +
        true_eff_label(x = 150, y = 0.02, label = "True vaccine protection = 30%") +
        scale_color_manual(values = cl) +
        scale_linetype_manual(values = lt) +
        coord_cartesian(xlim = c(0, 200), ylim = c(0, 0.3)) +
        scale_y_continuous(breaks = seq(0, 0.3, 0.1), labels = seq(0, 30, 10)) +
        theme_cowplot(14) +
        background_grid() +
        theme(legend.position = "none") +
        labs(x = "Time (days)", y = "Cumulative attack\nrate (%)")

# Panel C: final VE vs. true vaccine protection with shared distributions
het_plot_c <- final_ve_dt %>%
    filter(alpha_v == alpha_u) %>%
    ggplot() +
        aes(
            x = true_vax_prot,
            y = estd_ve,
            color = factor(alpha_u),
            group = interaction(alpha_u, alpha_v)
        ) +
        geom_hline(yintercept = 0, linetype = "44", color = "gray25") +
        geom_line(
            aes(y = true_vax_prot),
            linetype = "dashed",
            linewidth = 1,
            color = "black"
        ) +
        geom_line(linewidth = 1) +
        scale_color_manual(values = cl, guide = "none") +
        theme_cowplot() +
        background_grid() +
        theme(legend.position = "none") +
        labs(x = "True vaccine protection (%)", y = "Final VE (%)")

### MIDDLE ROW: unvaccinated distribution held constant (shape == 2)
# Panel D: vaccinated distributions with different shape parameters vs constant
#          unvaccinated distribution
het_plot_d <- ggplot() +
    geom_vline(
        xintercept = 1,
        color = "gray50",
        linetype = "12",
        linewidth = 1
    ) +
    draw_gamma(shape = 0.2, mean = 1, label = "Unvaccinated") +
    draw_gamma(shape = 2, mean = 1, label = "Unvaccinated") +
    draw_gamma(shape = 20, mean = 1, label = "Unvaccinated") +
    geom_function(
        fun = dgamma,
        args = list(shape = 2, rate = 2 / 1),
        xlim = c(0, 2),
        aes(linetype = "Vaccinated"),
        linewidth = 1,
        color = "black",
        n = 500
    ) +
    scale_color_manual(values = cl) +
    scale_linetype_manual(values = lt) +
    ylim(0, 2) +
    theme_cowplot(14) +
    theme(legend.position = "none") +
    labs(x = "Pre-vaccination susceptibility", y = "Density")

# Panel E: cumulative attack rates with constant unvaccinated distribution
het_plot_e <- cAR_dt %>%
    filter(alpha_v == 2) %>%
    ggplot() +
        aes(x = time, group = interaction(alpha_u, alpha_v)) +
        geom_line(
            aes(y = cAR_u, color = factor(alpha_u), linetype = "Unvaccinated"),
            linewidth = 1
        ) +
        geom_line(
            aes(y = cAR_v, linetype = "Vaccinated"),
            color = "black",
            linewidth = 1
        ) +
        true_eff_label(x = 150, y = 0.02, label = "True vaccine protection = 30%") +
        scale_color_manual(values = cl) +
        scale_linetype_manual(values = lt) +
        coord_cartesian(xlim = c(0, 200), ylim = c(0, 0.3)) +
        scale_y_continuous(breaks = seq(0, 0.3, 0.1), labels = seq(0, 30, 10)) +
        theme_cowplot(14) +
        background_grid() +
        theme(legend.position = "none") +
        labs(x = "Time (days)", y = "Cumulative attack\nrate (%)")

# Panel F: final VE vs. true vaccine protection with constant unvaccinated distribution
het_plot_f <- final_ve_dt %>%
    filter(alpha_v == 2) %>%
    ggplot() +
        aes(
            x = true_vax_prot,
            y = estd_ve,
            color = factor(alpha_u),
            group = interaction(alpha_u, alpha_v)
        ) +
        geom_hline(yintercept = 0, linetype = "44", color = "gray25") +
        geom_line(
            aes(y = true_vax_prot),
            linetype = "dashed",
            linewidth = 1,
            color = "black"
        ) +
        geom_line(linewidth = 1) +
        scale_color_manual(values = cl, guide = "none") +
        theme_cowplot() +
        background_grid() +
        theme(legend.position = "none") +
        labs(x = "True vaccine protection (%)", y = "Final VE (%)")

### BOTTOM ROW: vaccinated distribution held constant (shape == 2)
# Panel G: unvaccinated distributions with different shape parameters vs constant
#          vaccinated distribution
het_plot_g <- ggplot() +
    geom_vline(
        xintercept = 1,
        color = "gray50",
        linetype = "12",
        linewidth = 1
    ) +
    geom_function(
        fun = dgamma,
        args = list(shape = 2, rate = 2 / 1),
        xlim = c(0, 2),
        aes(linetype = "Unvaccinated"),
        linewidth = 1,
        color = "black",
        n = 500
    ) +
    draw_gamma(shape = 0.2, mean = 1, label = "Vaccinated") +
    draw_gamma(shape = 2, mean = 1, label = "Vaccinated") +
    draw_gamma(shape = 20, mean = 1, label = "Vaccinated") +
    scale_color_manual(values = cl) +
    scale_linetype_manual(values = lt) +
    ylim(0, 2) +
    theme_cowplot(14) +
    theme(legend.position = "none") +
    labs(x = "Pre-vaccination susceptibility", y = "Density")

# Panel H: cumulative attack rates with constant vaccinated distribution
het_plot_h <- cAR_dt %>%
    filter(alpha_u == 2) %>%
    ggplot() +
        aes(x = time, group = interaction(alpha_u, alpha_v)) +
        geom_line(
            aes(y = cAR_u, linetype = "Unvaccinated"),
            color = "black",
            linewidth = 1
        ) +
        geom_line(
            aes(y = cAR_v, color = factor(alpha_v), linetype = "Vaccinated"),
            linewidth = 1
        ) +
        true_eff_label(x = 150, y = 0.02, label = "True vaccine protection = 30%") +
        scale_color_manual(values = cl) +
        scale_linetype_manual(values = lt) +
        coord_cartesian(xlim = c(0, 200), ylim = c(0, 0.3)) +
        scale_y_continuous(breaks = seq(0, 0.3, 0.1), labels = seq(0, 30, 10)) +
        theme_cowplot(14) +
        background_grid() +
        theme(legend.position = "none") +
        labs(x = "Time (days)", y = "Cumulative attack\nrate (%)")

# Panel I: final VE vs. true vaccine protection with constant vaccinated distribution
het_plot_i <- final_ve_dt %>%
    filter(alpha_u == 2) %>%
    ggplot() +
        aes(
            x = true_vax_prot,
            y = estd_ve,
            color = factor(alpha_v),
            group = interaction(alpha_u, alpha_v)
        ) +
        geom_hline(yintercept = 0, linetype = "44", color = "gray25") +
        geom_line(
            aes(y = true_vax_prot),
            linetype = "dashed",
            linewidth = 1,
            color = "black"
        ) +
        geom_line(linewidth = 1) +
        scale_color_manual(values = cl, guide = "none") +
        theme_cowplot() +
        background_grid() +
        theme(legend.position = "none") +
        labs(x = "True vaccine protection (%)", y = "Final VE (%)")

# Combine panels into one plot and save
het_curve_plot <- plot_grid(
    het_plot_a, het_plot_b, het_plot_c,
    het_plot_d, het_plot_e, het_plot_f,
    het_plot_g, het_plot_h, het_plot_i,
    byrow = TRUE,
    labels = "AUTO",
    align = "hv",
    axis = "lb",
    ncol = 3,
    nrow = 3
)

fig_path <- here("plots")
dir.create(fig_path)

ggsave(
    here(fig_path, "fig2.png"),
    het_curve_plot,
    bg = "white",
    height = 9,
    width = 15,
    units = "in"
)

ggsave(
    here(fig_path, "fig2.pdf"),
    het_curve_plot,
    bg = "white",
    height = 9,
    width = 15,
    units = "in"
)