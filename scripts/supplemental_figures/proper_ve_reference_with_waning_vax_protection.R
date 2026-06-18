library(here)
library(tidyverse)
library(cowplot)
library(paletteer)

source(here("src", "full_VE_model.R"))
source(here("src", "data_generation.R"))

# Scenario options (no waning, continuous risk distributions, VE from cumulative attack rates)
opts = list(
    waning = TRUE,
    heterogeneity = FALSE,
    instantaneous = FALSE
)

# Heatmap scenario parameters (only mean pre-vaccination risk varies)
pars <- list(
    start_time = 0,
    end_time = 200,
    dt = 1,
    lambda = 0.005,
    theta_0 = 1 - 0.6,
    epsilon_v = 1,
    epsilon_u = 1,
    alpha_v = Inf,
    alpha_u = Inf,
    eta = c(30, 180, 720, Inf)
)

# Generate hazard and VE data
dt <- generate_par_sets(pars) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup() %>%
    mutate(
        vax_effect = map_vec(pars, function(x) vaccine_direct_effect(x$time, x, opts)),
        tot_vax_effect = map_vec(pars, function(x) total_vaccine_effect(x$time, x, opts)),
        estd_ve = map_vec(pars, function(x) estimate_math_ve(x, opts))
    ) %>%
    unnest(pars) %>%
    mutate(
        unvax_hazard = lambda * epsilon_u,
        tot_unvax_hazard = unvax_hazard * time,
        vax_hazard = lambda * epsilon_v * vax_effect,
        tot_vax_hazard = lambda * epsilon_v * tot_vax_effect,
        hazard_ratio = (1 - (tot_vax_hazard / tot_unvax_hazard)) * 100,
        avg_vax_protect = (1 - (tot_vax_effect / time)) * 100
    ) %>%
    select(eta, time, estd_ve, hazard_ratio, avg_vax_protect, ends_with("_hazard")) %>%
    pivot_longer(!c(eta, time))

# Consistent figure aesthetics
lt <- c(
  "unvax_hazard" = "solid",
  "vax_hazard" = "22",
  "estd_ve" = "solid",
  "hazard_ratio" = "22",
  "avg_vax_protect" = "solid"
)

cl <- c(
  "unvax_hazard" = "darkorange",
  "vax_hazard" = "dodgerblue",
  "estd_ve" = "#d95f02",
  "hazard_ratio" = "gray50",
  "avg_vax_protect" = "black"
)

fl <- c(
  "unvax_hazard" = "darkorange",
  "vax_hazard" = "dodgerblue",
  "estd_ve" = "#00000000",
  "hazard_ratio" = "#00000000",
  "avg_vax_protect" = "#00000000"
)

labs <- c(
  "unvax_hazard" = "Unvaccinated",
  "vax_hazard" = "Vaccinated",
  "estd_ve" = "VE estimate",
  "hazard_ratio" = "1 - Total hazard ratio",
  "avg_vax_protect" = "True average vaccine protection"
)

eta_labs <- c(
    "Inf" = "No vaccine waning",
    "720" = "Waning half life of 2 years",
    "180" = "Waning half life of 6 months",
    "30" = "Waning half life of 30 days"
)

common_plt <- list(
    scale_linetype_manual(
        name = NULL,
        values = lt,
        labels = labs
    ),
    scale_color_manual(
        name = NULL,
        values = cl,
        labels = labs
    ),
    theme_cowplot(),
    background_grid()
)

# Panel A: vaccinated and unvaccinated time-varying infeciton hazard vs. waning
hazard_names <- c("unvax_hazard", "vax_hazard")
hazard_plt <- ggplot(dt %>% filter(name %in% hazard_names)) +
    aes(
        x = time,
        y = value,
        ymax = value,
        ymin = 0,
        linetype = factor(name, levels = hazard_names),
        color = factor(name, levels = hazard_names),
        fill = factor(name, levels = hazard_names)
    ) +
    facet_grid(cols = vars(factor(
        eta,
        levels = c(Inf, 720, 180, 30),
        labels = eta_labs
    ))) +
    geom_ribbon(alpha = 0.5, color = NA) +
    geom_line(linewidth = 2) +
    common_plt +
    scale_fill_manual(
        name = NULL,
        values = fl,
        labels = labs
    ) +
    coord_cartesian(expand = FALSE, ylim = c(0, 0.0051)) +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.825, 0.15),
        legend.background = element_rect(fill = "#ffffff99"),
        legend.key.spacing.y = unit(0.01, "npc"),
        legend.key.width = unit(0.075, "npc"),
        legend.margin = margin(1, 1, 1, 1)
    ) +
    labs(x = "Time (days)", y = "Hazard")

# Panel B: VE estimates vs two reference values (hazard ratio, average vax protection)
ve_names <- c("avg_vax_protect", "hazard_ratio", "estd_ve")
ve_plt <- ggplot(dt %>% filter(name %in% ve_names)) +
    aes(
        x = time,
        y = value,
        linetype = factor(name, levels = ve_names),
        color = factor(name, levels = ve_names)
    ) +
    geom_line(linewidth = 2) +
    facet_grid(cols = vars(factor(
        eta,
        levels = c(Inf, 720, 180, 30),
        labels = eta_labs
    ))) +
    coord_cartesian(ylim = c(0, NA)) +
    common_plt +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.015, 0.215),
        legend.background = element_rect(fill = "#ffffff99"),
        legend.key.spacing.y = unit(0.01, "npc"),
        legend.key.width = unit(0.075, "npc")
    ) +
    labs(x = "Time (days)", y = "Vaccine effectiveness (%)")

# Combine panels and save plot
plt <- plot_grid(
  hazard_plt, ve_plt,
  ncol = 1,
  byrow = TRUE,
  labels = "AUTO",
  align = "hv",
  axis = "lb"
)

fig_dir <- here("plots", "supplemental_figs")
dir.create(fig_dir)

ggsave(
  here(fig_dir, "reference_for_ve_with_waning_vax_protection.png"),
  plt,
  width = 10,
  height = 6,
  units = "in",
  bg = "white"
)

ggsave(
  here(fig_dir, "reference_for_ve_with_waning_vax_protection.pdf"),
  plt,
  width = 10,
  height = 6,
  units = "in",
  bg = "white"
)