library(here)
library(tidyverse)
library(survival)
library(cowplot)
library(paletteer)

source(here("src", "full_VE_model.R"))
source(here("src", "data_generation.R"))

set.seed(0)

# Scenario options (no waning, continuous pre-vaccination risk, VE from cumulative attack rates)
opts = list(
    waning = FALSE,
    heterogeneity = TRUE,
    instantaneous = FALSE
)
cumul_ve_opts <- opts

# Create copy of options for estimating VE from instantaneous incidence rates
insnt_ve_opts <- opts
insnt_ve_opts$instantaneous <- TRUE

# Parameters to draw final VE curves
pars = list(
    start_time = 0,
    end_time = 200,
    dt = 14,
    lambda = 0.0015,
    lambda_negative = 0.0045,
    theta_0 = 1 - seq(0.01, 0.99, 0.01),
    epsilon_v = 1,
    epsilon_u = 1,
    alpha_v = c(0.2, 2, 20),
    alpha_u = c(0.2, 2, 20),
    mu_v = 1,
    mu_u = 1,
    pi_pos = 1,
    pi_neg = 1,
    c = 0.5,
    N = 5e3
)

pars$time <- pars$end_time

# Generate parameter sets
pars_dt <- crossing(!!!pars) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup()

# Generate VE data (analytical and logistic regression models)
final_ve_dt <- pars_dt %>%
    rowwise_mutate(counts = list(generate_counts(pars, opts))) %>%
    mutate(linelist = map(counts, function(x) generate_linelist(x))) %>%
    mutate(
        insnt_ve = map_vec(pars, function(x) estimate_math_ve(x, insnt_ve_opts)),
        cumul_ve = map_vec(pars, function(x) estimate_math_ve(x, cumul_ve_opts)),
        ulreg_ve = map_vec(linelist, function(x) estimate_ve_uncond_logreg(x)),
        clreg_ve = map_vec(linelist, function(x) estimate_ve_cond_logreg(x))
    ) %>%
    unnest(pars) %>%
    select(theta_0, starts_with("alpha_"), ends_with("_ve")) %>%
    pivot_longer(ends_with("_ve"), names_to = "method", values_to = "estd_ve") %>%
    mutate(
        true_vax_protect = (1 - theta_0) * 100,
        method_f = factor(
            method,
            levels = c("cumul_ve", "ulreg_ve", "clreg_ve", "insnt_ve"),
            labels = c("VE cumulative", "VE unconditional", "VE conditional", "VE instantaneous")
        )
    )

ve_comp <- final_ve_dt %>%
    select(!method_f) %>%
    pivot_wider(names_from = method, values_from = estd_ve) %>%
    select(ends_with("_ve")) %>%
    pivot_longer(!cumul_ve) %>%
    mutate(ve_diff = value - cumul_ve) %>%
    group_by(name) %>%
    arrange(ve_diff) %>%
    summarize(
        mean_diff = mean(ve_diff),
        lower_diff = quantile(ve_diff, probs = c(0.025)),
        upper_diff = quantile(ve_diff, probs = c(0.925))
    )

print("Difference in VE rel. to VE^cumulative (mean, 95% IQR)")
print(ve_comp)

# Maintain consistent linetypes
lt <- c("Unvaccinated" = "solid", "Vaccinated" = "dashed",
        "gamma_unvax_inf" = "solid", "gamma_vax_inf" = "dashed")

# Maintain consistent line colors
cl <- c("0.2" = "#d95f02", "2" = "#1b9e77", "20" = "#7570b3")

# Heterogeneity legend labels
labs <- c(
    "0.2" = expression(paste("High (", alpha, " = 0.2)")),
    "2" = expression(paste("Moderate (", alpha, " = 2)")),
    "20" = expression(paste("Low (", alpha, " = 20)"))
)

top <- final_ve_dt %>%
    filter(alpha_u == alpha_v) %>%
    ggplot() +
        aes(x = true_vax_protect, y = estd_ve, color = factor(alpha_u)) +
        geom_hline(yintercept = 0, linetype = "44", color = "gray25") +
        geom_line(
            aes(y = true_vax_protect),
            linetype = "dashed",
            linewidth = 1,
            color = "black"
        ) +
        geom_point() +
        scale_color_manual(
            name = "Heterogeneity",
            values = cl,
            labels = labs
        ) +
        facet_grid(cols = vars(method_f)) +
        theme_cowplot() +
        background_grid() +
        theme(
            legend.position = "inside",
            legend.position.inside = c(0.01, 0.75),
            legend.key.width = unit(35, "pt"),
            legend.spacing.y = unit(1, "pt")
        ) +
        labs(x = "True vaccine protection (%)", y = "Final VE (%)")

mid <- final_ve_dt %>%
    filter(alpha_v == 20, alpha_u != 20) %>%
    ggplot() +
        aes(x = true_vax_protect, y = estd_ve, color = factor(alpha_u)) +
        geom_hline(yintercept = 0, linetype = "44", color = "gray25") +
        geom_line(
            aes(y = true_vax_protect),
            linetype = "dashed",
            linewidth = 1,
            color = "black"
        ) +
        geom_point() +
        scale_color_manual(
            name = "Heterogeneity",
            values = cl,
            labels = labs
        ) +
        facet_grid(cols = vars(method_f)) +
        theme_cowplot() +
        background_grid() +
        theme(legend.position = "none") +
        labs(x = "True vaccine protection (%)", y = "Final VE (%)")

bot <- final_ve_dt %>%
    filter(alpha_u == 20, alpha_v != 20) %>%
    ggplot() +
        aes(x = true_vax_protect, y = estd_ve, color = factor(alpha_v)) +
        geom_hline(yintercept = 0, linetype = "44", color = "gray25") +
        geom_line(
            aes(y = true_vax_protect),
            linetype = "dashed",
            linewidth = 1,
            color = "black"
        ) +
        geom_point() +
        scale_color_manual(
            name = "Heterogeneity",
            values = cl,
            labels = labs
        ) +
        facet_grid(cols = vars(method_f)) +
        theme_cowplot() +
        background_grid() +
        theme(legend.position = "none") +
        labs(x = "True vaccine protection (%)", y = "Final VE (%)")

plt <- plot_grid(
    top, mid, bot,
    byrow = TRUE,
    align = "hv",
    axis = "lb",
    ncol = 1,
    nrow = 3
)

fig_path <- here("plots", "supplemental_figs")
dir.create(fig_path)

ggsave(
    here(fig_path, "final_ve_curves_across_methods_varying_heterogeneity.png"),
    plt,
    bg = "white",
    height = 9,
    width = 15,
    units = "in"
)

ggsave(
    here(fig_path, "final_ve_curves_across_methods_varying_heterogeneity.pdf"),
    plt,
    bg = "white",
    height = 9,
    width = 15,
    units = "in"
)
