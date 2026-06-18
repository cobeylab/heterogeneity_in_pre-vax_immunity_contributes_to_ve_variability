library(here)
library(tidyverse)
library(survival)
library(cowplot)
library(paletteer)

source(here("src", "full_VE_model.R"))
source(here("src", "data_generation.R"))

facet_labels = c(
    insnt_ve = "VE instantaneous",
    ulreg_ve = "VE unconditional",
    clreg_ve = "VE conditional"
)

### Final ve sweeping across infection hazard, mean pre-vaccination risk, and
### true vaccine protection (NO WANING)

set.seed(0)

# Scenario parameters
pars = list(
    start_time = 0,
    end_time = 200,
    dt = 14,
    lambda = c(0.0015, 0.005, 0.01),
    lambda_negative = 0.02,
    theta_0 = 1 - seq(0.1, 0.9, 0.1),
    eta = -1,
    epsilon_v = seq(0.1, 1.0, 0.1),
    epsilon_u = seq(0.1, 1.0, 0.1),
    alpha_v = 20,
    alpha_u = 20,
    mu_v = 1,
    mu_u = 1,
    pi_pos = 1,
    pi_neg = 1,
    c = 0.5,
    N = 1e3
)

pars$time <- pars$end_time

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

# Generate parameter sets
pars_dt <- crossing(!!!pars) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup()

print(paste0("generating ", nrow(pars_dt)," no waning scenario linelists"))

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
    unnest(pars)

print("generating no waning scenario plot")

# Calculate adjusted R-squared for each model comparison
no_waning_insnt_vs_cumul_lm <- lm(insnt_ve ~ cumul_ve, data = final_ve_dt)
no_waning_ulreg_vs_cumul_lm <- lm(ulreg_ve ~ cumul_ve, data = final_ve_dt)
no_waning_clreg_vs_cumul_lm <- lm(clreg_ve ~ cumul_ve, data = final_ve_dt)

no_waning_insnt_vs_cumul_adj_rsq <- summary(no_waning_insnt_vs_cumul_lm)$adj.r.squared
no_waning_ulreg_vs_cumul_adj_rsq <- summary(no_waning_ulreg_vs_cumul_lm)$adj.r.squared
no_waning_clreg_vs_cumul_adj_rsq <- summary(no_waning_clreg_vs_cumul_lm)$adj.r.squared

# Calculate absolute average difference between VE estimates
no_waning_ve_compare <- final_ve_dt %>%
    select(ends_with("_ve")) %>%
    pivot_longer(!cumul_ve) %>%
    group_by(name) %>%
    summarize(mean_diff = mean(abs(value - cumul_ve))) %>%
    ungroup() %>%
    arrange(name) %>%
    mutate(adj_rsq = c(no_waning_clreg_vs_cumul_adj_rsq,
                       no_waning_insnt_vs_cumul_adj_rsq,
                       no_waning_ulreg_vs_cumul_adj_rsq))

# Scatterplots comparing VE estimates from different models across scenarios
no_waning_plt <- final_ve_dt %>%
    pivot_longer(c(insnt_ve, ulreg_ve, clreg_ve)) %>%
    ggplot() +
        aes(x = cumul_ve, y = value) +
        geom_hline(yintercept = 0, linetype = "44", color = "gray50") +
        geom_vline(xintercept = 0, linetype = "44", color = "gray50") +
        geom_line(aes(y = cumul_ve), color = "black", linewidth = 1, linetype = "44") +
        geom_point(size = 2) +
        geom_smooth(
            method = "lm",
            linewidth = 2,
            color = "#ff9f9faa"
        ) +
        geom_text(
            data = no_waning_ve_compare,
            aes(
                x = -475,
                y = 100,
                label = paste0("Adj. r-squared: ", signif(adj_rsq, digits = 3))
            ),
            inherit.aes = FALSE
        ) +
        geom_text(
            data = no_waning_ve_compare,
            aes(
                x = -475,
                y = 85,
                label = paste0("Mean abs. difference: ", signif(mean_diff, digits = 3))
            ),
            inherit.aes = FALSE
        ) +
        facet_grid(cols = vars(factor(name, labels = facet_labels))) +
        theme_cowplot(20) +
        background_grid() +
        theme(legend.position = "top") +
        labs(x = "VE cumulative (%)", y = "Alternative VE estimate (%)")

rm(pars_dt, final_ve_dt)

### Final ve sweeping across infection hazard, mean pre-vaccination risk, and
### true vaccine protection (WITH WANING)

set.seed(0)

# Scenario parameters
pars = list(
    start_time = 0,
    end_time = 200,
    dt = 14,
    lambda = c(0.0015, 0.005, 0.01),
    lambda_negative = 0.02,
    theta_0 = 1 - c(0.3, 0.6, 0.9),
    eta = c(30, 180, 360, 1440),
    epsilon_v = seq(0.1, 1.0, 0.1),
    epsilon_u = seq(0.1, 1.0, 0.1),
    alpha_v = 20,
    alpha_u = 20,
    mu_v = 1,
    mu_u = 1,
    pi_pos = 1,
    pi_neg = 1,
    c = 0.5,
    N = 1e3
)

pars$time <- pars$end_time

# Scenario options (waning, continuous pre-vaccination risk, VE from cumulative attack rates)
opts = list(
    waning = TRUE,
    heterogeneity = TRUE,
    instantaneous = FALSE
)
cumul_ve_opts <- opts

# Create copy of options for estimating VE from instantaneous incidence rates
insnt_ve_opts <- opts
insnt_ve_opts$instantaneous <- TRUE

# Generate parameter sets
pars_dt <- crossing(!!!pars) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup()

print(paste0("generating ", nrow(pars_dt)," waning scenario linelists"))

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
    unnest(pars)

print("generating waning scenario plot")

# Calculate adjusted R-squared for each model comparison
waning_insnt_vs_cumul_lm <- lm(insnt_ve ~ cumul_ve, data = final_ve_dt)
waning_ulreg_vs_cumul_lm <- lm(ulreg_ve ~ cumul_ve, data = final_ve_dt)
waning_clreg_vs_cumul_lm <- lm(clreg_ve ~ cumul_ve, data = final_ve_dt)

waning_insnt_vs_cumul_adj_rsq <- summary(waning_insnt_vs_cumul_lm)$adj.r.squared
waning_ulreg_vs_cumul_adj_rsq <- summary(waning_ulreg_vs_cumul_lm)$adj.r.squared
waning_clreg_vs_cumul_adj_rsq <- summary(waning_clreg_vs_cumul_lm)$adj.r.squared

# Calculate absolute average difference between VE estimates
waning_ve_compare <- final_ve_dt %>%
    select(ends_with("_ve")) %>%
    pivot_longer(!cumul_ve) %>%
    group_by(name) %>%
    summarize(mean_diff = mean(abs(value - cumul_ve))) %>%
    ungroup() %>%
    arrange(name) %>%
    mutate(adj_rsq = c(waning_clreg_vs_cumul_adj_rsq,
                       waning_insnt_vs_cumul_adj_rsq,
                       waning_ulreg_vs_cumul_adj_rsq))

# Scatterplots comparing VE estimates from different models across scenarios
waning_plt <- final_ve_dt %>%
    pivot_longer(c(insnt_ve, ulreg_ve, clreg_ve)) %>%
    ggplot() +
        aes(x = cumul_ve, y = value) +
        geom_hline(yintercept = 0, linetype = "44", color = "gray50") +
        geom_vline(xintercept = 0, linetype = "44", color = "gray50") +
        geom_line(aes(y = cumul_ve), color = "black", linewidth = 1, linetype = "44") +
        geom_point(size = 2) +
        geom_smooth(
            method = "lm",
            linewidth = 2,
            color = "#ff9f9faa"
        ) +
        geom_text(
            data = waning_ve_compare,
            aes(
                x = -475,
                y = 100,
                label = paste0("Adj. r-squared: ", signif(adj_rsq, digits = 3))
            ),
            inherit.aes = FALSE
        ) +
        geom_text(
            data = waning_ve_compare,
            aes(
                x = -475,
                y = 85,
                label = paste0("Mean abs. difference: ", signif(mean_diff, digits = 3))
            ),
            inherit.aes = FALSE
        ) +
        facet_grid(cols = vars(factor(name, labels = facet_labels))) +
        theme_cowplot(20) +
        background_grid() +
        theme(legend.position = "top") +
        labs(x = "VE cumulative (%)", y = "Alternative VE estimate (%)")

rm(pars_dt, final_ve_dt)

### Final ve sweeping across infection hazard, and pre-vaccination risk heterogeneity (NO WANING)

set.seed(0)

# Scenario parameters
pars = list(
    start_time = 0,
    end_time = 200,
    dt = 14,
    lambda = 0.0015,
    lambda_negative = 0.0045,
    theta_0 = 1 - seq(0.1, 0.9, 0.1),
    eta = -1,
    epsilon_v = 1.0,
    epsilon_u = 1.0,
    alpha_v = seq(0.2, 2, 0.2),
    alpha_u = seq(0.2, 2, 0.2),
    mu_v = 1,
    mu_u = 1,
    pi_pos = 1,
    pi_neg = 1,
    c = 0.5,
    N = 1e3
)

pars$time <- pars$end_time

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

# Generate parameter sets
pars_dt <- crossing(!!!pars) %>%
    mutate(par_set_id = row_number()) %>%
    group_by(par_set_id) %>%
    nest(pars = -group_cols()) %>%
    ungroup()

print(paste0("generating ", nrow(pars_dt)," heterogeneity scenario linelists"))

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
    unnest(pars)

print("generating heterogeneity scenario plot")

# Calculate adjusted R-squared for each model comparison
het_insnt_vs_cumul_lm <- lm(insnt_ve ~ cumul_ve, data = final_ve_dt)
het_ulreg_vs_cumul_lm <- lm(ulreg_ve ~ cumul_ve, data = final_ve_dt)
het_clreg_vs_cumul_lm <- lm(clreg_ve ~ cumul_ve, data = final_ve_dt)

het_insnt_vs_cumul_adj_rsq <- summary(het_insnt_vs_cumul_lm)$adj.r.squared
het_ulreg_vs_cumul_adj_rsq <- summary(het_ulreg_vs_cumul_lm)$adj.r.squared
het_clreg_vs_cumul_adj_rsq <- summary(het_clreg_vs_cumul_lm)$adj.r.squared

# Calculate absolute average difference between VE estimates
het_ve_compare <- final_ve_dt %>%
    select(ends_with("_ve")) %>%
    pivot_longer(!cumul_ve) %>%
    group_by(name) %>%
    summarize(mean_diff = mean(abs(value - cumul_ve))) %>%
    ungroup() %>%
    arrange(name) %>%
    mutate(adj_rsq = c(het_clreg_vs_cumul_adj_rsq,
                       het_insnt_vs_cumul_adj_rsq,
                       het_ulreg_vs_cumul_adj_rsq))

# Scatterplots comparing VE estimates from different models across scenarios
het_plt <- final_ve_dt %>%
    pivot_longer(c(insnt_ve, ulreg_ve, clreg_ve)) %>%
    ggplot() +
        aes(x = cumul_ve, y = value) +
        geom_hline(yintercept = 0, linetype = "44", color = "gray50") +
        geom_vline(xintercept = 0, linetype = "44", color = "gray50") +
        geom_line(aes(y = cumul_ve), color = "black", linewidth = 1, linetype = "44") +
        geom_point(size = 2) +
        geom_smooth(
            method = "lm",
            linewidth = 2,
            color = "#ff9f9faa"
        ) +
        geom_text(
            data = het_ve_compare,
            aes(
                x = 0,
                y = 90,
                label = paste0("Adj. r-squared: ", signif(adj_rsq, digits = 3))
            ),
            inherit.aes = FALSE
        ) +
        geom_text(
            data = het_ve_compare,
            aes(
                x = 0,
                y = 85,
                label = paste0("Mean abs. difference: ", signif(mean_diff, digits = 3))
            ),
            inherit.aes = FALSE
        ) +
        facet_grid(cols = vars(factor(name, labels = facet_labels))) +
        theme_cowplot(20) +
        background_grid() +
        theme(legend.position = "top") +
        labs(x = "VE cumulative (%)", y = "Alternative VE estimate (%)")

rm(pars_dt, final_ve_dt)

# plot output
fig_path <- here("plots", "supplemental_figs")
dir.create(fig_path)

ggsave(
    here(fig_path, "ve_vs_mean_pre-vax_suscep_wo_waning.png"),
    no_waning_plt,
    bg = "white",
    height = 9,
    width = 12,
    units = "in"
)

ggsave(
    here(fig_path, "ve_vs_mean_pre-vax_suscep_w_waning.png"),
    waning_plt,
    bg = "white",
    height = 9,
    width = 12,
    units = "in"
)

ggsave(
    here(fig_path, "ve_vs_pre-vax_suscep_heterogeneity.png"),
    het_plt,
    bg = "white",
    height = 9,
    width = 12,
    units = "in"
)