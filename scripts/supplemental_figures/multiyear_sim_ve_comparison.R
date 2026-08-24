library(here)
library(tidyverse)
library(optparse)
library(RcppTOML)
library(survival)
library(cowplot)
library(paletteer)

source(here("src", "data_generation.R"))

# maintain consistent RNG seed
set.seed(0)

# parse command line arguments
if (interactive()) {
    parser <- OptionParser()
    parser <- add_option(parser, "--config", help = "path to config TOML file (rel path from this script)")
    parsed_args <- parse_args(parser, args = c("--config=fig3/fig3.toml"))
} else {
    parser <- OptionParser()
    parser <- add_option(parser, "--config", help = "path to config TOML file (rel path from this script)")
    parsed_args <- parse_args(parser)
}

if (is.null(parsed_args[["config"]])) {
  cat("ERROR: --config is required\n")
  quit("no", status = 1, runLast = FALSE)
}

# parse TOML parameter file
script_dir <- here("scripts")
config_path <- here(script_dir, parsed_args[["config"]])

config <- parseTOML(config_path)
pars <- config$parameters
exp_name <- pars$exp_name

# read in simulation data
print("READING DATA")

dat_dir <- here("data", exp_name)
exp_filepath <- here(dat_dir, "experiments.csv")
res_filepath <- here(dat_dir, "results.csv")
# ll_filepath <- here(dat_dir, "linelist.csv")
suscep_filepath <- here(dat_dir, "suscep", "suscep_1-1.csv")

experiments <- read_csv(exp_filepath, show_col_types = FALSE)
results <- read_csv(res_filepath, show_col_types = FALSE)
# sim_linelist <- read_csv(ll_filepath, show_col_types = FALSE)
suscep <- read_csv(suscep_filepath, show_col_types = FALSE)

# this will correspond to a test-negative cumulative attack rate of ~75%
beta_test_neg <- pars$beta * 0.75

# calculates the test-negative infection hazard following the same seasonal forcing
# used in the simulation model (i.e., shift = 0.5)
test_neg_infection_hazard <- function(time, beta, amplifier, shift = 0.5) {
    return(beta * (1 + (amplifier * cos(2 * pi * (time + shift)))))
}

# estimate vaccinated test-negative infections at a given time
# (taking into account the time step width)
vax_test_neg_infections <- function(time, beta, pars, shift = 0.5) {
    haz <- test_neg_infection_hazard(time, beta, pars$seasonality_amplifier, shift)
    return(pars$pop_size * pars$vax_coverage * haz * pars$dt)
}

# estimate unvaccinated test-negative infections at a given time
# (taking into account the time step width)
unvax_test_neg_infections <- function(time, beta, pars, shift = 0.5) {
    haz <- test_neg_infection_hazard(time, beta, pars$seasonality_amplifier, shift)
    return(pars$pop_size * (1 - pars$vax_coverage) * haz * pars$dt)
}

print("CLEANING INFECTION DATA AND ESTIMATING TEST-NEG INFECTIONS")

per_10k <- 1e3 / pars$pop_size

# gather test-positive and test-negative infections and adjust per 10k people
inf_dt <- results %>%
    filter(t > 0.005) %>%
    mutate(
        step = t / pars$dt,
        year = ceiling(t),
        time = t - (year - 1)
    ) %>%
    group_by(exp, t) %>%
    mutate(
        v_tp_inf = mean(vax_inf) * per_10k,
        u_tp_inf = mean(unvax_inf) * per_10k,
        v_tn_inf = vax_test_neg_infections(t, beta_test_neg, pars) * per_10k,
        u_tn_inf = unvax_test_neg_infections(t, beta_test_neg, pars) * per_10k
    ) %>%
    ungroup() %>%
    select(exp, t, year, ends_with("_tp_inf"), ends_with("_tn_inf")) %>%
    pivot_longer(!c(exp, t, year), values_to = "per_10k_inc") %>%
    mutate(count = as.integer(per_10k_inc / per_10k))

# generate linelists from the counts of cases and controls
vax_tp_linelist <- inf_dt %>%
    filter(name == "v_tp_inf") %>%
    mutate(
        count = as.integer(per_10k_inc / per_10k),
        vax = 1,
        inf = 1
    ) %>%
    select(-c(per_10k_inc, exp, name)) %>%
    uncount(count)

vax_tn_linelist <- inf_dt %>%
    filter(name == "v_tn_inf") %>%
    mutate(
        count = as.integer(per_10k_inc / per_10k),
        vax = 1,
        inf = 0
    ) %>%
    select(-c(per_10k_inc, exp, name)) %>%
    uncount(count)

unvax_tp_linelist <- inf_dt %>%
    filter(name == "u_tp_inf") %>%
    mutate(
        count = as.integer(per_10k_inc / per_10k),
        vax = 0,
        inf = 1
    ) %>%
    select(-c(per_10k_inc, exp, name)) %>%
    uncount(count)

unvax_tn_linelist <- inf_dt %>%
    filter(name == "u_tn_inf") %>%
    mutate(
        count = as.integer(per_10k_inc / per_10k),
        vax = 0,
        inf = 0
    ) %>%
    select(-c(per_10k_inc, exp, name)) %>%
    uncount(count)

# combine the linelists together
# a random sample of the total linelists for each year is used to minimize memory usage
# and generate sample sizes small enough for conditional logistic regression to succeeed
sample_prop <- 0.0015

linelist <- bind_rows(vax_tp_linelist, vax_tn_linelist,
                      unvax_tp_linelist, unvax_tn_linelist) %>%
    mutate(block = as.integer((t %/% 0.04) - ((year - 1) * 25))) %>%
    group_by(year) %>%
    nest() %>%
    mutate(sample_data = map(data, function(x) slice_sample(x, prop = sample_prop))) %>%
    ungroup()

# remove to reduce memory usage
rm(vax_tp_linelist, vax_tn_linelist, unvax_tp_linelist, unvax_tn_linelist)

print("CALCULATING VE ESTIMATES")

# estimate annual VE estimates and calculate average absolute difference between
# regression-based VE and VE from cumulative attack rates
sample_tot_pop_size <- pars$pop_size * sample_prop

final_ve_dt <- linelist %>%
    rowwise_mutate(ve_cohort = estimate_ve_cohort_sim(sample_data, year, sample_tot_pop_size, pars$vax_coverage)) %>%
    mutate(
        ve_clreg = map_vec(sample_data, function(x) estimate_ve_cond_logreg(x)),
        ve_ulreg = map_vec(sample_data, function(x) estimate_ve_uncond_logreg(x)),
        ve_cAR = map_vec(sample_data, function(x) estimate_ve_cumul_attack_rates_sim(x, sample_tot_pop_size, pars$vax_coverage))
    ) %>%
    select(year, starts_with("ve_")) %>%
    pivot_longer(!c(year, ve_cAR))

ve_comp <- final_ve_dt %>%
    mutate(ve_diff = value - ve_cAR) %>%
    group_by(name) %>%
    arrange(ve_diff) %>%
    summarize(
        mean_diff = mean(ve_diff, na.rm = TRUE),
        lower_diff = quantile(ve_diff, probs = c(0.025), na.rm = TRUE),
        upper_diff = quantile(ve_diff, probs = c(0.975), na.rm = TRUE)
    )

print("Difference in VE rel. to VE^cumulative (mean, 95% interval)")
print(ve_comp)

print("CLEANING SUSCEPTIBILITY DATA")

# calculate mean pre-vaccination susceptibilities each year
# a random sample of the total susceptibility data (5 million people each year)
# is used to draw distributions
suscep_dt <- suscep %>%
    mutate(
        year = ceiling(t) + 1,
        vax_status = ifelse(vax_status == 1, "Vaccinated", "Unvaccinated")
    ) %>%
    group_by(year, vax_status) %>%
    mutate(avg_suscep = mean(susceptibility)) %>%
    slice_sample(prop = 0.01) %>%
    ungroup()

print("CLEANING VE DATA")

# data used to draw background figure stipes for each simulated year
rect_df <- inf_dt %>%
    distinct(year) %>%
    mutate(
        xmin = as.numeric(factor(year)) - 0.5,
        xmax = as.numeric(factor(year)) + 0.5,
        fill = factor(as.numeric(factor(year)) %% 2)
    )

# helper function to draw and label true vaccine protection line
true_vax_protect_line <- function(xmin, xmax, y) {
    line <- annotate(
        "segment",
        x = xmin,
        xend = xmax,
        y = y,
        linetype = "22",
        color = "gray25",
        linewidth = 1
    )

    label <- annotate(
      "text",
      x = (xmin + xmax) / 2,
      y = y + 5,
      label = "True vaccine protection"
    )

    return(list(line, label))
}

# code to help draw split violin plots used to show pre-vaccination susceptibility distributions
# from https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- dplyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
        1))
        quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
        aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
        aesthetics$alpha <- rep(1, nrow(quantiles))
        both <- cbind(quantiles, aesthetics)
        quantile_grob <- GeomPath$draw_panel(both, ...)
        ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    }
    else {
        ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
# end of code from stack overflow

print("PLOTTING INFECTION DATA")

tmax <- as.numeric(pars$tmax)
true_vax_protection <- as.numeric(pars$vax_efficacy) * 100

inf_plt <- ggplot(inf_dt) +
    geom_rect(
        data = rect_df,
        aes(xmin = xmin - 0.5, xmax = xmax - 0.5, ymin = -Inf, ymax = Inf, fill = fill),
        alpha = 0.3,
        inherit.aes = FALSE
    ) +
    aes(
        x = t,
        y = per_10k_inc,
        color = name,
        linetype = name
    ) + 
    geom_line(linewidth = 1) +
    coord_cartesian(xlim = c(0, tmax + 0.1), ylim = c(0, 4), expand = FALSE) +
    scale_x_continuous(breaks = seq(0, tmax, 1), labels = seq(0, tmax, 1)) +
    scale_fill_manual(values = c("gray90", "gray60"), guide = "none") +
    scale_color_manual(
        name = NULL,
        breaks = c("v_tp_inf", "u_tp_inf", "v_tn_inf", "u_tn_inf"),
        values = c(
            "v_tp_inf" = "dodgerblue",
            "u_tp_inf" = "darkorange",
            "v_tn_inf" = "dodgerblue",
            "u_tn_inf" = "darkorange"
        ),
        labels = c("Vaccinated test-positive", "Unvaccinated test-positive",
                   "Vaccinated test-negative", "Unvaccinated test-negative")
    ) +
    scale_linetype_manual(
        name = NULL,
        breaks = c("v_tp_inf", "u_tp_inf", "v_tn_inf", "u_tn_inf"),
        values = c(
            "v_tp_inf" = "solid",
            "u_tp_inf" = "solid",
            "v_tn_inf" = "22",
            "u_tn_inf" = "44"
        ),
        labels = c("Vaccinated test-positive", "Unvaccinated test-positive",
                   "Vaccinated test-negative", "Unvaccinated test-negative")
    ) +
    theme_cowplot() +
    background_grid(major = "y") +
    theme(
        legend.position = "top",
        legend.position.inside = c(0.1, 0.9),
        legend.key.width = unit(25, "pt"),
        legend.spacing.y = unit(0.5, "pt"),
        legend.background = element_rect(fill = "#ffffffaa")
    ) +
    guides(
        color = guide_legend(ncol = 2, direction = "horizontal")
    ) +
    labs(
        x = "Time (years)",
        y = expression(atop("Infections", "(per "*10^3*" people)"))
    )

print("PLOTTING VE DATA")

ve_dt <- final_ve_dt %>%
    pivot_wider(names_from = name, values_from = value) %>%
    pivot_longer(!year)

final_ve_plt <- ggplot(ve_dt) +
geom_rect(
        data = rect_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
        alpha = 0.3,
        inherit.aes = FALSE
    ) +
    aes(
        x = year,
        y = value,
        color = name
    ) +
    true_vax_protect_line(0, tmax + 0.5, true_vax_protection) +
    geom_hline(yintercept = 0, color = "gray25", linewidth = 0.5) +
    geom_point(
        size = 2,
        position = position_dodge(width = 0.5)
    ) +
    coord_cartesian(xlim = c(0.5, tmax + 0.5), ylim = c(-5, true_vax_protection + 10), expand = FALSE) +
    scale_x_continuous(breaks = seq(0, tmax, 1), labels = seq(0, tmax, 1)) +
    scale_y_continuous(breaks = seq(-20, true_vax_protection, 10), labels = seq(-20, true_vax_protection, 10)) +
    scale_fill_manual(values = c("gray90", "gray60"), guide = "none") +
    scale_color_manual(
        name = NULL,
        breaks = c("ve_cAR", "ve_ulreg", "ve_clreg", "ve_cohort"),
        values = c("black", "darkorchid", "forestgreen", "maroon"),
        labels = c(
            expression(widehat(VE)^{cumulative}),
            expression(widehat(VE)^{unconditional}),
            expression(widehat(VE)^{conditional}),
            expression(widehat(VE)^{cohort})
        )
    ) +
    theme_cowplot() +
    background_grid(major = "y") +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.5, 0.2),
        legend.background = element_rect(fill = "#ffffff")
    ) +
    guides(color = guide_legend(nrow = 2)) +
    labs(
        x = "Year",
        y = expression("Final" ~ widehat(VE) ~ "estimate (%)")
    )

print("PLOTTING SUSCEPTIBILITY DATA")

suscep_plt <- ggplot(suscep_dt) +
    geom_rect(
        data = rect_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
        alpha = 0.3,
        inherit.aes = FALSE
    ) +
    aes(
        x = factor(t),
        y = susceptibility,
        fill = factor(vax_status)
    ) +
    geom_split_violin(scale = "width") +
    geom_point(
        aes(y = avg_suscep, color = "avg"),
        position = position_dodge(width = 1),
        size = 2,
        shape = 18
    ) +
    coord_cartesian(xlim = c(0.5, tmax + 0.6), ylim = c(0, 3.0), expand = FALSE) +
    scale_x_discrete(labels = NULL) +
    scale_fill_manual(
        name = NULL,
        breaks = c("Vaccinated", "Unvaccinated"),
        values = c(
            "0" = "gray90",
            "1" = "gray60",
            "Vaccinated" = "dodgerblue",
            "Unvaccinated" = "darkorange"),
            guide = "none"
    ) +
    scale_color_manual(
        name = NULL,
        breaks = c("avg"),
        labels = c("Mean pre-vaccination\nsusceptibility"),
        values = c("black")
    ) +
    theme_cowplot() +
        background_grid(major = "y") +
        theme(
        legend.position = "inside",
        legend.position.inside = c(0.6, 0.95),
        legend.background = element_rect(fill = "#ffffff")
    ) +
    labs(
        x = "Density",
        y = "Pre-vaccination\nsusceptibility"
    )

print("SAVING FINAL PLOTS")

plt <- plot_grid(
    inf_plt, final_ve_plt, suscep_plt,
    nrow = 3,
    labels = "AUTO",
    align = "v",
    axis = "lr"
)

fig_path <- here("plots", "supplemental_figs")
dir.create(fig_path)

ggsave(
  here(fig_path, "multiyear_simulation_ve_comparison.png"),
  plt,
  width = 6,
  height = 9,
  units = "in",
  bg = "white"
)

ggsave(
  here(fig_path, "multiyear_simulation_ve_comparison.pdf"),
  plt,
  width = 6,
  height = 9,
  units = "in",
  bg = "white"
)
