library(here)
library(tidyverse)
library(optparse)
library(RcppTOML)
library(cowplot)
library(paletteer)

set.seed(0)

# parse command line arguments
if (interactive()) {
    parser <- OptionParser()
    parser <- add_option(parser, "--config", help = "path to config TOML file (rel path from this script)")
    parsed_args <- parse_args(parser, args = c("--config=fig3.toml"))
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
script_dir <- here("scripts", "fig3")
config_path <- here(script_dir, parsed_args[["config"]])

config <- parseTOML(config_path)
pars <- config$parameters
exp_name <- pars$exp_name

# read in simulation data
print("READING DATA")

dat_dir <- here("data", exp_name)
exp_filepath <- here(dat_dir, "experiments.csv")
res_filepath <- here(dat_dir, "results.csv")
suscep_filepath <- here(dat_dir, "suscep", "suscep_1-1.csv")

experiments <- read_csv(exp_filepath, show_col_types = FALSE)
results <- read_csv(res_filepath, show_col_types = FALSE)
suscep <- read_csv(suscep_filepath, show_col_types = FALSE)

print("CLEANING INFECTION DATA")

# grab the vaccinated and unvaccinated population sizes for each simualted year
vax_coverage_by_year <- results %>%
    mutate(year = ceiling(t)) %>%
    group_by(year) %>%
    distinct(vaxd) %>%
    rename(vaxd_pop_size = vaxd) %>%
    mutate(unvaxd_pop_size = pars$pop_size - vaxd_pop_size)

# calculate mean daily infections per 10k people by vaccination status
inf_dt <- results %>%
    mutate(
        step = t / pars$dt,
        year = ceiling(t),
        time = t - (year - 1)
    ) %>%
    left_join(vax_coverage_by_year, by = "year") %>%
    group_by(exp, t) %>%
    summarize(
        v_inf = mean(vax_inf) / vaxd_pop_size * 1e3,
        u_inf = mean(unvax_inf) / unvaxd_pop_size * 1e3
    ) %>%
    ungroup() %>%
    pivot_longer(!c(exp, t))

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

# ve is calculated from time-varying cumulative attack rates each year
ve_dt <- results %>%
    mutate(
        step = t / pars$dt,
        year = ceiling(t),
        time = t - (year - 1)
    ) %>%
    group_by(exp, rep, year) %>%
    mutate(
        cumul_vax_inf = cumsum(vax_inf),
        cumul_unvax_inf = cumsum(unvax_inf)
    ) %>%
    ungroup() %>%
    select(exp, rep, year, time, cumul_unvax_inf, cumul_vax_inf) %>%
    left_join(vax_coverage_by_year, by = "year") %>%
    mutate(
        cumul_vax_inf = cumul_vax_inf / vaxd_pop_size,
        cumul_unvax_inf = cumul_unvax_inf / unvaxd_pop_size,
        estd_ve = (1 - (cumul_vax_inf / cumul_unvax_inf)) * 100
    )

# data used to draw background figure stipes for each simulated year
rect_df <- ve_dt %>%
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

inf_plt <- ggplot(inf_dt) +
    geom_rect(
        data = rect_df,
        aes(xmin = xmin - 0.5, xmax = xmax - 0.5, ymin = -Inf, ymax = Inf, fill = fill),
        alpha = 0.3,
        inherit.aes = FALSE
    ) +
    aes(
        x = t,
        y = value,
        color = name
    ) + 
    geom_line(linewidth = 1) +
    coord_cartesian(xlim = c(0, tmax + 0.1), ylim = c(0, 8), expand = FALSE) +
    scale_x_continuous(breaks = seq(0, tmax, 1), labels = seq(0, tmax, 1)) +
    scale_fill_manual(values = c("gray90", "gray60"), guide = "none") +
    scale_color_manual(
        name = NULL,
        breaks = c("v_inf", "u_inf"),
        values = c(
            "v_inf" = "dodgerblue",
            "u_inf" = "darkorange"
        ),
        labels = c("Vaccinated", "Unvaccinated")
    ) +
    theme_cowplot() +
    background_grid(major = "y") +
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.725, 0.9),
        legend.background = element_rect(fill = "#ffffffaa")
    ) +
    labs(
        x = "Time (years)",
        y = expression(atop("Infections", "(per "*10^3*" people)"))
    )

print("PLOTTING VE DATA")

true_vax_protection <- as.numeric(pars$true_vax_protection) * 100

# remove early noisy ve trajectory
ve_dt <- ve_dt %>%
  mutate(estd_ve = ifelse(time < 0.1, NA, estd_ve))

final_ve_plt <- ggplot(ve_dt) +
    geom_rect(
        data = rect_df,
        aes(xmin = xmin - 0.5, xmax = xmax - 0.5, ymin = -Inf, ymax = Inf, fill = fill),
        alpha = 0.3,
        inherit.aes = FALSE
    ) +
    aes(
        x = time + (year - 1),
        y = estd_ve
    ) +
    true_vax_protect_line(0, tmax + 0.1, true_vax_protection) +
    geom_hline(yintercept = 0, color = "gray25", linewidth = 0.5) +
    geom_line(color = "gray50", linewidth = 1) +
    geom_point(
        data = ve_dt %>% filter(time == 1.0),
        size = 2
    ) +
    coord_cartesian(xlim = c(0, tmax + 0.1), ylim = c(-5, true_vax_protection + 10), expand = FALSE) +
    scale_x_continuous(breaks = seq(0, tmax, 1), labels = seq(0, tmax, 1)) +
    scale_y_continuous(breaks = seq(-20, true_vax_protection, 10), labels = seq(-20, true_vax_protection, 10)) +
    scale_fill_manual(values = c("gray90", "gray60"), guide = "none") +
    scale_color_manual(
        name = NULL,
        breaks = c(0.5, 0.85),
        values = c("dodgerblue", "darkorange"),
        labels = c("Random vaccination", "Auto-correlated vaccination")
    ) +
    theme_cowplot() +
    background_grid(major = "y") +
    theme(legend.position = "top") +
    labs(
        x = "Time (years)",
        y = "Estimated vaccine\neffectiveness (%)"
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

fig_dir <- ifelse(exp_name == "fig3", here("plots"), here("plots", "supplemental_figs"))
dir.create(fig_dir, recursive = TRUE)

ggsave(
  here(fig_dir, paste0(exp_name, ".png")),
  plt,
  width = 6,
  height = 9,
  units = "in",
  bg = "white"
)

ggsave(
  here(fig_dir, paste0(exp_name, ".pdf")),
  plt,
  width = 6,
  height = 9,
  units = "in",
  bg = "white"
)
