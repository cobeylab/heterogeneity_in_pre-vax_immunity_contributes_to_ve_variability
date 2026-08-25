# Helper function to perform a mutate for each row of a tibble.
rowwise_mutate <- function(data, ...) {
    data %>%
        rowwise() %>%
        mutate(...) %>%
        ungroup()
}

# Generate combinations of parameters after generating the sequence of time points.
# p --> list of parameters
# include_early --> if TRUE, include a non-zero early time point
# early_time --> early time point to include (used in case calculations are NA exactly at zero)
generate_par_sets <- function(p, include_early = TRUE, early_time = 1e-3) {
    times <- seq(p$start_time, p$end_time, p$dt)
    p$time = if (include_early) c(early_time, times) else times
    return(crossing(!!!p))
}

# Generate TND case/control cumulative and short-term case/control counts.
# p --> list of parameters
# o --> list of options (o$instantaneous must be FALSE)
generate_counts <- function(p, o) {
    if (o$instantaneous == TRUE) {
        stop("ERROR: generate_counts must have o$instantaneous = FALSE")
    }

    # Generate time points
    times <- seq(p$start_time, p$end_time, p$dt)

    tibble(time = times) %>%
        # Generate cumulative cases and controls across time
        mutate(
            cumul_1_1 = test_positive_cases(times, p, o, vaccinated = TRUE),
            cumul_0_1 = test_positive_cases(times, p, o, vaccinated = FALSE),
            cumul_1_0 = test_negative_controls(times, p, o, vaccinated = TRUE),
            cumul_0_0 = test_negative_controls(times, p, o, vaccinated = FALSE)
        ) %>%
        # These steps calculate the a short-term case/control counts based on the
        # number of cumulative cases/controls up to each time point. The length
        # of the short-term time intervals is set by p$dt, which is the same time
        # interval that will be used when performing logstic regression.
        pivot_longer(starts_with("cumul_"), values_to = "cumul_count") %>%
        separate_wider_delim(name, delim = "_", names = c("metric", "vax", "inf")) %>%
        select(-metric) %>%
        group_by(vax, inf) %>%
        mutate(
            short_term_count = c(0, diff(cumul_count)),
            block = time %/% p$dt
        ) %>%
        ungroup()
}

# Generate a linelist according to the short-term case/control counts.
# counts --> tibble of case/control counts
generate_linelist <- function(counts) {
    counts %>%
        filter(block > 0) %>%
        mutate(n_infs = round(short_term_count)) %>%
        unite(block_vax_inf, c(block, vax, inf)) %>%
        select(-ends_with("_count")) %>%
        uncount(n_infs) %>%
        separate_wider_delim(block_vax_inf, delim = "_", names = c("block", "vax", "inf")) %>%
        mutate(across(c(block, vax, inf), as.integer))
}

# Helper function to calculate VE from the results of logistic regression.
ve_from_logreg <- function(coef) {
    odds_ratio <- as.numeric(exp(coef))
    return((1 - odds_ratio) * 100)
}

# Calculate VE using unconditional logistic regression adjusting for time.
estimate_ve_uncond_logreg <- function(data) {
    reg <- glm(inf ~ vax + block, data = data, family = binomial)
    return(ve_from_logreg(reg$coefficients["vax"]))
}

# Calculate VE using conditional logistic regression matching on time.
estimate_ve_cond_logreg <- function(data) {
    reg <- clogit(inf ~ vax + strata(block), data = data)
    return(ve_from_logreg(reg$coefficients["vax"]))
}

# Calculate VE using analytical models.
estimate_math_ve <- function(pars, opts) {
    analytical_VE(pars$time, pars, opts)
}

# Estimate starting VE from vaccine direct effects and mean pre-vaccination risk.
estimate_ve_starting <- function(x) {
    return((1 - (x$theta_0 * (x$epsilon_v / x$epsilon_u))) * 100)
}

# Estimate VE using analytical cohort-like VE model.
estimate_ve_cohort <- function(pars, opts) {
    v <- average_instantaneous_incidence_rate(pars$time, pars, opts, vaccinated = TRUE)
    u <- average_instantaneous_incidence_rate(pars$time, pars, opts, vaccinated = FALSE)
    return((1 - (v / u)) * 100)
}

# Helper function to calculate vaccinated and unvaccinated population size
calculate_pop_size <- function(tot_pop_size, vax_coverage, vax_status) {
    prob_vaccinated <- ifelse(vax_status == 1, vax_coverage, 1 - vax_coverage)
    return(tot_pop_size * prob_vaccinated)
}

# Estimate VE using cumulative attack rates from linelist data
estimate_ve_cumul_attack_rates_sim <- function(x, tot_pop_size, vax_coverage) {
    x %>%
        # only need test-positive infections
        filter(inf == 1) %>%
        # count total infections by vaccination status
        group_by(vax) %>%
        summarize(tot_infs = n()) %>%
        ungroup() %>%
        # calculate cumulative attack rates
        mutate(
            cAR = tot_infs / calculate_pop_size(tot_pop_size, vax_coverage, vax)
        ) %>%
        # calculate ve from risk ratio using cumulative attack rates
        select(-tot_infs) %>%
        pivot_wider(
            names_from = vax,
            values_from = cAR,
            names_glue = "{vax}_{.value}"
        ) %>%
        summarize((1 - (`1_cAR` / `0_cAR`)) * 100) %>%
        as.numeric()
}

# Estimate VE using average instantaneous incidence rates (i.e., cohort-based estimate)
# from linelist data
estimate_ve_cohort_sim <- function(x, year, tot_pop_size, vax_coverage) {
    x %>%
        # only need test-positive infections
        filter(inf == 1) %>%
        # calculate time-to-infection (pyar) and the remainder of the year after
        # censoring at infection (pyar_avoided)
        mutate(
            pyar = t - (year - 1),
            pyar_avoided = 1 - pyar
        ) %>%
        # count total infections and pyar avoided by vaccination status
        group_by(vax) %>%
        summarize(
            tot_infs = n(),
            tot_pyar_avoided = sum(pyar_avoided)
        ) %>%
        ungroup() %>%
        # calculate sampled population size by vaccination status
        # calculate total pyar by subtracting total pyar avoided from (1 year * pop size)
        # calculate average instantaneous incidence rates
        mutate(
            sample_pop_size = calculate_pop_size(tot_pop_size, vax_coverage, vax),
            tot_pyar = sample_pop_size - tot_pyar_avoided,
            avg_insnt_inc_rate = tot_infs / tot_pyar
        ) %>%
        # calculate ve from ratio of avg instantaneous incidence rates
        select(vax, avg_insnt_inc_rate) %>%
        pivot_wider(names_from = vax, names_glue = "{vax}_{.value}", values_from = avg_insnt_inc_rate) %>%
        summarize((1 - (`1_avg_insnt_inc_rate` / `0_avg_insnt_inc_rate`)) * 100) %>%
        as.numeric()
}