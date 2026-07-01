# p := expected parameter list
# p = list(
#     start_time = start time for data generation,
#     end_time = end time for data generation,
#     dt = time intervatal used to generate sequence of time points,
#     lambda = exogenous test-positive infection hazard,
#     lambda_negative = exogenous test-negative infection hazard,
#     theta_0 = starting vaccine direct effect,
#     eta = vaccine protection half life,
#     epsilon_v = mean vaccinated baseline susceptibility,
#     epsilon_u = mean unvaccinated baseline susceptibility,
#     alpha_v = vaccinated baseline susceptibility distribution shape parameter,
#     alpha_u = unvaccinated baseline susceptibility distribution shape parameter,
#     mu_v = probability of seeking care given vaccinated,
#     mu_u = probability of seeking care given unvaccinated,
#     pi_pos = probability of symptoms given test-positive infection,
#     pi_neg = probability of symptoms given test-negative infection,
#     c = vaccination coverage,
#     N = total population size
# )

# o := expected options
# o = list(
#     waning = waning vaccine direct effects,
#     heterogeneity = include full gamma-distributed baseline susceptibility (FALSE for point densities),
#     instantaneous = calculate using instantaneous incidence rates
# )

# helper function to check if p is a valid probability value
# p --> numeric arg
is_probability <- function(p) {
    stopifnot(p >= 0, p <= 1)
    return(TRUE)
}

# calculates vaccine direct effects (i.e., true individual-level vaccine protection) at a given time
# time --> numeric arg
# p --> list of parameters (see above)
# o --> list of options (see above)
#   o$waning = TRUE --> use exponential vaccine protection waning model
#   o$waning = FALSE --> use constant vaccine protection model
vaccine_direct_effect <- function(time, p, o) {
    stopifnot(time >= 0, is_probability(p$theta_0))

    if (o$waning) {
        stopifnot(p$eta > 0)
        return(1 - ((1 - p$theta_0) * exp(-time * (log(2) / p$eta))))
    } else {
        return(rep(p$theta_0, times = length(time)))
    }
}

# calculates the total vaccine direct effects over time interval [0,t] (i.e., integral of vaccine direct effects)
# time --> numeric arg
# p --> list of parameters (see above)
# o --> list of options (see above)
#   o$waning = TRUE --> integral of exponential vaccine protection waning model
#   o$waning = FALSE --> integral of constant vaccine protection model
total_vaccine_effect <- function(time, p, o) {
    stopifnot(time >= 0, is_probability(p$theta_0))

    if (o$waning) {
        stopifnot(p$eta > 0)
        if (is.infinite(p$eta)) {
            return(p$theta_0 * time)
        } else {
            numer <- (2^(-time / p$eta)) * (2^(time / p$eta) - 1) * p$eta * (p$theta_0 - 1)
            return(time + (numer / log(2)))
        }
    } else {
        return(p$theta_0 * time)
    }
}

# calculates the susceptible fraction at a given time
# time --> numeric arg
# p --> list of parameters (see above)
# o --> list of options (see above)
#   o$hetergeneity = TRUE --> calculate using gamma-distributed base suscep model
#   o$hetergeneity = FALSE --> calculate using point-density base suscep model
# vaccinated --> arg controlling if calculating for vaccinated (TRUE) or unvaccinated (FALSE) population
susceptible_fraction <- function(time, p, o, vaccinated = FALSE) {
    stopifnot(time >= 0, p$epsilon_v > 0, p$epsilon_u > 0, p$alpha_v > 0,
              p$alpha_u > 0, p$lambda > 0)

    if (vaccinated) {
        Theta <- total_vaccine_effect(time, p, o)
        if (o$heterogeneity) {
            base <- p$alpha_v / (p$alpha_v + (p$epsilon_v * p$lambda * Theta))
            return(base ^ p$alpha_v)
        } else {
            return(exp(-p$epsilon_v * p$lambda * Theta))
        }
    } else {
        if (o$heterogeneity) {
            base <- p$alpha_u / (p$alpha_u + (p$epsilon_u * p$lambda * time))
            return(base ^ p$alpha_u)
        } else {
            return(exp(-p$epsilon_u * p$lambda * time))
        }
    }
}

# calculates cumulative attack rate at a given time
# time --> numeric arg
# p --> list of parameters (see above)
# o --> list of options (see above)
# vaccinated --> arg controlling if calculating for vaccinated (TRUE) or unvaccinated (FALSE) population
cumulative_attack_rate <- function(time, p, o, vaccinated = FALSE) {
    return(1 - susceptible_fraction(time, p, o, vaccinated))
}

# calculates instantaneous incidence rate at a given time
# time --> numeric arg
# p --> list of parameters (see above)
# o --> list of options (see above)
# vaccinated --> arg controlling if calculating for vaccinated (TRUE) or unvaccinated (FALSE) population
instantaneous_incidence_rate <- function(time, p, o, vaccinated = FALSE) {
    stopifnot(time >= 0, p$epsilon_v > 0, p$epsilon_u > 0, p$alpha_v > 0,
              p$alpha_u > 0, p$lambda > 0)

    if (vaccinated) {
        Theta <- total_vaccine_effect(time, p, o)
        theta <- vaccine_direct_effect(time, p, o)
        if (o$heterogeneity) {
            base <- p$alpha_v / (p$alpha_v + (p$epsilon_v * p$lambda * Theta))
            return(p$lambda * p$epsilon_v * theta * (base ^ (p$alpha_v + 1)))
        } else {
            return(p$lambda * p$epsilon_v * theta * exp(-p$epsilon_v * p$lambda * Theta))
        }
    } else {
        if (o$heterogeneity) {
            base <- p$alpha_u / (p$alpha_u + (p$epsilon_u * p$lambda * time))
            return(p$lambda * p$epsilon_u * (base ^ (p$alpha_u + 1)))
        } else {
            return(p$lambda * p$epsilon_u * exp(-p$epsilon_u * p$lambda * time))
        }
    }
}

# calculates expected count of test-positive cases at a given time
# time --> numeric arg
# p --> list of parameters (see above)
# o --> list of options (see above)
#   o$instantaneous = TRUE --> calculate using instantaneous incidence rates
#   o$instantaneous = FALSE --> calculate using cumulative attack rates
# vaccinated --> arg controlling if calculating for vaccinated (TRUE) or unvaccinated (FALSE) population
test_positive_cases <- function(time, p, o, vaccinated = FALSE) {
    stopifnot(p$N > 0, is_probability(p$mu_v), is_probability(p$mu_u), is_probability(p$pi_pos), 
              is_probability(p$pi_neg), is_probability(p$c), p$dt > 0)

    # to approximate the number of cases that occur at this time, multiply the instantaneous
    # incidence rate by the time interval p$dt (NOTE: this is still an approximation and
    # these approximate coutns will not exactly sum to the counts from the cumulative
    # attack rates)
    pr_infection <- if (o$instantaneous)
        instantaneous_incidence_rate(time, p, o, vaccinated) * p$dt else
        cumulative_attack_rate(time, p, o, vaccinated)
    
    if (vaccinated) {
        return(p$N * p$mu_v * p$pi_pos * p$c * pr_infection)
    } else {
        return(p$N * p$mu_u * p$pi_pos * (1 - p$c) * pr_infection)
    }
}

# calculates expected count of test-negative controls at a given time
# time --> numeric arg
# p --> list of parameters (see above)
# o --> list of options (see above)
#   o$instantaneous = TRUE --> calculate using instantaneous incidence rates
#   o$instantaneous = FALSE --> calculate using cumulative attack rates
# vaccinated --> arg controlling if calculating for vaccinated (TRUE) or unvaccinated (FALSE) population
test_negative_controls <- function(time, p, o, vaccinated = FALSE) {
    stopifnot(time >= 0, p$lambda_negative > 0, p$N > 0, is_probability(p$mu_v), is_probability(p$mu_u), 
              is_probability(p$pi_pos), is_probability(p$pi_neg), is_probability(p$c), p$dt > 0)

    # to approximate the number of cases that occur at this time, multiply the instantaneous
    # incidence rate by the time interval p$dt (NOTE: this is still an approximation and
    # these approximate coutns will not exactly sum to the counts from the cumulative
    # attack rates)
    pr_infection <- if (o$instantaneous) p$lambda_negative * p$dt else p$lambda_negative * time
    
    if (vaccinated) {
        return(p$N * p$mu_v * p$pi_neg * p$c * pr_infection)
    } else {
        return(p$N * p$mu_u * p$pi_neg * (1 - p$c) * pr_infection)
    }
}

# calculates VE from the TND odds ratio
# time --> numeric arg
# p --> list of parameters (see above)
# o --> list of options (see above)
TND_VE <- function(time, p, o) {
    X_11 <- test_positive_cases(time, p, o, vaccinated = TRUE)
    X_01 <- test_positive_cases(time, p, o, vaccinated = FALSE)
    X_12 <- test_negative_controls(time, p, o, vaccinated = TRUE)
    X_02 <- test_negative_controls(time, p, o, vaccinated = FALSE)

    odds_ratio <- (X_11 / X_01) / (X_12 / X_02)
    return((1 - odds_ratio) * 100)
}

# calculates VE from the analytical models
# time --> numeric arg
# p --> list of parameters (see above)
# o --> list of options (see above)
#   o$instantaneous = TRUE --> calculate using instantaneous incidence rates
#   o$instantaneous = FALSE --> calculate using cumulative attack rates
analytical_VE <- function(time, p, o) {
    if (o$instantaneous) {
        v <- instantaneous_incidence_rate(time, p, o, vaccinated = TRUE)
        u <- instantaneous_incidence_rate(time, p, o, vaccinated = FALSE)
        return((1 - (v / u)) * 100)
    } else {
        v <- cumulative_attack_rate(time, p, o, vaccinated = TRUE)
        u <- cumulative_attack_rate(time, p, o, vaccinated = FALSE)
        return((1 - (v / u)) * 100)
    }
}
