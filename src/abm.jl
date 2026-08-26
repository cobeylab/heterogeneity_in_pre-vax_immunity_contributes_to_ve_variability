# Define possible disease states for an individual
@enum DiseaseState begin
    Susceptible
    Infected
end

# Define possible vaccination status states for an individual
@enum VaccinationStatus begin
    Unvaccinated
    Vaccinated
end

# Define model parameters and their data types (must match TOML parameter file)
@kwdef struct Parameters{F<:AbstractFloat, I<:Integer, S<:AbstractString, B<:Bool}
    exp_name::S                   # experiment name
    exp_idx::I                    # experiment ID number

    save_linelist::B              # should the infection linelist be saved and reported
    save_current_suscep::B        # should pre-vax susceptibilities be saved and reported

    tmp_dir::S                    # path to directory for temporary data
    sim_dir::S                    # path to directory for simulation data

    nrep::I                       # number of simulation replicates

    tmax::F                       # maximum simulation time
    t_boost_transmission::F       # time when increased transmission occurs
    dt::F                         # simulation time step size

    pop_size::I                   # synthetic population size

    beta::F                       # exogenous infection hazard
    beta_boosted::F               # infection hazard during year of increased transmission

    seasonality_amplifier::F      # amplitude modifier for seasonal forcing
    seasonality_shift::F          # phase shift for seasonal forcing

    gamma::F                      # 1 / average infection duration

    vax_mean_pre_vax_suscep::F    # initial vaccinaed pre-vax suscep distribution mean
    unvax_mean_pre_vax_suscep::F  # initial unvaccinaed pre-vax suscep distribution mean

    vax_pre_vax_suscep_shape::F   # initial vaccinaed pre-vax suscep distribution shape
    unvax_pre_vax_suscep_shape::F # initial unvaccinaed pre-vax suscep distribution shape

    inf_imm_halflife::F           # half life for exponential waning of infection immunity

    allow_reinf::B                # can individuals be re-infected

    vax_coverage::F               # approximate vaccination coverage
    pr_vax_to_vax::F              # probability that a vaccinated person re-vaccinates
    vaccination_checkpoint::F     # time point when re-vaccination occurs
    t_vax_start::F                # time when vaccination first occurs

    true_vax_protection::F               # true vaccine protection
    vax_imm_halflife::F           # half life for exponential waning of vaccine immunity
end

function warn_if_invalid(valid, message, value)
    if !valid
        println(message, " (current value: ", value, ")")
    end
end

function validate_parameters(par::Parameters)
    checks::Vector{Bool} = []

    # time step must end in 1 or 5 so that it creates an integral number of steps
    dt_valid = string(par.dt)[end] == '1' || string(par.dt)[end] == '5'
    push!(checks, dt_valid)
    warn_if_invalid(
        dt_valid,
        "dt_valid must end in 1 or 5 to create an integral number of steps",
        par.dt
    )

    # population must have at least one person
    pop_size_valid = par.pop_size > 0
    push!(checks, pop_size_valid)
    warn_if_invalid(pop_size_valid, "pop_size must be larger than zero", par.pop_size)

    # infection hazards must be zero or postive
    beta_valid = par.beta >= 0
    push!(checks, beta_valid)
    warn_if_invalid(
        beta_valid,
        "beta must be larger then or equal to zero",
        par.beta
    )

    beta_boosted_valid = par.beta_boosted >= 0
    push!(checks, beta_boosted_valid)
    warn_if_invalid(
        beta_boosted_valid,
        "beta_boosted must be larger then or equal to zero",
        par.beta_boosted
    )

    # seasonal forcing amplitude must be between 0 and 1
    seasonality_amplitude_valid = par.seasonality_amplifier >= 0.0 && par.seasonality_amplifier <= 1.0
    push!(checks, seasonality_amplitude_valid)
    warn_if_invalid(
        seasonality_amplitude_valid,
        "seasonality_amplitude must be between 0.0 and 1.0",
        par.seasonality_amplifier
    )

    # recovery rate must be greater than zero
    gamma_valid = par.gamma > 0
    push!(checks, gamma_valid)
    warn_if_invalid(gamma_valid, "gamma must be greater than zero", par.gamma)

    # pre-vaccination susceptibility means must be greater than zero
    vax_mean_pre_vax_suscep_valid = par.vax_mean_pre_vax_suscep > 0
    push!(checks, vax_mean_pre_vax_suscep_valid)
    warn_if_invalid(
        vax_mean_pre_vax_suscep_valid,
        "vax_mean_pre_vax_suscep must be greater than zero",
        par.vax_mean_pre_vax_suscep
    )

    unvax_mean_pre_vax_suscep_valid = par.unvax_mean_pre_vax_suscep > 0
    push!(checks, unvax_mean_pre_vax_suscep_valid)
    warn_if_invalid(
        unvax_mean_pre_vax_suscep_valid,
        "unvax_mean_pre_vax_suscep must be greater than zero",
        par.unvax_mean_pre_vax_suscep
    )

    # pre-vaccination susceptibility distribution shape parameters must be greater than zero
    vax_pre_vax_suscep_shape_valid = par.vax_pre_vax_suscep_shape >= 0
    push!(checks, vax_pre_vax_suscep_shape_valid)
    warn_if_invalid(
        vax_pre_vax_suscep_shape_valid,
        "vax_pre_vax_suscep_shape must be greater than or equal to zero",
        par.vax_pre_vax_suscep_shape
    )

    unvax_pre_vax_suscep_shape_valid = par.unvax_pre_vax_suscep_shape >= 0
    push!(checks, unvax_pre_vax_suscep_shape_valid)
    warn_if_invalid(
        unvax_pre_vax_suscep_shape_valid,
        "unvax_pre_vax_suscep_shape must be greater than or equal to zero",
        par.unvax_pre_vax_suscep_shape
    )

    # infection-derived immunity waning rate must be greater than zero
    inf_imm_halflife_valid = par.inf_imm_halflife >= 0
    push!(checks, inf_imm_halflife_valid)
    warn_if_invalid(
        inf_imm_halflife_valid,
        "inf_imm_halflife must be greater than or equal to zero",
        par.inf_imm_halflife
    )

    # vaccination coverage must be between 0.0 and 1.0
    vax_coverage_valid = par.vax_coverage >= 0.0 && par.vax_coverage <= 1.0
    push!(checks, vax_coverage_valid)
    warn_if_invalid(vax_coverage_valid, "vax_coverage must be between 0.0 and 1.0", par.vax_coverage)

    # vaccination status auto-correlation must be between 0.0 and 1.0
    pr_vax_to_vax_valid = par.pr_vax_to_vax >= 0.0 && par.pr_vax_to_vax <= 1.0
    push!(checks, pr_vax_to_vax_valid)
    warn_if_invalid(pr_vax_to_vax_valid, "pr_vax_to_vax must be between 0.0 and 1.0", par.pr_vax_to_vax)

    # true vaccine protection must be between 0.0 and 1.0
    true_vax_protection_valid = par.true_vax_protection >= 0.0 && par.true_vax_protection <= 1.0
    push!(checks, true_vax_protection_valid)
    warn_if_invalid(
        true_vax_protection_valid,
        "true_vax_protection must be between 0.0 and 1.0",
        par.true_vax_protection
    )

    # vaccine-derived immunity waning rate must be greater than zero
    vax_imm_halflife_valid = par.vax_imm_halflife >= 0
    push!(checks, vax_imm_halflife_valid)
    warn_if_invalid(
        vax_imm_halflife_valid,
        "vax_imm_halflife must be greater than or equal to zero",
        par.vax_imm_halflife
    )

    return sum(checks) == length(checks)
end

# Helper function to create parameter object from TOML file
Parameters(d::Dict{String, Any}) = Parameters(; (Symbol.(keys(d)) .=> values(d))... )

# Object to store infection data
@kwdef mutable struct Infection
    t_infection::Float64 = -1.0     # time of infection
    t_recovery::Float64 = -1.0      # time of recovery

    suscep_at_inf::Float64 = -1.0   # infectee's susceptibility when infected

    vax_status::VaccinationStatus   # infectee's vaccination status when infected
    vax_eff_at_inf::Float64 = -1.0  # infectee's vaccine protection when infected
end

# Object to store vaccination data
@kwdef struct Vaccination
    t_vaccination::Float64 = -1.0   # time of vaccination
    true_vax_protection::Float64 = -1.0    # starting true vaccine protection
end

# Object to store data for individuals
@kwdef mutable struct Person
    id::Int64 = -1                               # person ID
    state::DiseaseState = Susceptible            # current disease state
    nextstate::DiseaseState = Susceptible        # next disease state
    vax_status::VaccinationStatus = Unvaccinated # current vaccination status

    tnow::Float64 = 0.0                          # current simulation time
    tnext::Float64 = 0.0                         # simulation time for next event

    base_suscep::Float64 = 1.0                   # initial sampled susceptibility
    current_suscep::Float64 = 1.0                # current susceptibility

    foi::Float64 = 0.0                           # total force of infection this person faces
    gamma::Float64 = 0.0                         # recovery rate

    inf_hist::Vector{Infection} = []             # stores this person's infections
    vax_hist::Vector{Vaccination} = []           # stores this person's vaccinations
end

# Check if a person has ever been vaccinated during the simulation
function has_been_vaccinated(p::Person)
    return length(p.vax_hist) > 0
end

# Get a person's vaccination status
function vax_status(p::Person)
    return p.vax_status
end

# Check if a person is currently vaccinated
function is_vaccinated(p::Person)
    return vax_status(p) == Vaccinated
end

# Get a vaccinated person's current vaccine protection
function current_true_vax_protection(par::Parameters, p::Person)
    @assert is_vaccinated(p) == true
    if par.vax_imm_halflife == zero(par.vax_imm_halflife)
        return p.vax_hist[end].true_vax_protection
    else
        t_delta = p.tnow - p.vax_hist[end].t_vaccination
        waning_mult = exp(-t_delta * (log(2) / par.vax_imm_halflife))
        return p.vax_hist[end].true_vax_protection * waning_mult
    end
end

# Check if a person has ever been infected during the simulation
function has_been_infected(p::Person)
    return length(p.inf_hist) > 0
end

# Get a previously infected person's last recovery time
function last_recovery_time(p::Person)
    @assert has_been_infected(p) == true
    return p.inf_hist[end].t_recovery
end

# Get a previously infected person's last infection time
function last_infection_time(p::Person)
    @assert has_been_infected(p) == true
    return p.inf_hist[end].t_infection
end

# Person maintains a person's susceptible state
function susceptible!(p::Person, time)
    p.tnext = time
    p.nextstate = Susceptible::DiseaseState
end

# Person transitions to infected state (i.e., is infected)
function infect!(par::Parameters, p::Person, time)
    inf = Infection(
        vax_status = vax_status(p),
        vax_eff_at_inf = is_vaccinated(p) ? current_true_vax_protection(par, p) : -1.0,
        t_infection = time,
        suscep_at_inf = p.current_suscep
    )
    push!(p.inf_hist, inf)

    p.tnext = time
    p.nextstate = Infected::DiseaseState
    p.current_suscep = zero(p.current_suscep)
end

# Person transitions from infected to susceptible state
function recover!(p::Person, time)
    inf = p.inf_hist[end]
    inf.t_recovery = time

    p.tnext = time
    p.nextstate = Susceptible::DiseaseState
end

# Sample gamma distributed susceptibility value
function sample_base_suscep(par::Parameters, vax_status)
    mean = vax_status == Vaccinated ? par.vax_mean_pre_vax_suscep : par.unvax_mean_pre_vax_suscep
    shape = vax_status == Vaccinated ? par.vax_pre_vax_suscep_shape : par.unvax_pre_vax_suscep_shape

    if shape == zero(shape)
        # if the shape parameter is zero, assume point density and return the mean
        return mean
    else
        # otherwise, sample susceptibility from gamma distribution
        scale = mean / shape
        distr = Gamma(shape, scale)
        suscep = rand(distr)
        return suscep
    end
end

# Initialize synthetic population
function init_pop(par::Parameters)
    people = Vector{Person}(undef, par.pop_size)
    for i in eachindex(people)
        # assign each person an initial vaccination status (if vaccination starts at t=0)
        coverage = par.t_vax_start == zero(par.t_vax_start) ? par.vax_coverage : 0.0
        vs = rand() < coverage ? Vaccinated : Unvaccinated

        # sample pre-vaccination susceptibility based on vaccination status
        bs = sample_base_suscep(par, vs)

        people[i] = Person(
            id = i,
            base_suscep = bs,
            foi = par.beta,
            gamma = par.gamma,
            state = Susceptible::DiseaseState
        )

        if vs == Vaccinated
            people[i].vax_status = Vaccinated::VaccinationStatus
            vax = Vaccination(t_vaccination = 0.0, true_vax_protection = par.true_vax_protection)
            push!(people[i].vax_hist, vax)
        end

        susceptible!(people[i], 0.0)
    end

    return people
end

# Helper function for an individual to become vaccinated
function become_vaccinated!(par::Parameters, p::Person, time)
    p.vax_status = Vaccinated::VaccinationStatus
    vax = Vaccination(t_vaccination = time, true_vax_protection = par.true_vax_protection)
    push!(p.vax_hist, vax)
end

# Helper function for an individual to become unvaccinated
function become_unvaccinated!(p::Person)
    p.vax_status = Unvaccinated::VaccinationStatus
end

# Update vaccination status of all individuals
function revaccinate!(par::Parameters, people::Vector{Person}, time)
    for p in people
        # calculate unvax-to-vax transition probability based on the parameterized
        # vax-to-vax transition probability and vaccination coverage
        pr_unvax_to_vax = (par.vax_coverage * (1 - par.pr_vax_to_vax)) / (1 - par.vax_coverage)
        vax_prob = p.vax_status == Vaccinated ? par.pr_vax_to_vax : pr_unvax_to_vax

        if rand() < vax_prob
            become_vaccinated!(par, p, time)
        else
            become_unvaccinated!(p)
        end
    end
end

# Helper function to return current susceptibility accounting for waning infection immunity
function current_suscep(par::Parameters, p::Person, time)
    if par.inf_imm_halflife == zero(par.inf_imm_halflife)
        # if the infection immunity half life is zero, there is no waning
        return p.base_suscep
    else
        if has_been_infected(p)
            t_recov = last_recovery_time(p)
            if t_recov == -1.0
                # this person is mid-infection, and has susceptibility = 0
                return zero(p.base_suscep)
            else
                # this person is not currently infected but has been previously infected
                t_delta = time - last_recovery_time(p)
                # calculate waned infection immunity
                waning_mult = 1 - exp(-t_delta * (log(2) / par.inf_imm_halflife))
                return p.base_suscep * waning_mult
            end
        else
            # otherwise return initial susceptibility
            return p.base_suscep
        end
    end
end

function sample_exponential(scale)
    distr = Exponential(scale)
    return rand(distr)
end

# Determine if currently susceptible individual becomes infected during this time step
function update_susceptible!(par::Parameters, p::Person, step_tmax::Real)
    # update individual's time, state, and susceptibility
    p.tnow = p.tnext
    p.state = p.nextstate
    p.current_suscep = current_suscep(par, p, p.tnow)

    # calculate force of infection accoutning for susceptibility and any vaccine protection
    foi = p.current_suscep * p.foi
    foi *= is_vaccinated(p) ? 1 - current_true_vax_protection(par, p) : 1.0

    # sample time to next infection based on current force of infection
    tau = p.tnow + sample_exponential(1 / foi)

    if tau < step_tmax
        # if the sampled infection time is before the end of this time step, the
        # individual becomes infected
        infect!(par, p, tau)
    else
        # otherwise, the individual remains susceptible
        susceptible!(p, step_tmax + eps())
    end
end

# Determine when an individual transitions from infected to recovered
function update_infected!(par::Parameters, p::Person)
    # update individual's time, state
    p.tnow = p.tnext
    p.state = p.nextstate

    if par.allow_reinf
        # sample when the person recovers (transitions from infected to susceptible)
        t_rec = p.tnow + sample_exponential(1 / p.gamma)
        recover!(p, t_rec)
    else
        p.tnext = Inf
        p.nextstate = Infected::DiseaseState
    end
end

# Simulate potential events for each person based on their current state until the end
# of the time step is reached
function sim_person!(par::Parameters, p::Person, step_tmax::Real)
    while p.tnext < step_tmax
        if p.nextstate == Susceptible::DiseaseState
            update_susceptible!(par, p, step_tmax)
        else
            update_infected!(par, p)
        end
    end
end

# Calculate seasonally forced infection hazard
function calculate_beta(par::Parameters, time)
    # decide if this is the time period when the increased hazard is used
    use_beta_boosted = (time > par.t_boost_transmission) && (time < (par.t_boost_transmission + 1.0))
    beta = use_beta_boosted ? par.beta_boosted : par.beta

    # calculate seasonally forced infection hazard
    cos_term = cos(2*pi*(time + par.seasonality_shift))
    return beta * (1 + (par.seasonality_amplifier * cos_term))
end

# Calculate baseline infection hazard for each person
function update_foi(par::Parameters, people::Vector{Person}, time)
    # calculate seasonally forced infection hazard
    beta = calculate_beta(par, time)

    # assign current infection hazard to each person
    for p in people
        p.foi = beta
    end

    return beta
end

# Stores simulation data to save to disc
@kwdef mutable struct MetricReport
    # simulation time at each time step
    times::Vector{Float64} = Float64[]

    # current number of people in different states in each time step
    num_susceptible::Vector{Int64} = Int64[]
    num_vaccinated::Vector{Int64}  = Int64[]

    # number of people infected during each time step
    num_vax_infected::Vector{Int64}   = Int64[]
    num_unvax_infected::Vector{Int64} = Int64[]
end

MetricReport(n_steps) = MetricReport(
    times = zeros(Float64, n_steps),
    num_susceptible = zeros(Int64, n_steps),
    num_vax_infected = zeros(Int64, n_steps),
    num_unvax_infected = zeros(Int64, n_steps),
    num_vaccinated = zeros(Int64, n_steps)
)

function generate_report(exp_idx, rep, report::MetricReport)
    df = DataFrame(
        exp = exp_idx,
        rep = rep,
        t = report.times,
        s = report.num_susceptible,
        vax_inf = report.num_vax_infected,
        unvax_inf = report.num_unvax_infected,
        vaxd = report.num_vaccinated
    )

    return df
end

# Helper function to shrink a vector to a given length, and apply a funtion to it
function safe_shrink_and_summarize!(vec, len, fun)
    resize!(vec, len)
    return length(vec) > 0 ? fun(vec) : -1.0
end

# Helper function to determin if a person was infected in the current time step
function infected_this_step(p::Person, step, dt)
    t_min = (step - 1) * dt
    t_max = step * dt
    return has_been_infected(p) && last_infection_time(p) >= t_min && last_infection_time(p) < t_max
end

# Calculate simulation metrics to save to disc
function reporter!(report::MetricReport, par::Parameters, people::Vector{Person}, t)
    step, tnow = t
    report.times[step] = tnow

    for p in people
        if p.vax_status == Vaccinated
            report.num_vaccinated[step] += 1
            report.num_vax_infected[step] += infected_this_step(p, step, par.dt) ? 1 : 0
        else
            report.num_unvax_infected[step] += infected_this_step(p, step, par.dt) ? 1 : 0
        end

        if p.state == Susceptible
            report.num_susceptible[step] += 1
        end
    end
end

# Checks if the simulation has reached the parameterized checkpoint when re-vaccination occurs
function should_revaccinate(par::Parameters, time)
    sim_has_started = time > zero(time)
    revaccination_enabled = par.vaccination_checkpoint != 0.0
    checkpoint_reached = (time % par.vaccination_checkpoint) == 0
    vax_campaign_started = time >= par.t_vax_start

    return sim_has_started && revaccination_enabled && checkpoint_reached &&
           vax_campaign_started
end

# Saves all individuals' susceptibility to a CSV file
function save_current_suscep(par::Parameters, people::Vector{Person}, time, rep, init)
    saveto = datadir(par.exp_name, "suscep")
    fname = string("suscep_", par.exp_idx, "-", rep, ".csv")
    fullpath = joinpath(saveto, fname)

    vax_stat = zeros(Int8, par.pop_size)
    suscep = zeros(Float64, par.pop_size)

    for i in 1:par.pop_size
        p = people[i]
        vax_stat[i] = p.vax_status == Vaccinated ? 1 : 0
        suscep[i] = current_suscep(par, p, time)
    end

    df = DataFrame(
        exp = par.exp_idx,
        t = time,
        vax_status = vax_stat,
        susceptibility = suscep
    )

    mkpath(saveto)

    if init
        # if this is the first time the data is being saved, overwrite any existing file
        CSV.write(fullpath, df, append = false)
    else
        # otherwise, only append to the file
        CSV.write(fullpath, df, append = true)
    end
end

# Main simulation function
function simulate(par::Parameters, rep::Real)
    # only continue if parameters are valid
    params_are_valid = validate_parameters(par)
    if !params_are_valid
        println("ERROR: parameter validation failed")
        exit(1)
    end

    # initialize the RNG and time-keeping variables
    Random.seed!(rep)
    tnow = 0.0
    step = 0

    # initialize the synthetic population
    people = init_pop(par)

    # initialize the simulation metric report
    n_steps = Int(par.tmax / par.dt)
    report = MetricReport(n_steps)

    # store initial susceptibilities if desired
    if par.save_current_suscep
        save_current_suscep(par, people, tnow, rep, true)
    end

    # main simulation loop
    while tnow < par.tmax
        print("Simulating... ", round(Int, (step/n_steps)*100), "%")

        # if a re-vaccination checkpoint is reached, re-vaccinate and save current susceptibility values
        if should_revaccinate(par, tnow)
            revaccinate!(par, people, tnow)
            if par.save_current_suscep
                save_current_suscep(par, people, tnow, rep, false)
            end
        end

        # update individual force of infection
        update_foi(par, people, tnow)

        # step to the next time point
        step += 1
        tnow = round(par.dt * step, digits = 5)

        # simulate events that occur this time step for all people
        for p in people
            sim_person!(par, p, tnow)
        end

        # save simulation metrics
        reporter!(report, par, people, (step, tnow))
        print("\r")
    end

    println("\rSimulating... completed.")

    # generate data frame with simulation metrics
    df = generate_report(par.exp_idx, rep, report)

    return (df = df, people = people)
end

# Generate linelist data of each infection that occurred during the simulation
function gen_line_list(exp_idx::Real, rep::Integer, people::Vector{Person})
    person_ids = Integer[]
    vax_status = VaccinationStatus[]
    suscep_at_inf = AbstractFloat[]
    vax_eff_at_inf = AbstractFloat[]
    inf_times = AbstractFloat[]
    rec_times = AbstractFloat[]

    for p in people
        if has_been_infected(p)
            for i in p.inf_hist
                push!(person_ids, p.id)
                push!(vax_status, i.vax_status)
                push!(suscep_at_inf, i.suscep_at_inf)
                push!(vax_eff_at_inf, i.vax_eff_at_inf)
                push!(inf_times, i.t_infection)
                push!(rec_times, i.t_recovery)
            end
        end
    end

    return DataFrame(
        exp = exp_idx,
        sim_rep = rep,
        person_id = person_ids,
        vax_status = vax_status,
        suscep_at_inf = suscep_at_inf,
        vax_eff_at_inf = vax_eff_at_inf,
        inf_time = inf_times,
        rec_time = rec_times
    )
end

# Save linelist of all simulated infections to CSV
function save_linelist(par::Parameters, rep::Real, people)
    ll = gen_line_list(par.exp_idx, rep, people)

    saveto = par.sim_dir == "default" ? datadir(par.exp_name, "sims") : par.sim_dir
    fname = string("sim_", par.exp_idx, "-", rep, ".csv")
    fullpath = joinpath(saveto, fname)

    mkpath(saveto)
    println("saving sim linelist to ", saveto)
    CSV.write(fullpath, ll)
    println("done")

    return ll
end

# Save simulation metrics to CSV
function save_report(par::Parameters, rep::Real, df)
    saveto = par.sim_dir == "default" ? datadir(par.exp_name, "sims") : par.sim_dir
    fname = string("sim_", par.exp_idx, "-", rep, ".csv")
    fullpath = joinpath(saveto, fname)

    mkpath(saveto)
    println("saving sim metrics to ", saveto)
    CSV.write(fullpath, df)
    println("done")

    return df

end