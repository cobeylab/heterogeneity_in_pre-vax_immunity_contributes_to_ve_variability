@enum DiseaseState begin
    Susceptible
    Exposed
    Infected
    Recovered
end

@enum VaccinationStatus begin
    Unvaccinated
    Vaccinated
end

@kwdef struct Parameters{F<:AbstractFloat, I<:Integer, S<:AbstractString, B<:Bool}
    exp_name::S
    exp_idx::I

    save_linelist::B
    save_prevax_suscep::B

    tmp_dir::S
    sim_dir::S

    nrep::I

    tmax::F
    tburn::F
    dt::F
    pop_size::I
    initial_infected::I
    initial_recovered::I

    simulate_transmission::B
    pr_seeded_inf::F

    beta_burn::F
    beta_boosted::F
    # beta1::F
    seasonality_amplifier::F
    seasonality_shift::F

    sigma::F
    gamma::F
    omega::F

    vax_base_suscep_mean::F
    unvax_base_suscep_mean::F
    vax_base_suscep_shape::F
    unvax_base_suscep_shape::F
    inf_imm_halflife::F

    allow_reinf::B

    vax_coverage::F
    pr_vax_to_vax::F
    # pr_unvax_to_vax::F
    vax_efficacy::F
    vax_imm_halflife::F
    vaccination_checkpoint::F
    t_vax_start::F
end

Parameters(d::Dict{String, Any}) = Parameters(; (Symbol.(keys(d)) .=> values(d))... )

@kwdef mutable struct Infection
    t_infection::Float64 = -1.0
    t_recovery::Float64 = -1.0

    suscep_at_inf::Float64 = -1.0

    vax_status::VaccinationStatus
    vax_eff_at_inf::Float64 = -1.0
end

@kwdef struct Vaccination
    t_vaccination::Float64 = -1.0
    vax_efficacy::Float64 = -1.0
end

@kwdef mutable struct Person
    id::Int64 = -1
    state::DiseaseState = Susceptible
    nextstate::DiseaseState = Susceptible
    vax_status::VaccinationStatus = Unvaccinated

    tnow::Float64 = 0.0
    tnext::Float64 = 0.0

    base_suscep::Float64 = 1.0
    current_suscep::Float64 = 1.0

    foi::Float64 = 0.0
    sigma::Float64 = 0.0
    gamma::Float64 = 0.0
    omega::Float64 = 0.0

    inf_hist::Vector{Infection} = []
    vax_hist::Vector{Vaccination} = []
end

function has_been_vaccinated(p::Person)
    return length(p.vax_hist) > 0
end

function is_vaccinated(p::Person)
    return p.vax_status == Vaccinated
end

function vax_status(p::Person)
    return p.vax_status
end

function current_vax_efficacy(par::Parameters, p::Person)
    if par.vax_imm_halflife == zero(par.vax_imm_halflife)
        return p.vax_hist[end].vax_efficacy
    else
        t_delta = p.tnow - p.vax_hist[end].t_vaccination
        waning_mult = exp(-t_delta * (log(2) / par.vax_imm_halflife))
        return p.vax_hist[end].vax_efficacy * waning_mult
    end
end

function has_been_infected(p::Person)
    return length(p.inf_hist) > 0
end

function last_recovery_time(p::Person)
    @assert has_been_infected(p) == true
    return p.inf_hist[end].t_recovery
end

function last_infection_time(p::Person)
    @assert has_been_infected(p) == true
    return p.inf_hist[end].t_infection
end

function susceptible!(p::Person, time)
    p.tnext = time
    p.nextstate = Susceptible::DiseaseState
end

function expose!(par::Parameters, p::Person, time)
    inf = Infection(
        vax_status = vax_status(p),
        vax_eff_at_inf = is_vaccinated(p) ? current_vax_efficacy(par, p) : -1.0,
        t_infection = time,
        suscep_at_inf = p.current_suscep
    )
    push!(p.inf_hist, inf)

    p.tnext = time
    p.nextstate = Exposed::DiseaseState
    p.current_suscep = zero(p.current_suscep)
end

function infect!(p::Person, time)
    p.tnext = time
    p.nextstate = Infected::DiseaseState
end

function recover!(p::Person, time)
    inf = p.inf_hist[end]
    inf.t_recovery = time

    p.tnext = time
    p.nextstate = Recovered::DiseaseState
end

function sample_base_suscep(par::Parameters, vax_status)
    mean = vax_status == Vaccinated ? par.vax_base_suscep_mean : par.unvax_base_suscep_mean
    shape = vax_status == Vaccinated ? par.vax_base_suscep_shape : par.unvax_base_suscep_shape

    if shape == zero(shape)
        return mean
    else
        scale = mean / shape
        distr = Gamma(shape, scale)
        suscep = rand(distr)
        return suscep
    end
end

function init_pop(par::Parameters)
    people = Vector{Person}(undef, par.pop_size)
    for i in eachindex(people)
        coverage = par.t_vax_start == zero(par.t_vax_start) ? par.vax_coverage : 0.0
        vs = rand() < coverage ? Vaccinated : Unvaccinated
        bs = sample_base_suscep(par, vs)

        people[i] = Person(
            id = i,
            base_suscep = bs,
            foi = par.beta_burn,
            sigma = par.sigma,
            gamma = par.gamma,
            omega = par.omega
        )

        if vs == Vaccinated
            people[i].vax_status = Vaccinated::VaccinationStatus
            vax = Vaccination(t_vaccination = 0.0, vax_efficacy = par.vax_efficacy)
            push!(people[i].vax_hist, vax)
        end

        if i <= par.initial_infected
            people[i].state = Exposed::DiseaseState
            expose!(par, people[i], 0.0)
        elseif (i > par.initial_infected) && (i <= (par.initial_infected + par.initial_recovered))
            people[i].state = Recovered::DiseaseState
            expose!(par, people[i], 0.0)
            recover!(people[i], 0.0)
        else
            people[i].state = Susceptible::DiseaseState
            susceptible!(people[i], 0.0)
        end
    end

    return people
end

function become_vaccinated!(par::Parameters, p::Person, time)
    p.vax_status = Vaccinated::VaccinationStatus
    vax = Vaccination(t_vaccination = time, vax_efficacy = par.vax_efficacy)
    push!(p.vax_hist, vax)
end

function become_unvaccinated!(par::Parameters, p::Person)
    p.vax_status = Unvaccinated::VaccinationStatus
end

function revaccinate!(par::Parameters, people::Vector{Person}, time)
    for p in people
        pr_unvax_to_vax = (par.vax_coverage * (1 - par.pr_vax_to_vax)) / (1 - par.vax_coverage)
        vax_prob = p.vax_status == Vaccinated ? par.pr_vax_to_vax : pr_unvax_to_vax

        init_vax = (time >= par.t_vax_start) && (time < (par.t_vax_start + par.vaccination_checkpoint))
        vax_prob = init_vax ? 0.5 : vax_prob

        if rand() < vax_prob
            become_vaccinated!(par, p, time)
        else
            become_unvaccinated!(par, p)
        end
    end
end

function current_suscep(par::Parameters, p::Person, time)
    if par.inf_imm_halflife == zero(par.inf_imm_halflife)
        return p.base_suscep
    else
        if has_been_infected(p)
            t_delta = time - last_recovery_time(p)
            waning_mult = 1 - exp(-t_delta * (log(2) / par.inf_imm_halflife))
            return p.base_suscep * waning_mult
        else
            return p.base_suscep
        end
    end
end

function sample_exponential(scale)
    distr = Exponential(scale)
    return rand(distr)
end

function update_susceptible!(par::Parameters, p::Person, step_tmax::Real)
    p.tnow = p.tnext
    p.state = p.nextstate
    p.current_suscep = current_suscep(par, p, p.tnow)

    foi = p.current_suscep * p.foi
    foi *= is_vaccinated(p) ? 1 - current_vax_efficacy(par, p) : 1.0
    tau = p.tnow + sample_exponential(1 / foi)

    if tau < step_tmax
        expose!(par, p, tau)
    else
        susceptible!(p, step_tmax + eps())
    end
end

function update_exposed!(p::Person)
    p.tnow = p.tnext
    p.state = p.nextstate

    t_inf = p.tnow + sample_exponential(1 / p.sigma)
    infect!(p, t_inf)
end

function update_infected!(p::Person)
    p.tnow = p.tnext
    p.state = p.nextstate

    t_rec = p.tnow + sample_exponential(1 / p.gamma)
    recover!(p, t_rec)
end

function update_recovered!(par::Parameters, p::Person)
    p.tnow = p.tnext
    p.state = p.nextstate

    if par.allow_reinf
        t_wane = p.tnow + sample_exponential(1 / p.omega)
        susceptible!(p, t_wane)
    else
        p.tnext = Inf
        p.nextstate = Recovered::DiseaseState
    end
end

function sim_person!(par::Parameters, p::Person, step_tmax::Real)
    while p.tnext < step_tmax
        if p.nextstate == Susceptible::DiseaseState
            update_susceptible!(par, p, step_tmax)
        elseif p.nextstate == Exposed::DiseaseState
            update_exposed!(p)
        elseif p.nextstate == Infected::DiseaseState
            update_infected!(p)
        else
            update_recovered!(par, p)
        end
    end
end

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

function calculate_beta(par::Parameters, time)
    first_year_after_burnin = (time > par.tburn) && (time < (par.tburn + 1.0))
    beta = first_year_after_burnin ? par.beta_boosted : par.beta_burn
    cos_term = cos(2*pi*(time + par.seasonality_shift))
    return beta * (1 + (par.seasonality_amplifier * cos_term))
end

function update_foi(par::Parameters, people::Vector{Person}, time)
    beta = calculate_beta(par, time)

    if par.simulate_transmission
        num_infected = 0
        for p in people
            if p.state == Infected::DiseaseState
                num_infected += 1
            end
        end

        beta *= (num_infected / par.pop_size)
    end

    for p in people
        p.foi = beta
    end

    return beta
end

function seed_infections!(par::Parameters, people::Vector{Person}, time)
    for p in people
        if (p.state == Susceptible::DiseaseState) && (rand() < par.pr_seeded_inf)
            p.state = Exposed::DiseaseState
            expose!(par, p, time)
        end
    end
end

@kwdef mutable struct MetricReport
    times::Vector{Float64}          = Float64[]
    num_susceptible::Vector{Int64}  = Int64[]
    num_exposed::Vector{Int64}      = Int64[]
    num_infected::Vector{Int64}     = Int64[]
    num_recovered::Vector{Int64}    = Int64[]

    num_vax_infected::Vector{Int64} = Int64[]
    num_unvax_infected::Vector{Int64} = Int64[]

    avg_susceptibility::Vector{Float64} = Float64[]
    avg_unvax_susceptibility::Vector{Float64} = Float64[]
    avg_vax_susceptibility::Vector{Float64}   = Float64[]

    var_susceptibility::Vector{Float64} = Float64[]
    var_unvax_susceptibility::Vector{Float64} = Float64[]
    var_vax_susceptibility::Vector{Float64}   = Float64[]

    min_quant_unvax_suscep::Vector{Float64} = Float64[]
    min_quant_vax_suscep::Vector{Float64} = Float64[]

    low_quant_unvax_suscep::Vector{Float64} = Float64[]
    low_quant_vax_suscep::Vector{Float64} = Float64[]

    med_quant_unvax_suscep::Vector{Float64} = Float64[]
    med_quant_vax_suscep::Vector{Float64} = Float64[]

    upp_quant_unvax_suscep::Vector{Float64} = Float64[]
    upp_quant_vax_suscep::Vector{Float64} = Float64[]

    max_quant_unvax_suscep::Vector{Float64} = Float64[]
    max_quant_vax_suscep::Vector{Float64} = Float64[]

    beta::Vector{Float64} = Float64[]
end

MetricReport(n_steps) = MetricReport(
    times = zeros(Float64, n_steps),
    num_susceptible = zeros(Int64, n_steps),
    num_exposed = zeros(Int64, n_steps),
    num_infected = zeros(Int64, n_steps),
    num_recovered = zeros(Int64, n_steps),
    num_vax_infected = zeros(Int64, n_steps),
    num_unvax_infected = zeros(Int64, n_steps),
    avg_susceptibility = zeros(Float64, n_steps),
    avg_unvax_susceptibility = zeros(Float64, n_steps),
    avg_vax_susceptibility = zeros(Float64, n_steps),
    var_susceptibility = zeros(Float64, n_steps),
    var_unvax_susceptibility = zeros(Float64, n_steps),
    var_vax_susceptibility = zeros(Float64, n_steps),
    min_quant_unvax_suscep = zeros(Float64, n_steps),
    min_quant_vax_suscep = zeros(Float64, n_steps),
    low_quant_unvax_suscep = zeros(Float64, n_steps),
    low_quant_vax_suscep = zeros(Float64, n_steps),
    med_quant_unvax_suscep = zeros(Float64, n_steps),
    med_quant_vax_suscep = zeros(Float64, n_steps),
    upp_quant_unvax_suscep = zeros(Float64, n_steps),
    upp_quant_vax_suscep = zeros(Float64, n_steps),
    max_quant_unvax_suscep = zeros(Float64, n_steps),
    max_quant_vax_suscep = zeros(Float64, n_steps),
    beta = zeros(Float64, n_steps)
)

function generate_report(exp_idx, rep, report::MetricReport)
    df = DataFrame(
        exp = exp_idx,
        rep = rep,
        t = report.times,
        s = report.num_susceptible,
        e = report.num_exposed,
        i = report.num_infected,
        r = report.num_recovered,
        vax_inf = report.num_vax_infected,
        unvax_inf = report.num_unvax_infected,
        tt_avg = report.avg_susceptibility,
        tt_var = report.var_susceptibility,
        us_avg = report.avg_unvax_susceptibility,
        us_var = report.var_unvax_susceptibility,
        vs_avg = report.avg_vax_susceptibility,
        vs_var = report.var_vax_susceptibility,
        min_quant_unvax_suscep = report.min_quant_unvax_suscep,
        min_quant_vax_suscep = report.min_quant_vax_suscep,
        low_quant_unvax_suscep = report.low_quant_unvax_suscep,
        low_quant_vax_suscep = report.low_quant_vax_suscep,
        med_quant_unvax_suscep = report.med_quant_unvax_suscep,
        med_quant_vax_suscep = report.med_quant_vax_suscep,
        upp_quant_unvax_suscep = report.upp_quant_unvax_suscep,
        upp_quant_vax_suscep = report.upp_quant_vax_suscep,
        max_quant_unvax_suscep = report.max_quant_unvax_suscep,
        max_quant_vax_suscep = report.max_quant_vax_suscep,
        beta = report.beta
    )

    return df
end

function safe_shrink_and_summarize!(vec, len, fun)
    resize!(vec, len)
    return length(vec) > 0 ? fun(vec) : -1.0
end

function infected_this_step(p::Person, step, dt)
    t_min = (step - 1) * dt
    t_max = step * dt
    return has_been_infected(p) && last_infection_time(p) >= t_min && last_infection_time(p) < t_max
end

function reporter!(report::MetricReport, par::Parameters, people::Vector{Person}, t)
    step, tnow = t
    report.times[step] = tnow
    report.beta[step] = calculate_beta(par, tnow)

    tot_suscep = zeros(Float64, length(people))
    tot_idx = 0

    vax_suscep = zeros(Float64, length(people))
    vs_idx = 0

    unvax_suscep = zeros(Float64, length(people))
    us_idx = 0
    for p in people
        suscep = current_suscep(par, p, tnow)
        tot_idx += 1
        tot_suscep[tot_idx] = suscep

        if p.vax_status == Vaccinated
            vs_idx += 1
            vax_suscep[vs_idx] = suscep

            report.num_vax_infected[step] += infected_this_step(p, step, par.dt) ? 1 : 0
        else
            us_idx += 1
            unvax_suscep[us_idx] = suscep

            report.num_unvax_infected[step] += infected_this_step(p, step, par.dt) ? 1 : 0
        end

        if p.state == Susceptible
            report.num_susceptible[step] += 1
        elseif p.state == Exposed
            report.num_exposed[step] += 1
        elseif p.state == Infected
            report.num_infected[step] += 1
        else
            report.num_recovered[step] += 1
        end
    end

    report.avg_susceptibility[step] = safe_shrink_and_summarize!(tot_suscep, tot_idx, mean)
    report.avg_unvax_susceptibility[step] = safe_shrink_and_summarize!(unvax_suscep, us_idx, mean)
    report.avg_vax_susceptibility[step] = safe_shrink_and_summarize!(vax_suscep, vs_idx, mean)

    report.var_susceptibility[step] = safe_shrink_and_summarize!(tot_suscep, tot_idx, var)
    report.var_unvax_susceptibility[step] = safe_shrink_and_summarize!(unvax_suscep, us_idx, var)
    report.var_vax_susceptibility[step] = safe_shrink_and_summarize!(vax_suscep, vs_idx, var)

    unvax_susep_quants = length(unvax_suscep) > 0 ?
        quantile(unvax_suscep, [0.1, 0.25, 0.5, 0.75, 0.9]) :
        [-1.0, -1.0, -1.0, -1.0, -1.0]

    vax_susep_quants = length(vax_suscep) > 0 ?
        quantile(vax_suscep, [0.1, 0.25, 0.5, 0.75, 0.9]) :
        [-1.0, -1.0, -1.0, -1.0, -1.0]

    report.min_quant_unvax_suscep[step] = unvax_susep_quants[1]
    report.low_quant_unvax_suscep[step] = unvax_susep_quants[2]
    report.med_quant_unvax_suscep[step] = unvax_susep_quants[3]
    report.upp_quant_unvax_suscep[step] = unvax_susep_quants[4]
    report.max_quant_unvax_suscep[step] = unvax_susep_quants[5]

    report.min_quant_vax_suscep[step] = vax_susep_quants[1]
    report.low_quant_vax_suscep[step] = vax_susep_quants[2]
    report.med_quant_vax_suscep[step] = vax_susep_quants[3]
    report.upp_quant_vax_suscep[step] = vax_susep_quants[4]
    report.max_quant_vax_suscep[step] = vax_susep_quants[5]
end

function vax_coverage(par::Parameters, people::Vector{Person})
    num_vaxd = 0
    for p in people
        num_vaxd += p.vax_status == Vaccinated ? 1 : 0
    end

    return num_vaxd / par.pop_size
end

function should_revaccinate(par::Parameters, time)
    sim_has_started = time > zero(time)
    revaccination_enabled = par.vaccination_checkpoint != 0.0
    checkpoint_reached = (time % par.vaccination_checkpoint) == 0
    vax_campaign_started = time >= par.t_vax_start

    return sim_has_started && revaccination_enabled && checkpoint_reached &&
           vax_campaign_started
end

function save_prevax_suscep(par::Parameters, people::Vector{Person}, time, rep, init)
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
        CSV.write(fullpath, df, append = false)
    else
        CSV.write(fullpath, df, append = true)
    end
end

function simulate(par::Parameters, rep::Real)
    Random.seed!(rep)

    tnow = 0.0
    step = 0
    people = init_pop(par)

    n_steps = Int(par.tmax / par.dt)
    report = MetricReport(n_steps)

    if par.save_prevax_suscep
        save_prevax_suscep(par, people, tnow, rep, true)
    end

    while tnow < par.tmax
        print("Simulating... ", round(Int, (step/n_steps)*100), "%")
        if should_revaccinate(par, tnow)
            revaccinate!(par, people, tnow)
            if par.save_prevax_suscep
                save_prevax_suscep(par, people, tnow, rep, false)
            end
        end

        seed_infections!(par, people, tnow)

        update_foi(par, people, tnow)

        step += 1
        tnow = round(par.dt * step, digits = 5)

        for p in people
            sim_person!(par, p, tnow)
        end

        reporter!(report, par, people, (step, tnow))
        print("\r")
    end

    println("\rSimulating... completed.")

    df = generate_report(par.exp_idx, rep, report)

    return (df = df, people = people)
end

function save_linelist(par::Parameters, rep::Real, people)
    ll = gen_line_list(par.exp_idx, rep, people)

    saveto = par.sim_dir == "default" ? datadir(par.exp_name, "sims") : par.sim_dir
    fname = string("sim_", par.exp_idx, "-", rep, ".csv")
    fullpath = joinpath(saveto, fname)

    mkpath(saveto)
    println("saveing to ", saveto)
    CSV.write(fullpath, ll)
    println("done")

    return ll
end

function save_report(par::Parameters, rep::Real, df)
    saveto = par.sim_dir == "default" ? datadir(par.exp_name, "sims") : par.sim_dir
    fname = string("sim_", par.exp_idx, "-", rep, ".csv")
    fullpath = joinpath(saveto, fname)

    mkpath(saveto)
    println("saveing to ", saveto)
    CSV.write(fullpath, df)
    println("done")

    return df

end