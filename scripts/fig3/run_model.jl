using DrWatson
@quickactivate :PreVaxImmunityVsMeasuredVE

using TOML, ArgParse

function parse_cmdline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--config"
            required = true
            help = "path to config TOML file (abs path or rel path from this script)"
    end

    return parse_args(s)
end

function generate_experiment_csv(par_sets, exp_name)
    experiments = Vector{DataFrame}(undef, length(par_sets))
    for (i, d) in enumerate(par_sets)
        d["exp_idx"] = i
        experiments[i] = DataFrame(d)
    end

    mkpath(datadir(exp_name))
    CSV.write(datadir(exp_name, "experiments.csv"), vcat(experiments...))
end

function generate_parameter_sets(all_pars)
    par_sets = dict_list(all_pars)
    return par_sets
end

function main()
    # parse commandline arguments
    parsed_args = parse_cmdline()
    usr_path = parsed_args["config"]
    config_path = isabspath(usr_path) ? usr_path : abspath(usr_path)

    # parse config and divide into jobs
    config = TOML.parsefile(config_path)
    all_pars = config["parameters"]

    par_sets = generate_parameter_sets(all_pars)
    generate_experiment_csv(par_sets, all_pars["exp_name"])

    if length(par_sets) > 1
        println("ERROR: more than one scenario is parameterized in TOML file")
        exit()
    end

    pars = Parameters(par_sets[1])
    for rep in 1:pars.nrep
        df, people = simulate(pars, rep)
        if pars.save_linelist
            save_linelist(pars, rep, people)
        else
            save_report(pars, rep, df)
        end
    end
end

main()