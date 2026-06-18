using DrWatson
@quickactivate :PreVaxImmunityVsMeasuredVE

# config deps
using TOML, ArgParse
# data manip deps
using DataFramesMeta, CSV

function parse_cmdline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--config"
            required = true
            help = "path to config TOML file (abs path or rel path from this script)"
        "--clean"
            default = true
            arg_type = Bool
            help = "clean up (delete) sim datafiles after collection into results.csv"
    end

    return parse_args(s)
end

function collect_simulation_data(sim_dir)
    csv_vec = readdir(sim_dir)
    results = Vector{DataFrame}(undef, length(csv_vec))

    for i in eachindex(csv_vec)
        f = joinpath(sim_dir, csv_vec[i])
        results[i] = DataFrame(CSV.File(f))
    end

    return results
end

function main()
    # parse commandline arguments
    parsed_args = parse_cmdline()

    delete_after_read = parsed_args["clean"]

    usr_path = parsed_args["config"]
    config_path = isabspath(usr_path) ? usr_path : abspath(usr_path)

    # parse config
    config = TOML.parsefile(config_path)
    all_pars = config["parameters"]
    exp_name = all_pars["exp_name"]

    println("slurping sim csv files")
    sim_dir = all_pars["sim_dir"] == "default" ? datadir(exp_name, "sims") : all_pars["sim_dir"]
    results = collect_simulation_data(sim_dir)

    println("writing results.csv")
    mkpath(datadir(exp_name))
    CSV.write(datadir(exp_name, "results.csv"), vcat(results...))

    if delete_after_read
        println("cleaning sim data dir")
        rm(sim_dir, recursive = true)
    end
end

main()
