using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "compare_dispersive_slab_analytical.jl"))

function benchmark_dispersive_models(configs::Vector{String})
    summaries = NamedTuple[]
    for config in configs
        println("Running dispersive benchmark for $(config)")
        result = compare_dispersive_slab_analytical(config)
        push!(summaries, (
            config = config,
            model = result.slab.model,
            rmse_reflection = result.rmse_reflection,
            rmse_transmission = result.rmse_transmission,
            max_reflection_error = result.max_reflection_error,
            max_transmission_error = result.max_transmission_error,
        ))
        println()
    end

    println("Benchmark summary:")
    for summary in summaries
        println("  $(summary.model) [$(summary.config)]")
        println("    reflection RMSE = $(summary.rmse_reflection)")
        println("    transmission RMSE = $(summary.rmse_transmission)")
        println("    reflection max error = $(summary.max_reflection_error)")
        println("    transmission max error = $(summary.max_transmission_error)")
    end
    return summaries
end

configs = isempty(ARGS) ? [
    joinpath(@__DIR__, "..", "config", "tfsf_debye_slab_analytical_compare.toml"),
    joinpath(@__DIR__, "..", "config", "tfsf_drude_slab_analytical_compare.toml"),
    joinpath(@__DIR__, "..", "config", "tfsf_lorentz_slab_analytical_compare.toml"),
] : collect(ARGS)
benchmark_dispersive_models(configs)
