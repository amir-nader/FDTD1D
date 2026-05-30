using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using FDTD1D

function run_case(name, source; nx = 401, dx = 1e-3, nsteps = 1400, source_position = 201)
    config = SimulationConfig(
        nx = nx,
        dx = dx,
        nsteps = nsteps,
        source_position = source_position,
        source = source,
        left_boundary = PECBoundary(),
        right_boundary = PECBoundary(),
        save_every = 4,
    )

    result = run_fdtd(config)
    output = joinpath("outputs", "$(name)_E.gif")
    animate_field(result; field = :E, output = output, title_prefix = name)
    println("Saved animation to $(abspath(output))")
    return result
end

sources = available_source_templates()

run_case("gaussian_source", sources[:gaussian])
run_case("cosine_source", sources[:cosine])
run_case("modulated_cosine_source", sources[:gaussian_modulated_cosine])
run_case("ricker_source", sources[:ricker])
