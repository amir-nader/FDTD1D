using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using FDTD1D
using Plots

function halfspace_from_config(config)
    material = config.material
    material isa GridMaterial || throw(ArgumentError("Half-space comparison currently requires GridMaterial."))

    start_index = findfirst(!=(1.0), material.eps_r_e)
    isnothing(start_index) && throw(ArgumentError("No dielectric half-space found."))

    eps_r = material.eps_r_e[start_index]
    mu_r = start_index <= length(material.mu_r_h) ? material.mu_r_h[min(start_index, length(material.mu_r_h))] : 1.0

    all(i -> isapprox(material.eps_r_e[i], eps_r; atol = 1e-12), start_index:length(material.eps_r_e)) ||
        throw(ArgumentError("Half-space example expects a constant dielectric from the interface to the right edge."))

    return (
        start_index = start_index,
        interface_x = (start_index - 1) * config.dx,
        eps_r = eps_r,
        mu_r = mu_r,
    )
end

function analytical_halfspace_reflection(frequencies; eps_r::Real, mu_r::Real = 1.0, incident_eps_r::Real = 1.0, incident_mu_r::Real = 1.0)
    η1 = sqrt(FDTD1D.mu0 * incident_mu_r / (FDTD1D.eps0 * incident_eps_r))
    η2 = sqrt(FDTD1D.mu0 * mu_r / (FDTD1D.eps0 * eps_r))
    Γ = (η2 - η1) / (η2 + η1)
    return (
        frequencies = collect(frequencies),
        reflection_amplitude = fill(Γ, length(frequencies)),
        reflection = fill(abs2(Γ), length(frequencies)),
    )
end

function compare_halfspace_reflection(path::AbstractString)
    params = load_simulation_parameters(path)
    result = run_fdtd(params.config)
    diagnostics_cfg = params.diagnostics

    frequencies = collect(range(
        Float64(get(diagnostics_cfg, "frequency_min", 2.0e8));
        stop = Float64(get(diagnostics_cfg, "frequency_max", 3.0e9)),
        length = Int(get(diagnostics_cfg, "frequency_count", 121)),
    ))

    fdtd = compute_frequency_scattering_diagnostics(
        result;
        frequencies = frequencies,
        incident_monitor = String(get(diagnostics_cfg, "incident_monitor", "incident")),
        reflected_monitor = String(get(diagnostics_cfg, "reflected_monitor", "incident")),
        transmitted_monitor = String(get(diagnostics_cfg, "transmitted_monitor", "transmitted")),
        window = String(get(diagnostics_cfg, "window", "hann")),
        gate_start = haskey(diagnostics_cfg, "gate_start") ? Float64(diagnostics_cfg["gate_start"]) : nothing,
        gate_end = haskey(diagnostics_cfg, "gate_end") ? Float64(diagnostics_cfg["gate_end"]) : nothing,
    )

    halfspace = halfspace_from_config(params.config)
    analytical = analytical_halfspace_reflection(
        frequencies;
        eps_r = halfspace.eps_r,
        mu_r = halfspace.mu_r,
        incident_eps_r = params.config.excitation.incident_eps_r,
        incident_mu_r = params.config.excitation.incident_mu_r,
    )

    output = String(get(diagnostics_cfg, "comparison_plot", joinpath("outputs", "halfspace_reflection_fdtd_vs_analytical.png")))
    plot(
        frequencies ./ 1e9,
        fdtd.reflection;
        label = "FDTD reflection",
        xlabel = "Frequency (GHz)",
        ylabel = "Reflection coefficient",
        linewidth = 2,
    )
    plot!(frequencies ./ 1e9, analytical.reflection; label = "Analytical reflection", linewidth = 2, linestyle = :dash)
    mkpath(dirname(abspath(output)))
    savefig(output)

    center = cld(length(frequencies), 2)
    println("Saved half-space reflection plot to $(abspath(output))")
    println("Interface at x = $(halfspace.interface_x) m")
    println("Half-space eps_r=$(halfspace.eps_r), mu_r=$(halfspace.mu_r)")
    println("FFT window: $(fdtd.window)")
    if !isnothing(fdtd.gate_start) || !isnothing(fdtd.gate_end)
        println("FFT gate: [$(fdtd.gate_start), $(fdtd.gate_end)] s")
    end
    println("At $(frequencies[center]) Hz:")
    println("  FDTD reflection = $(fdtd.reflection[center])")
    println("  Analytical reflection = $(analytical.reflection[center])")
    return (fdtd = fdtd, analytical = analytical, result = result, halfspace = halfspace)
end

config_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "config", "tfsf_halfspace_reflection.toml")
compare_halfspace_reflection(config_path)
