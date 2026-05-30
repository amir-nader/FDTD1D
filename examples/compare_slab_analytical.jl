using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using FDTD1D
using Plots

function first_slab_from_config(config)
    material = config.material
    material isa GridMaterial || throw(ArgumentError("Analytical slab comparison requires GridMaterial."))

    start_index = findfirst(!=(1.0), material.eps_r_e)
    isnothing(start_index) && throw(ArgumentError("No non-vacuum slab found."))

    stop_index = start_index
    while stop_index < length(material.eps_r_e) && material.eps_r_e[stop_index + 1] == material.eps_r_e[start_index]
        stop_index += 1
    end

    return (
        start_index = start_index,
        stop_index = stop_index,
        eps_r = material.eps_r_e[start_index],
        mu_r = stop_index > start_index ? material.mu_r_h[start_index] : 1.0,
        thickness = (stop_index - start_index + 1) * config.dx,
    )
end

function compare_slab_analytical(path::AbstractString)
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

    slab = first_slab_from_config(params.config)
    analytical = analytical_slab_rt(
        frequencies;
        slab_eps_r = slab.eps_r,
        slab_mu_r = slab.mu_r,
        thickness = slab.thickness,
        incident_eps_r = params.config.excitation.incident_eps_r,
        incident_mu_r = params.config.excitation.incident_mu_r,
    )

    output = String(get(diagnostics_cfg, "comparison_plot", joinpath("outputs", "slab_fdtd_vs_analytical.png")))
    plot(
        frequencies ./ 1e9,
        fdtd.reflection;
        label = "FDTD R",
        xlabel = "Frequency (GHz)",
        ylabel = "Power coefficient",
        linewidth = 2,
    )
    plot!(frequencies ./ 1e9, analytical.reflection; label = "Analytical R", linewidth = 2, linestyle = :dash)
    plot!(frequencies ./ 1e9, fdtd.transmission; label = "FDTD T", linewidth = 2)
    plot!(frequencies ./ 1e9, analytical.transmission; label = "Analytical T", linewidth = 2, linestyle = :dash)
    mkpath(dirname(abspath(output)))
    savefig(output)

    center = cld(length(frequencies), 2)
    println("Saved comparison plot to $(abspath(output))")
    println("FFT window: $(fdtd.window)")
    if !isnothing(fdtd.gate_start) || !isnothing(fdtd.gate_end)
        println("FFT gate: [$(fdtd.gate_start), $(fdtd.gate_end)] s")
    end
    println("Slab eps_r=$(slab.eps_r), mu_r=$(slab.mu_r), thickness=$(slab.thickness) m")
    println("At $(frequencies[center]) Hz:")
    println("  FDTD R=$(fdtd.reflection[center]), analytical R=$(analytical.reflection[center])")
    println("  FDTD T=$(fdtd.transmission[center]), analytical T=$(analytical.transmission[center])")
    return (fdtd = fdtd, analytical = analytical, result = result)
end

config_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "config", "tfsf_slab_analytical_compare.toml")
compare_slab_analytical(config_path)
