using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using FDTD1D
using Plots

function first_dispersive_slab_from_config(config)
    material = config.material

    if material isa DebyeMaterial
        start_index = findfirst(i -> material.eps_static_e[i] != 1.0 || material.eps_inf_e[i] != 1.0, eachindex(material.eps_inf_e))
        isnothing(start_index) && throw(ArgumentError("No Debye slab found."))
        stop_index = start_index
        while stop_index < length(material.eps_inf_e) &&
              material.eps_inf_e[stop_index + 1] == material.eps_inf_e[start_index] &&
              material.eps_static_e[stop_index + 1] == material.eps_static_e[start_index] &&
              material.tau_e[stop_index + 1] == material.tau_e[start_index]
            stop_index += 1
        end
        return (
            model = :debye,
            start_index = start_index,
            stop_index = stop_index,
            thickness = (stop_index - start_index + 1) * config.dx,
            eps_inf = material.eps_inf_e[start_index],
            eps_static = material.eps_static_e[start_index],
            tau = material.tau_e[start_index],
            mu_r = stop_index > start_index ? material.mu_r_h[start_index] : 1.0,
        )
    elseif material isa DrudeMaterial
        start_index = findfirst(i -> material.omega_p_e[i] != 0.0 || material.eps_inf_e[i] != 1.0, eachindex(material.eps_inf_e))
        isnothing(start_index) && throw(ArgumentError("No Drude slab found."))
        stop_index = start_index
        while stop_index < length(material.eps_inf_e) &&
              material.eps_inf_e[stop_index + 1] == material.eps_inf_e[start_index] &&
              material.omega_p_e[stop_index + 1] == material.omega_p_e[start_index] &&
              material.gamma_e[stop_index + 1] == material.gamma_e[start_index]
            stop_index += 1
        end
        return (
            model = :drude,
            start_index = start_index,
            stop_index = stop_index,
            thickness = (stop_index - start_index + 1) * config.dx,
            eps_inf = material.eps_inf_e[start_index],
            omega_p = material.omega_p_e[start_index],
            gamma = material.gamma_e[start_index],
            mu_r = stop_index > start_index ? material.mu_r_h[start_index] : 1.0,
        )
    elseif material isa LorentzMaterial
        start_index = findfirst(i -> material.delta_eps_e[i] != 0.0 || material.eps_inf_e[i] != 1.0, eachindex(material.eps_inf_e))
        isnothing(start_index) && throw(ArgumentError("No Lorentz slab found."))
        stop_index = start_index
        while stop_index < length(material.eps_inf_e) &&
              material.eps_inf_e[stop_index + 1] == material.eps_inf_e[start_index] &&
              material.delta_eps_e[stop_index + 1] == material.delta_eps_e[start_index] &&
              material.omega_0_e[stop_index + 1] == material.omega_0_e[start_index] &&
              material.gamma_e[stop_index + 1] == material.gamma_e[start_index]
            stop_index += 1
        end
        return (
            model = :lorentz,
            start_index = start_index,
            stop_index = stop_index,
            thickness = (stop_index - start_index + 1) * config.dx,
            eps_inf = material.eps_inf_e[start_index],
            delta_eps = material.delta_eps_e[start_index],
            omega_0 = material.omega_0_e[start_index],
            gamma = material.gamma_e[start_index],
            mu_r = stop_index > start_index ? material.mu_r_h[start_index] : 1.0,
        )
    end

    throw(ArgumentError("Dispersive analytical comparison requires DebyeMaterial, DrudeMaterial, or LorentzMaterial."))
end

function compare_dispersive_slab_analytical(path::AbstractString)
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

    slab = first_dispersive_slab_from_config(params.config)
    analytical = analytical_dispersive_slab_rt(
        frequencies;
        model = slab.model,
        thickness = slab.thickness,
        slab_mu_r = slab.mu_r,
        incident_eps_r = params.config.excitation.incident_eps_r,
        incident_mu_r = params.config.excitation.incident_mu_r,
        eps_inf = haskey(slab, :eps_inf) ? slab.eps_inf : 1.0,
        eps_static = haskey(slab, :eps_static) ? slab.eps_static : nothing,
        tau = haskey(slab, :tau) ? slab.tau : nothing,
        omega_p = haskey(slab, :omega_p) ? slab.omega_p : nothing,
        gamma = haskey(slab, :gamma) ? slab.gamma : nothing,
        delta_eps = haskey(slab, :delta_eps) ? slab.delta_eps : nothing,
        omega_0 = haskey(slab, :omega_0) ? slab.omega_0 : nothing,
    )

    abs_reflection_error = abs.(fdtd.reflection .- analytical.reflection)
    abs_transmission_error = abs.(fdtd.transmission .- analytical.transmission)
    rmse_reflection = sqrt(sum(abs2, abs_reflection_error) / length(abs_reflection_error))
    rmse_transmission = sqrt(sum(abs2, abs_transmission_error) / length(abs_transmission_error))
    max_reflection_error = maximum(abs_reflection_error)
    max_transmission_error = maximum(abs_transmission_error)

    output = String(get(diagnostics_cfg, "comparison_plot", joinpath("outputs", "dispersive_slab_fdtd_vs_analytical.png")))
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
    println("Model: $(slab.model), thickness=$(slab.thickness) m")
    println("FFT window: $(fdtd.window)")
    if !isnothing(fdtd.gate_start) || !isnothing(fdtd.gate_end)
        println("FFT gate: [$(fdtd.gate_start), $(fdtd.gate_end)] s")
    end
    println("At $(frequencies[center]) Hz:")
    println("  FDTD R=$(fdtd.reflection[center]), analytical R=$(analytical.reflection[center])")
    println("  FDTD T=$(fdtd.transmission[center]), analytical T=$(analytical.transmission[center])")
    println("Spectral error summary:")
    println("  Reflection RMSE: $rmse_reflection")
    println("  Transmission RMSE: $rmse_transmission")
    println("  Reflection max abs error: $max_reflection_error")
    println("  Transmission max abs error: $max_transmission_error")
    return (
        fdtd = fdtd,
        analytical = analytical,
        result = result,
        slab = slab,
        rmse_reflection = rmse_reflection,
        rmse_transmission = rmse_transmission,
        max_reflection_error = max_reflection_error,
        max_transmission_error = max_transmission_error,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    config_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "config", "tfsf_debye_slab_analytical_compare.toml")
    compare_dispersive_slab_analytical(config_path)
end
