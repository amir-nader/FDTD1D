using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using FDTD1D
using Plots

function energy_history(result)
    history = zeros(length(result.times))
    for i in eachindex(result.times)
        history[i] = sum(abs2, result.e_history[:, i]) + sum(abs2, result.h_history[:, i])
    end
    return history
end

function source_x(config)
    idx = isnothing(config.source_position) ? fld(config.nx + 1, 2) : config.source_position
    return (idx - 1) * config.dx
end

function monitor_x(config, name::AbstractString)
    for monitor in config.monitors
        if monitor.name == name
            return (monitor.index - 1) * config.dx
        end
    end
    throw(ArgumentError("Monitor '$name' not found."))
end

function right_boundary_reflection_gate_time(config; monitor = "left_probe")
    x_src = source_x(config)
    x_mon = monitor_x(config, monitor)
    boundary_x = (config.nx - 1) * config.dx
    travel_time = (abs(x_src - boundary_x) + abs(x_mon - boundary_x)) / FDTD1D.c0
    pulse_delay = hasproperty(config.source, :t0) ? getproperty(config.source, :t0) : 0.0
    return pulse_delay + travel_time
end

function reflected_peak(trace, times; time_window = nothing)
    values = isnothing(time_window) ? trace : trace[times .>= time_window]
    return isempty(values) ? 0.0 : maximum(abs, values)
end

function gated_signal_energy(trace, times; time_window = nothing, side = :after)
    mask = if isnothing(time_window)
        trues(length(times))
    elseif side === :after
        times .>= time_window
    elseif side === :before
        times .< time_window
    else
        throw(ArgumentError("side must be :before or :after."))
    end
    return sum(abs2, trace[mask])
end

function extend_right_domain_config(config; extra_cells = config.nx)
    material = config.material
    if material isa GridMaterial
        throw(ArgumentError("Boundary comparison reference extension currently supports vacuum material only."))
    end

    return SimulationConfig(
        nx = config.nx + extra_cells,
        dx = config.dx,
        courant_factor = config.courant_factor,
        nsteps = config.nsteps,
        source_position = config.source_position,
        source = config.source,
        excitation = config.excitation,
        left_boundary = config.left_boundary,
        right_boundary = config.right_boundary,
        material = config.material,
        monitors = config.monitors,
        save_every = config.save_every,
    )
end

function extracted_right_boundary_reflection(result, monitor::AbstractString)
    reference = run_fdtd(extend_right_domain_config(result.config))
    signal = result.monitor_traces[monitor]
    reference_signal = reference.monitor_traces[monitor]
    length(signal) == length(reference_signal) ||
        throw(ArgumentError("Reference signal length mismatch in boundary comparison."))
    return signal .- reference_signal
end

function compare_boundary_absorption(
    abc_path::AbstractString,
    pml_path::AbstractString;
    output::AbstractString = "boundary_absorption_compare.png",
)
    abc_params = load_simulation_parameters(abc_path)
    pml_params = load_simulation_parameters(pml_path)

    abc_result = run_fdtd(abc_params.config)
    pml_result = run_fdtd(pml_params.config)

    abc_energy = energy_history(abc_result)
    pml_energy = energy_history(pml_result)

    gate_time = right_boundary_reflection_gate_time(abc_params.config; monitor = "left_probe")
    abc_reflected_trace = extracted_right_boundary_reflection(abc_result, "left_probe")
    pml_reflected_trace = extracted_right_boundary_reflection(pml_result, "left_probe")
    incident_energy = gated_signal_energy(abc_result.monitor_traces["left_probe"], abc_result.times; time_window = gate_time, side = :before)
    abc_peak = reflected_peak(abc_reflected_trace, abc_result.times; time_window = gate_time)
    pml_peak = reflected_peak(pml_reflected_trace, pml_result.times; time_window = gate_time)
    abc_reflected_energy = gated_signal_energy(abc_reflected_trace, abc_result.times; time_window = gate_time, side = :after)
    pml_reflected_energy = gated_signal_energy(pml_reflected_trace, pml_result.times; time_window = gate_time, side = :after)
    abc_reflection_ratio = incident_energy > 0 ? abc_reflected_energy / incident_energy : 0.0
    pml_reflection_ratio = incident_energy > 0 ? pml_reflected_energy / incident_energy : 0.0

    p1 = plot(
        abc_result.times .* 1e9,
        abc_reflected_trace;
        label = "Mur ABC",
        linewidth = 2,
        xlabel = "Time (ns)",
        ylabel = "Extracted reflected E",
        title = "Right-boundary reflection at left probe",
    )
    plot!(p1, pml_result.times .* 1e9, pml_reflected_trace; label = "PML", linewidth = 2)
    vline!(p1, [gate_time * 1e9]; label = "Reflection gate", linestyle = :dash, color = :black)

    p2 = plot(
        abc_result.times .* 1e9,
        abc_energy;
        label = "Mur ABC",
        linewidth = 2,
        xlabel = "Time (ns)",
        ylabel = "Discrete field energy",
        title = "Residual energy in domain",
        yscale = :log10,
    )
    plot!(p2, pml_result.times .* 1e9, pml_energy; label = "PML", linewidth = 2)
    vline!(p2, [gate_time * 1e9]; label = "Reflection gate", linestyle = :dash, color = :black)

    comparison = plot(p1, p2; layout = (2, 1), size = (900, 800))
    mkpath(dirname(abspath(output)))
    savefig(comparison, output)

    println("Saved boundary comparison plot to $(abspath(output))")
    println("Extracted right-boundary reflection comparison after $(gate_time) s:")
    println("  Mur ABC peak |E|: $abc_peak")
    println("  PML peak |E|: $pml_peak")
    println("Reflected energy ratio at left probe:")
    println("  Incident energy proxy: $incident_energy")
    println("  Mur ABC reflected energy: $abc_reflected_energy")
    println("  PML reflected energy: $pml_reflected_energy")
    println("  Mur ABC reflected/incident: $abc_reflection_ratio")
    println("  PML reflected/incident: $pml_reflection_ratio")
    println("Final discrete energy:")
    println("  Mur ABC: $(abc_energy[end])")
    println("  PML: $(pml_energy[end])")

    return (
        abc = abc_result,
        pml = pml_result,
        abc_peak = abc_peak,
        pml_peak = pml_peak,
        incident_energy = incident_energy,
        abc_reflected_energy = abc_reflected_energy,
        pml_reflected_energy = pml_reflected_energy,
        abc_reflection_ratio = abc_reflection_ratio,
        pml_reflection_ratio = pml_reflection_ratio,
    )
end

abc_config = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "config", "abc_vs_pml_abc.toml")
pml_config = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "..", "config", "abc_vs_pml_pml.toml")
output_file = length(ARGS) >= 3 ? ARGS[3] : joinpath("outputs", "boundary_absorption_compare.png")
compare_boundary_absorption(abc_config, pml_config; output = output_file)
