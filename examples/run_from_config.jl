using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using FDTD1D

function optional_plane_from_cfg(cfg, key_x, key_index, dx)
    if haskey(cfg, key_x) && haskey(cfg, key_index)
        throw(ArgumentError("Specify only one of '$key_x' or '$key_index'."))
    elseif haskey(cfg, key_x)
        return Float64(cfg[key_x])
    elseif haskey(cfg, key_index)
        return (Int(cfg[key_index]) - 1) * dx
    end
    return nothing
end

function optional_string_from_cfg(cfg, keys...)
    for key in keys
        if haskey(cfg, key)
            return String(cfg[key])
        end
    end
    return nothing
end

function run_from_config(path::AbstractString)
    params = load_simulation_parameters(path)
    config = params.config
    output_cfg = params.output
    diagnostics_cfg = params.diagnostics
    use_run_directory = Bool(get(output_cfg, "use_run_directory", true))
    manager = use_run_directory ? create_output_manager(
        path;
        root = String(get(output_cfg, "output_root", "outputs")),
        case_name = String(get(output_cfg, "case_name", splitext(basename(path))[1])),
        timestamped = Bool(get(output_cfg, "timestamped", true)),
    ) : nothing
    managed_path(requested; default_name) =
        isnothing(manager) ? String(requested) : output_path(manager, String(requested); default_name = default_name)

    files = Dict{String,Any}()
    if !isnothing(manager)
        files["input_config"] = copy_input_config!(manager)
    end

    result = run_fdtd(config)

    field = Symbol(get(output_cfg, "field", "E"))
    output_file = managed_path(
        get(output_cfg, "file", joinpath("outputs", "fdtd_1d.gif"));
        default_name = "animation.gif",
    )
    fps = Int(get(output_cfg, "fps", 20))
    title_prefix = String(get(output_cfg, "title_prefix", "1D FDTD"))
    show_material = Bool(get(output_cfg, "show_material", true))
    label_materials = Bool(get(output_cfg, "label_materials", true))
    label_boundaries = Bool(get(output_cfg, "label_boundaries", true))
    diagnostics = nothing
    spectrum = nothing
    sparams = nothing

    animate_field(
        result;
        field = field,
        output = output_file,
        fps = fps,
        title_prefix = title_prefix,
        show_material = show_material,
        label_materials = label_materials,
        label_boundaries = label_boundaries,
    )
    files["animation"] = output_file
    if !isnothing(manager)
        files["monitor_traces"] = write_monitor_traces_csv(joinpath(manager.run_dir, "monitor_traces.csv"), result)
        files["material_profile"] = write_material_profile_csv(joinpath(manager.run_dir, "material_profile.csv"), config)
    end

    if Bool(get(diagnostics_cfg, "enabled", false))
        diagnostics = compute_scattering_diagnostics(
            result;
            incident_monitor = String(get(diagnostics_cfg, "incident_monitor", "incident")),
            reflected_monitor = String(get(diagnostics_cfg, "reflected_monitor", "incident")),
            transmitted_monitor = String(get(diagnostics_cfg, "transmitted_monitor", "transmitted")),
        )
        reflection = diagnostics["reflection_coefficient"]
        transmission = diagnostics["transmission_coefficient"]
        println("Scattering diagnostics:")
        println("  Reflection coefficient: $reflection")
        println("  Transmission coefficient: $transmission")

        if haskey(diagnostics_cfg, "frequency_min") && haskey(diagnostics_cfg, "frequency_max")
            frequency_count = Int(get(diagnostics_cfg, "frequency_count", 101))
            frequencies = collect(range(
                Float64(diagnostics_cfg["frequency_min"]);
                stop = Float64(diagnostics_cfg["frequency_max"]),
                length = frequency_count,
            ))
            spectrum = compute_frequency_scattering_diagnostics(
                result;
                frequencies = frequencies,
                incident_monitor = String(get(diagnostics_cfg, "incident_monitor", "incident")),
                reflected_monitor = String(get(diagnostics_cfg, "reflected_monitor", "incident")),
                transmitted_monitor = String(get(diagnostics_cfg, "transmitted_monitor", "transmitted")),
                window = String(get(diagnostics_cfg, "window", "hann")),
                gate_start = haskey(diagnostics_cfg, "gate_start") ? Float64(diagnostics_cfg["gate_start"]) : nothing,
                gate_end = haskey(diagnostics_cfg, "gate_end") ? Float64(diagnostics_cfg["gate_end"]) : nothing,
            )
            mid = cld(length(frequencies), 2)
            println("  FFT window: $(spectrum.window)")
            if !isnothing(spectrum.gate_start) || !isnothing(spectrum.gate_end)
                println("  FFT gate: [$(spectrum.gate_start), $(spectrum.gate_end)] s")
            end
            println("  Frequency sample at $(frequencies[mid]) Hz:")
            println("    R(f): $(spectrum.reflection[mid])")
            println("    T(f): $(spectrum.transmission[mid])")
            if !isnothing(manager)
                files["spectrum_csv"] = write_spectrum_csv(joinpath(manager.run_dir, "spectrum.csv"), spectrum)
            end

            if Bool(get(diagnostics_cfg, "sparameters", false))
                plot_file = managed_path(
                    get(diagnostics_cfg, "sparameter_plot", joinpath("outputs", "sparameters.png"));
                    default_name = "sparameters.png",
                )
                s2p_file = managed_path(
                    get(diagnostics_cfg, "s2p_file", joinpath("outputs", "sparameters.s2p"));
                    default_name = "sparameters.s2p",
                )
                title = String(get(diagnostics_cfg, "sparameter_title", "S-parameters"))
                reference_impedance = Float64(get(diagnostics_cfg, "s2p_reference_impedance", 50.0))
                touchstone_format = String(get(diagnostics_cfg, "touchstone_format", "RI"))
                touchstone_ports = Int(get(diagnostics_cfg, "touchstone_ports", 2))
                port1_monitor = optional_string_from_cfg(diagnostics_cfg, "port1_monitor", "input_monitor")
                port1_reflected_monitor = optional_string_from_cfg(diagnostics_cfg, "port1_reflected_monitor", "reflected_monitor")
                port2_monitor = optional_string_from_cfg(diagnostics_cfg, "port2_monitor", "output_monitor")
                port1_reference_plane = optional_plane_from_cfg(diagnostics_cfg, "port1_reference_plane_x", "port1_reference_plane_index", config.dx)
                port2_reference_plane = optional_plane_from_cfg(diagnostics_cfg, "port2_reference_plane_x", "port2_reference_plane_index", config.dx)
                left_reference_plane = optional_plane_from_cfg(diagnostics_cfg, "left_reference_plane_x", "left_reference_plane_index", config.dx)
                right_reference_plane = optional_plane_from_cfg(diagnostics_cfg, "right_reference_plane_x", "right_reference_plane_index", config.dx)
                sparams = compute_sparameters(
                    result;
                    frequencies = frequencies,
                    incident_monitor = String(get(diagnostics_cfg, "incident_monitor", "incident")),
                    reflected_monitor = String(get(diagnostics_cfg, "reflected_monitor", "incident")),
                    transmitted_monitor = String(get(diagnostics_cfg, "transmitted_monitor", "transmitted")),
                    port1_monitor = port1_monitor,
                    port1_reflected_monitor = port1_reflected_monitor,
                    port2_monitor = port2_monitor,
                    window = String(get(diagnostics_cfg, "window", "hann")),
                    gate_start = haskey(diagnostics_cfg, "gate_start") ? Float64(diagnostics_cfg["gate_start"]) : nothing,
                    gate_end = haskey(diagnostics_cfg, "gate_end") ? Float64(diagnostics_cfg["gate_end"]) : nothing,
                    left_reference_plane = left_reference_plane,
                    right_reference_plane = right_reference_plane,
                    port1_reference_plane = port1_reference_plane,
                    port2_reference_plane = port2_reference_plane,
                )
                plot_sparameters(sparams; output = plot_file, title_prefix = title)
                files["sparameter_plot"] = plot_file
                if touchstone_ports == 1
                    write_touchstone_s1p(
                        sparams,
                        s2p_file;
                        parameter = Symbol(get(diagnostics_cfg, "s1p_parameter", "s11")),
                        format = touchstone_format,
                        reference_impedance = reference_impedance,
                    )
                else
                    write_touchstone_s2p(
                        sparams,
                        s2p_file;
                        format = touchstone_format,
                        reference_impedance = reference_impedance,
                    )
                end
                files["touchstone"] = s2p_file
                if !isnothing(manager)
                    files["sparameters_csv"] = write_sparameters_csv(joinpath(manager.run_dir, "sparameters.csv"), sparams)
                end
                println("  S-parameter sample at $(frequencies[mid]) Hz:")
                println("    S11: $(sparams.s11[mid])")
                println("    S21: $(sparams.s21[mid])")
                println("    S12: $(sparams.s12[mid])")
                println("    S22: $(sparams.s22[mid])")
                println("  Reference-plane shifts:")
                println("    Left shift: $(haskey(sparams, :left_shift) ? sparams.left_shift : 0.0) m")
                println("    Right shift: $(haskey(sparams, :right_shift) ? sparams.right_shift : 0.0) m")
                if haskey(sparams, :metadata)
                    println("  Port definitions:")
                    println("    Port 1 monitor: $(sparams.metadata.port1_monitor), reflected monitor: $(sparams.metadata.port1_reflected_monitor)")
                    println("    Port 2 monitor: $(sparams.metadata.port2_monitor)")
                    println("    Port 1 reference plane: $(sparams.metadata.port1_reference_plane) m")
                    println("    Port 2 reference plane: $(sparams.metadata.port2_reference_plane) m")
                end
                if Bool(get(diagnostics_cfg, "sparameter_delays", false))
                    delay_parameter = Symbol(get(diagnostics_cfg, "sparameter_delay_parameter", "s21"))
                    delays = compute_sparameter_delays(sparams; parameter = delay_parameter)
                    delay_plot = managed_path(
                        get(diagnostics_cfg, "sparameter_delay_plot", joinpath("outputs", "sparameter_delays.png"));
                        default_name = "sparameter_delays.png",
                    )
                    delay_title = String(get(diagnostics_cfg, "sparameter_delay_title", "S-parameter delays"))
                    plot_sparameter_delays(delays; output = delay_plot, title_prefix = delay_title)
                    files["sparameter_delay_plot"] = delay_plot
                    println("  Delay sample at $(frequencies[mid]) Hz:")
                    println("    Phase delay ($(uppercase(String(delay_parameter)))): $(delays.phase_delay[mid]) s")
                    println("    Group delay ($(uppercase(String(delay_parameter)))): $(delays.group_delay[mid]) s")
                    println("  Saved delay plot to $(abspath(delay_plot))")
                end
                println("  Saved S-parameter plot to $(abspath(plot_file))")
                println("  Saved Touchstone file to $(abspath(s2p_file))")
            end
        end
    end

    if !isnothing(manager)
        files["summary"] = write_run_summary(manager, result; diagnostics = diagnostics, spectrum = spectrum, sparams = sparams, files = files)
        summary_file = files["summary"]
        println("Run directory: $(abspath(manager.run_dir))")
        println("Saved run summary to $(abspath(summary_file))")
    end
    println("Saved animation to $(abspath(output_file))")
    return result
end

config_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "config", "default.toml")
run_from_config(config_path)
