using Test
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using FDTD1D

include(joinpath(@__DIR__, "..", "examples", "compare_dispersive_slab_analytical.jl"))

function rms_error(a, b)
    return sqrt(sum(abs2, a .- b) / length(a))
end

function field_energy_history(result)
    return vec(sum(abs2, result.e_history; dims = 1)) .+
           vec(sum(abs2, result.h_history; dims = 1))
end

function dielectric_slab_validation()
    params = load_simulation_parameters(joinpath(@__DIR__, "..", "config", "tfsf_slab_analytical_compare.toml"))
    result = run_fdtd(params.config)
    diagnostics_cfg = params.diagnostics
    frequencies = collect(range(
        Float64(diagnostics_cfg["frequency_min"]);
        stop = Float64(diagnostics_cfg["frequency_max"]),
        length = Int(diagnostics_cfg["frequency_count"]),
    ))

    fdtd = compute_frequency_scattering_diagnostics(
        result;
        frequencies = frequencies,
        incident_monitor = String(diagnostics_cfg["incident_monitor"]),
        reflected_monitor = String(diagnostics_cfg["reflected_monitor"]),
        transmitted_monitor = String(diagnostics_cfg["transmitted_monitor"]),
        window = String(diagnostics_cfg["window"]),
    )
    analytical = analytical_slab_rt(frequencies; slab_eps_r = 4.0, thickness = 0.061)

    return (
        rmse_reflection = rms_error(fdtd.reflection, analytical.reflection),
        rmse_transmission = rms_error(fdtd.transmission, analytical.transmission),
        max_reflection_error = maximum(abs.(fdtd.reflection .- analytical.reflection)),
        max_transmission_error = maximum(abs.(fdtd.transmission .- analytical.transmission)),
    )
end

function halfspace_reflection_validation()
    params = load_simulation_parameters(joinpath(@__DIR__, "..", "config", "tfsf_halfspace_reflection.toml"))
    config = params.config
    monitors = [FieldMonitor("reflected", 61); config.monitors]
    config = SimulationConfig(
        nx = config.nx,
        dx = config.dx,
        courant_factor = config.courant_factor,
        nsteps = config.nsteps,
        source_position = config.source_position,
        source = config.source,
        excitation = config.excitation,
        left_boundary = config.left_boundary,
        right_boundary = config.right_boundary,
        material = config.material,
        monitors = monitors,
        save_every = config.save_every,
    )

    result = run_fdtd(config)
    reflected = result.monitor_traces["reflected"]
    reflected_window = reflected[result.times .>= 1.5e-9]
    peak = reflected_window[argmax(abs.(reflected_window))]

    η1 = sqrt(FDTD1D.mu0 / FDTD1D.eps0)
    η2 = sqrt(FDTD1D.mu0 / (FDTD1D.eps0 * 4.0))
    Γ = (η2 - η1) / (η2 + η1)
    return (measured_peak = peak, analytical_amplitude = Γ)
end

function pml_thickness_validation()
    function config_with_pml_cells(cells)
        return SimulationConfig(
            nx = 401,
            dx = 1e-3,
            courant_factor = 0.99,
            nsteps = 1600,
            source_position = 181,
            source = GaussianCurrent(1.0, 4.0e-11, 1.2e-11),
            excitation = SoftCurrentExcitation(),
            left_boundary = PMLBoundary(ncells = cells, grading_order = 3, target_reflection = 1.0e-6),
            right_boundary = PMLBoundary(ncells = cells, grading_order = 3, target_reflection = 1.0e-6),
            material = Vacuum(),
            monitors = FieldMonitor[],
            save_every = 4,
        )
    end

    thin = run_fdtd(config_with_pml_cells(10))
    thick = run_fdtd(config_with_pml_cells(40))
    thin_energy = field_energy_history(thin)
    thick_energy = field_energy_history(thick)

    return (
        thin_final = thin_energy[end],
        thick_final = thick_energy[end],
        thin_residual_ratio = thin_energy[end] / maximum(thin_energy),
        thick_residual_ratio = thick_energy[end] / maximum(thick_energy),
    )
end

function abc_pml_absorption_validation()
    abc = run_fdtd(load_simulation_parameters(joinpath(@__DIR__, "..", "config", "abc_vs_pml_abc.toml")).config)
    pml = run_fdtd(load_simulation_parameters(joinpath(@__DIR__, "..", "config", "abc_vs_pml_pml.toml")).config)
    abc_energy = field_energy_history(abc)
    pml_energy = field_energy_history(pml)

    return (
        abc_final = abc_energy[end],
        pml_final = pml_energy[end],
        abc_residual_ratio = abc_energy[end] / maximum(abc_energy),
        pml_residual_ratio = pml_energy[end] / maximum(pml_energy),
    )
end

@testset "Nondispersive analytical validations" begin
    @testset "Vacuum TF/SF propagation" begin
        config = SimulationConfig(
            nx = 301,
            dx = 1e-3,
            courant_factor = 0.99,
            nsteps = 800,
            source = GaussianCurrent(1.0, 4.0e-10, 8.0e-11),
            excitation = TFSFExcitation(60, 250, 1.0, 1.0),
            left_boundary = PMLBoundary(ncells = 30),
            right_boundary = PMLBoundary(ncells = 30),
            material = Vacuum(),
            monitors = [FieldMonitor("probe", 160)],
            save_every = 1,
        )
        result = run_fdtd(config)
        expected = FDTD1D.incident_field_trace(config, result.times, "probe")
        measured = result.monitor_traces["probe"]

        @test maximum(abs.(measured .- expected)) < 1.0e-3
        @test rms_error(measured, expected) < 1.0e-4
        @test maximum(abs.(measured)) ≈ 1.0 atol = 1.0e-3
    end

    @testset "Dielectric half-space reflection" begin
        result = halfspace_reflection_validation()

        @test abs(result.measured_peak) ≈ abs(result.analytical_amplitude) atol = 2.0e-2
    end

    @testset "Lossless dielectric slab spectra" begin
        result = dielectric_slab_validation()

        @test result.rmse_reflection < 0.10
        @test result.rmse_transmission < 0.12
        @test result.max_reflection_error < 0.20
        @test result.max_transmission_error < 0.25
    end

    @testset "PML thickness improves residual absorption" begin
        result = pml_thickness_validation()

        @test result.thick_final < result.thin_final
        @test result.thick_final / result.thin_final < 0.75
        @test result.thin_residual_ratio < 5.0e-7
        @test result.thick_residual_ratio < 2.5e-7
    end

    @testset "ABC and PML absorb outgoing waves" begin
        result = abc_pml_absorption_validation()

        @test result.abc_residual_ratio < 5.0e-7
        @test result.pml_residual_ratio < 2.5e-7
        @test result.pml_final < result.abc_final
        @test result.pml_final / result.abc_final < 0.80
    end
end

@testset "S-parameter port metadata" begin
    params = load_simulation_parameters(joinpath(@__DIR__, "..", "config", "tfsf_slab_sparameters.toml"))
    config = params.config
    short_config = SimulationConfig(
        nx = config.nx,
        dx = config.dx,
        courant_factor = config.courant_factor,
        nsteps = 500,
        source_position = config.source_position,
        source = config.source,
        excitation = config.excitation,
        left_boundary = config.left_boundary,
        right_boundary = config.right_boundary,
        material = config.material,
        monitors = config.monitors,
        save_every = 2,
    )
    result = run_fdtd(short_config)
    sparams = compute_sparameters(
        result;
        frequencies = [1.0e9, 1.5e9],
        port1_monitor = "incident",
        port1_reflected_monitor = "reflected",
        port2_monitor = "transmitted",
        port1_reference_plane = 0.18,
        port2_reference_plane = 0.24,
        window = "hann",
    )

    @test sparams.metadata.port1_monitor == "incident"
    @test sparams.metadata.port1_reflected_monitor == "reflected"
    @test sparams.metadata.port2_monitor == "transmitted"
    @test sparams.metadata.port1_reference_plane ≈ 0.18
    @test sparams.metadata.port2_reference_plane ≈ 0.24
    @test sparams.left_shift ≈ 0.05
    @test sparams.right_shift ≈ 0.05

    path = joinpath(mktempdir(), "metadata_check.s2p")
    write_touchstone_s2p(sparams, path)
    contents = read(path, String)
    @test occursin("! port1_monitor: incident", contents)
    @test occursin("! port2_reference_plane: 0.24", contents)
end

@testset "Output manager artifacts" begin
    manager = create_output_manager(nothing; root = mktempdir(), case_name = "manager test", timestamped = false)
    config = SimulationConfig(
        nx = 101,
        dx = 1e-3,
        nsteps = 20,
        source = GaussianCurrent(1.0, 4.0e-11, 1.2e-11),
        source_position = 51,
        monitors = [FieldMonitor("center", 51)],
        save_every = 2,
    )
    result = run_fdtd(config)

    monitors_file = write_monitor_traces_csv(joinpath(manager.run_dir, "monitor_traces.csv"), result)
    material_file = write_material_profile_csv(joinpath(manager.run_dir, "material_profile.csv"), config)
    summary_file = write_run_summary(
        manager,
        result;
        files = Dict("monitor_traces" => monitors_file, "material_profile" => material_file),
    )

    @test isdir(manager.run_dir)
    @test basename(manager.run_dir) == "manager_test"
    @test isfile(monitors_file)
    @test isfile(material_file)
    @test isfile(summary_file)
    @test startswith(readline(monitors_file), "time_s,center")
    @test startswith(readline(material_file), "index,x_m,eps_r_e")
    @test occursin("case_name = \"manager_test\"", read(summary_file, String))
end

const DISPERSIVE_BASELINES = Dict(
    "config/tfsf_debye_slab_analytical_compare.toml" => (
        rmse_reflection = 0.045015358048504754,
        rmse_transmission = 0.02830345876497214,
        max_reflection_error = 0.13018779963927923,
        max_transmission_error = 0.09768822195628435,
    ),
    "config/tfsf_drude_slab_analytical_compare.toml" => (
        rmse_reflection = 0.15435432213167338,
        rmse_transmission = 0.032856244603649276,
        max_reflection_error = 0.28650661113988296,
        max_transmission_error = 0.09436367084229713,
    ),
    "config/tfsf_lorentz_slab_analytical_compare.toml" => (
        rmse_reflection = 0.16574935409428881,
        rmse_transmission = 0.0668779319482676,
        max_reflection_error = 0.34124390212956146,
        max_transmission_error = 0.26753968773522363,
    ),
)

@testset "Dispersive analytical regressions" begin
    for (relpath, baseline) in DISPERSIVE_BASELINES
        config = joinpath(@__DIR__, "..", relpath)
        result = compare_dispersive_slab_analytical(config)

        @test isfinite(result.rmse_reflection)
        @test isfinite(result.rmse_transmission)
        @test isfinite(result.max_reflection_error)
        @test isfinite(result.max_transmission_error)

        @test result.rmse_reflection ≈ baseline.rmse_reflection rtol = 0.25 atol = 1e-3
        @test result.rmse_transmission ≈ baseline.rmse_transmission rtol = 0.25 atol = 1e-3
        @test result.max_reflection_error ≈ baseline.max_reflection_error rtol = 0.30 atol = 1e-3
        @test result.max_transmission_error ≈ baseline.max_transmission_error rtol = 0.30 atol = 1e-3

        @test result.rmse_reflection < 0.25
        @test result.rmse_transmission < 0.10
        @test result.max_reflection_error < 0.40
        @test result.max_transmission_error < 0.30
    end
end
