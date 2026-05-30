module FDTD1D

using FFTW
using Plots
using TOML

export AbstractMaterial,
       Vacuum,
       GridMaterial,
       DebyeMaterial,
       DrudeMaterial,
       LorentzMaterial,
       MaterialLayer,
       AbstractExcitation,
       SoftCurrentExcitation,
       TFSFExcitation,
       FieldMonitor,
       SimulationConfig,
       PECBoundary,
       MurABCBoundary,
       PMLBoundary,
       NoBoundaryCondition,
       GaussianCurrent,
       CosineCurrent,
       GaussianModulatedCosineCurrent,
       RickerCurrent,
       SimulationResult,
       available_source_templates,
       load_simulation_parameters,
       compute_scattering_diagnostics,
       compute_frequency_scattering_diagnostics,
       compute_sparameters,
       compute_sparameter_delays,
       analytical_slab_rt,
       analytical_dispersive_slab_rt,
       source_from_dict,
       boundary_from_name,
       OutputManager,
       create_output_manager,
       output_path,
       copy_input_config!,
       write_run_summary,
       write_monitor_traces_csv,
       write_material_profile_csv,
       write_spectrum_csv,
       write_sparameters_csv,
       run_fdtd,
       animate_field,
       plot_sparameter_delays,
       plot_sparameters,
       write_touchstone_s1p,
       write_touchstone_s2p

const c0 = 299_792_458.0
const mu0 = 4e-7 * pi
const eps0 = 1.0 / (mu0 * c0^2)

abstract type AbstractMaterial end
struct Vacuum <: AbstractMaterial end

struct MaterialLayer
    start_index::Int
    end_index::Int
    eps_r::Float64
    mu_r::Float64
    sigma_e::Float64
end

struct GridMaterial <: AbstractMaterial
    eps_r_e::Vector{Float64}
    mu_r_h::Vector{Float64}
    sigma_e::Vector{Float64}
end

struct DebyeMaterial <: AbstractMaterial
    eps_inf_e::Vector{Float64}
    eps_static_e::Vector{Float64}
    tau_e::Vector{Float64}
    mu_r_h::Vector{Float64}
    sigma_e::Vector{Float64}
end

struct DrudeMaterial <: AbstractMaterial
    eps_inf_e::Vector{Float64}
    omega_p_e::Vector{Float64}
    gamma_e::Vector{Float64}
    mu_r_h::Vector{Float64}
    sigma_e::Vector{Float64}
end

struct LorentzMaterial <: AbstractMaterial
    eps_inf_e::Vector{Float64}
    delta_eps_e::Vector{Float64}
    omega_0_e::Vector{Float64}
    gamma_e::Vector{Float64}
    mu_r_h::Vector{Float64}
    sigma_e::Vector{Float64}
end

abstract type BoundaryCondition end
struct PECBoundary <: BoundaryCondition end
struct MurABCBoundary <: BoundaryCondition end
Base.@kwdef struct PMLBoundary <: BoundaryCondition
    ncells::Int = 24
    grading_order::Int = 3
    sigma_max::Float64 = 0.0
    target_reflection::Float64 = 1.0e-6
end
struct NoBoundaryCondition <: BoundaryCondition end

abstract type CurrentSource end
abstract type AbstractExcitation end

struct SoftCurrentExcitation <: AbstractExcitation end

struct TFSFExcitation <: AbstractExcitation
    start_index::Int
    end_index::Int
    incident_eps_r::Float64
    incident_mu_r::Float64
end

struct FieldMonitor
    name::String
    index::Int
end

"""
Simple Gaussian pulse in time.
"""
struct GaussianCurrent <: CurrentSource
    amplitude::Float64
    t0::Float64
    spread::Float64
end

"""
Continuous cosine source.
"""
struct CosineCurrent <: CurrentSource
    amplitude::Float64
    frequency::Float64
    phase::Float64
end

"""
Cosine carrier inside a Gaussian envelope.
"""
struct GaussianModulatedCosineCurrent <: CurrentSource
    amplitude::Float64
    frequency::Float64
    t0::Float64
    spread::Float64
    phase::Float64
end

"""
Ricker wavelet, often useful for broadband transient simulations.
"""
struct RickerCurrent <: CurrentSource
    amplitude::Float64
    frequency::Float64
    t0::Float64
end

"""
Core simulation configuration.

The boundary fields are explicit types so additional boundary treatments
such as higher-order ABCs or more advanced PML variants can be added later via new methods.
"""
Base.@kwdef struct SimulationConfig
    nx::Int = 401
    dx::Float64 = 1e-3
    courant_factor::Float64 = 0.99
    nsteps::Int = 1_200
    source_position::Union{Int,Nothing} = nothing
    source::CurrentSource = GaussianCurrent(1.0, 40e-12, 12e-12)
    excitation::AbstractExcitation = SoftCurrentExcitation()
    left_boundary::BoundaryCondition = PECBoundary()
    right_boundary::BoundaryCondition = PECBoundary()
    material::AbstractMaterial = Vacuum()
    monitors::Vector{FieldMonitor} = FieldMonitor[]
    save_every::Int = 2
end

struct SimulationResult
    x::Vector{Float64}
    times::Vector{Float64}
    e_history::Matrix{Float64}
    h_history::Matrix{Float64}
    monitor_traces::Dict{String,Vector{Float64}}
    config::SimulationConfig
end

material_permittivity(::AbstractMaterial) = eps0
material_permeability(::AbstractMaterial) = mu0

function electric_permittivity_profile(::Vacuum, nx::Int)
    return fill(eps0, nx)
end

function magnetic_permeability_profile(::Vacuum, nx::Int)
    return fill(mu0, nx - 1)
end

function electric_conductivity_profile(::Vacuum, nx::Int)
    return zeros(nx)
end

function electric_permittivity_profile(material::GridMaterial, nx::Int)
    length(material.eps_r_e) == nx ||
        throw(ArgumentError("GridMaterial eps_r_e length must match nx."))
    return eps0 .* material.eps_r_e
end

function magnetic_permeability_profile(material::GridMaterial, nx::Int)
    length(material.mu_r_h) == nx - 1 ||
        throw(ArgumentError("GridMaterial mu_r_h length must be nx - 1."))
    return mu0 .* material.mu_r_h
end

function electric_conductivity_profile(material::GridMaterial, nx::Int)
    length(material.sigma_e) == nx ||
        throw(ArgumentError("GridMaterial sigma_e length must match nx."))
    return material.sigma_e
end

function electric_permittivity_profile(material::DebyeMaterial, nx::Int)
    length(material.eps_inf_e) == nx ||
        throw(ArgumentError("DebyeMaterial eps_inf_e length must match nx."))
    return eps0 .* material.eps_inf_e
end

function magnetic_permeability_profile(material::DebyeMaterial, nx::Int)
    length(material.mu_r_h) == nx - 1 ||
        throw(ArgumentError("DebyeMaterial mu_r_h length must be nx - 1."))
    return mu0 .* material.mu_r_h
end

function electric_conductivity_profile(material::DebyeMaterial, nx::Int)
    length(material.sigma_e) == nx ||
        throw(ArgumentError("DebyeMaterial sigma_e length must match nx."))
    return material.sigma_e
end

function electric_permittivity_profile(material::DrudeMaterial, nx::Int)
    length(material.eps_inf_e) == nx ||
        throw(ArgumentError("DrudeMaterial eps_inf_e length must match nx."))
    return eps0 .* material.eps_inf_e
end

function magnetic_permeability_profile(material::DrudeMaterial, nx::Int)
    length(material.mu_r_h) == nx - 1 ||
        throw(ArgumentError("DrudeMaterial mu_r_h length must be nx - 1."))
    return mu0 .* material.mu_r_h
end

function electric_conductivity_profile(material::DrudeMaterial, nx::Int)
    length(material.sigma_e) == nx ||
        throw(ArgumentError("DrudeMaterial sigma_e length must match nx."))
    return material.sigma_e
end

function electric_permittivity_profile(material::LorentzMaterial, nx::Int)
    length(material.eps_inf_e) == nx ||
        throw(ArgumentError("LorentzMaterial eps_inf_e length must match nx."))
    return eps0 .* material.eps_inf_e
end

function magnetic_permeability_profile(material::LorentzMaterial, nx::Int)
    length(material.mu_r_h) == nx - 1 ||
        throw(ArgumentError("LorentzMaterial mu_r_h length must be nx - 1."))
    return mu0 .* material.mu_r_h
end

function electric_conductivity_profile(material::LorentzMaterial, nx::Int)
    length(material.sigma_e) == nx ||
        throw(ArgumentError("LorentzMaterial sigma_e length must match nx."))
    return material.sigma_e
end

function material_regions(::Vacuum, x, dx)
    return NamedTuple[]
end

function material_regions(material::GridMaterial, x, dx)
    regions = NamedTuple[]
    start_idx = nothing

    for i in eachindex(material.eps_r_e)
        is_background = isapprox(material.eps_r_e[i], 1.0; atol = 1e-12) &&
                        isapprox(material.sigma_e[i], 0.0; atol = 1e-12)

        if !is_background && isnothing(start_idx)
            start_idx = i
        elseif is_background && !isnothing(start_idx)
            stop_idx = i - 1
            push!(regions, (
                x_center = (x[start_idx] + x[stop_idx]) / 2,
                x0 = x[start_idx] - dx / 2,
                x1 = x[stop_idx] + dx / 2,
                eps_r = material.eps_r_e[start_idx],
                sigma_e = material.sigma_e[start_idx],
            ))
            start_idx = nothing
        end
    end

    if !isnothing(start_idx)
        push!(regions, (
            x_center = (x[start_idx] + x[end]) / 2,
            x0 = x[start_idx] - dx / 2,
            x1 = x[end] + dx / 2,
            eps_r = material.eps_r_e[start_idx],
            sigma_e = material.sigma_e[start_idx],
        ))
    end

    return regions
end

function material_regions(material::DebyeMaterial, x, dx)
    regions = NamedTuple[]
    start_idx = nothing

    for i in eachindex(material.eps_inf_e)
        is_background = isapprox(material.eps_inf_e[i], 1.0; atol = 1e-12) &&
                        isapprox(material.eps_static_e[i], 1.0; atol = 1e-12) &&
                        isapprox(material.sigma_e[i], 0.0; atol = 1e-12)

        if !is_background && isnothing(start_idx)
            start_idx = i
        elseif is_background && !isnothing(start_idx)
            stop_idx = i - 1
            push!(regions, (
                x_center = (x[start_idx] + x[stop_idx]) / 2,
                x0 = x[start_idx] - dx / 2,
                x1 = x[stop_idx] + dx / 2,
                eps_r = material.eps_static_e[start_idx],
                eps_inf = material.eps_inf_e[start_idx],
                tau = material.tau_e[start_idx],
                sigma_e = material.sigma_e[start_idx],
            ))
            start_idx = nothing
        end
    end

    if !isnothing(start_idx)
        push!(regions, (
            x_center = (x[start_idx] + x[end]) / 2,
            x0 = x[start_idx] - dx / 2,
            x1 = x[end] + dx / 2,
            eps_r = material.eps_static_e[start_idx],
            eps_inf = material.eps_inf_e[start_idx],
            tau = material.tau_e[start_idx],
            sigma_e = material.sigma_e[start_idx],
        ))
    end

    return regions
end

function material_regions(material::DrudeMaterial, x, dx)
    regions = NamedTuple[]
    start_idx = nothing

    for i in eachindex(material.eps_inf_e)
        is_background = isapprox(material.eps_inf_e[i], 1.0; atol = 1e-12) &&
                        isapprox(material.omega_p_e[i], 0.0; atol = 1e-12) &&
                        isapprox(material.sigma_e[i], 0.0; atol = 1e-12)

        if !is_background && isnothing(start_idx)
            start_idx = i
        elseif is_background && !isnothing(start_idx)
            stop_idx = i - 1
            push!(regions, (
                x_center = (x[start_idx] + x[stop_idx]) / 2,
                x0 = x[start_idx] - dx / 2,
                x1 = x[stop_idx] + dx / 2,
                eps_r = material.eps_inf_e[start_idx],
                omega_p = material.omega_p_e[start_idx],
                gamma = material.gamma_e[start_idx],
                sigma_e = material.sigma_e[start_idx],
            ))
            start_idx = nothing
        end
    end

    if !isnothing(start_idx)
        push!(regions, (
            x_center = (x[start_idx] + x[end]) / 2,
            x0 = x[start_idx] - dx / 2,
            x1 = x[end] + dx / 2,
            eps_r = material.eps_inf_e[start_idx],
            omega_p = material.omega_p_e[start_idx],
            gamma = material.gamma_e[start_idx],
            sigma_e = material.sigma_e[start_idx],
        ))
    end

    return regions
end

function material_regions(material::LorentzMaterial, x, dx)
    regions = NamedTuple[]
    start_idx = nothing

    for i in eachindex(material.eps_inf_e)
        is_background = isapprox(material.eps_inf_e[i], 1.0; atol = 1e-12) &&
                        isapprox(material.delta_eps_e[i], 0.0; atol = 1e-12) &&
                        isapprox(material.sigma_e[i], 0.0; atol = 1e-12)

        if !is_background && isnothing(start_idx)
            start_idx = i
        elseif is_background && !isnothing(start_idx)
            stop_idx = i - 1
            push!(regions, (
                x_center = (x[start_idx] + x[stop_idx]) / 2,
                x0 = x[start_idx] - dx / 2,
                x1 = x[stop_idx] + dx / 2,
                eps_r = material.eps_inf_e[start_idx] + material.delta_eps_e[start_idx],
                eps_inf = material.eps_inf_e[start_idx],
                delta_eps = material.delta_eps_e[start_idx],
                omega_0 = material.omega_0_e[start_idx],
                gamma = material.gamma_e[start_idx],
                sigma_e = material.sigma_e[start_idx],
            ))
            start_idx = nothing
        end
    end

    if !isnothing(start_idx)
        push!(regions, (
            x_center = (x[start_idx] + x[end]) / 2,
            x0 = x[start_idx] - dx / 2,
            x1 = x[end] + dx / 2,
            eps_r = material.eps_inf_e[start_idx] + material.delta_eps_e[start_idx],
            eps_inf = material.eps_inf_e[start_idx],
            delta_eps = material.delta_eps_e[start_idx],
            omega_0 = material.omega_0_e[start_idx],
            gamma = material.gamma_e[start_idx],
            sigma_e = material.sigma_e[start_idx],
        ))
    end

    return regions
end

function region_color(region, eps_min, eps_max)
    loss_weight = region.sigma_e > 0 ? 0.35 : 0.0
    if isapprox(eps_max, eps_min; atol = 1e-12)
        t = 0.5
    else
        t = clamp((region.eps_r - eps_min) / (eps_max - eps_min), 0.0, 1.0)
    end

    # Blend from cool blue to warm orange and shift toward red for lossy media.
    r = clamp(0.15 + 0.75 * t + 0.25 * loss_weight, 0.0, 1.0)
    g = clamp(0.45 + 0.25 * (1 - t) - 0.25 * loss_weight, 0.0, 1.0)
    b = clamp(0.9 - 0.55 * t - 0.45 * loss_weight, 0.0, 1.0)
    return Plots.RGB(r, g, b)
end

function region_label(region)
    label = "eps_r=$(round(region.eps_r; digits = 2))"
    if :eps_inf in keys(region)
        label *= ", eps_inf=$(round(region.eps_inf; digits = 2))"
    end
    if :tau in keys(region)
        label *= ", tau=$(round(region.tau * 1e12; digits = 3)) ps"
    end
    if :delta_eps in keys(region)
        label *= ", delta_eps=$(round(region.delta_eps; digits = 2))"
    end
    if :omega_p in keys(region)
        label *= ", wp=$(round(region.omega_p / (2π * 1e9); digits = 3)) GHz"
    end
    if :omega_0 in keys(region)
        label *= ", w0=$(round(region.omega_0 / (2π * 1e9); digits = 3)) GHz"
    end
    if :gamma in keys(region)
        label *= ", gamma=$(round(region.gamma / (2π * 1e9); digits = 3)) GHz"
    end
    if region.sigma_e > 0
        label *= ", sigma=$(round(region.sigma_e; digits = 4))"
    end
    return label
end

boundary_label(::PECBoundary) = "PEC"
boundary_label(::MurABCBoundary) = "ABC"
boundary_label(boundary::PMLBoundary) = "PML($(boundary.ncells))"
boundary_label(::NoBoundaryCondition) = "Open"
excitation_label(::SoftCurrentExcitation) = "Soft current"
excitation_label(excitation::TFSFExcitation) =
    "TF/SF $(excitation.start_index):$(excitation.end_index) eps=$(round(excitation.incident_eps_r; digits = 2)) mu=$(round(excitation.incident_mu_r; digits = 2))"

function excitation_markers(config::SimulationConfig, x)
    if config.excitation isa SoftCurrentExcitation
        return [x[resolve_source_position(config)]]
    elseif config.excitation isa TFSFExcitation
        excitation = config.excitation
        return [x[excitation.start_index], x[excitation.end_index + 1]]
    end
    return Float64[]
end

function resolve_grid_index(
    data::AbstractDict,
    index_key::AbstractString,
    x_key::AbstractString,
    dx::Real,
    nx::Int;
    default = nothing,
)
    has_index = haskey(data, index_key)
    has_x = haskey(data, x_key)

    if has_index && has_x
        throw(ArgumentError("Specify only one of '$index_key' or '$x_key'."))
    elseif has_index
        return Int(data[index_key])
    elseif has_x
        x = Float64(data[x_key])
        x < 0 && throw(ArgumentError("'$x_key' must be nonnegative."))
        return round(Int, x / dx) + 1
    end

    return default
end

function resolve_required_grid_index(
    data::AbstractDict,
    index_key::AbstractString,
    x_key::AbstractString,
    dx::Real,
    nx::Int,
)
    idx = resolve_grid_index(data, index_key, x_key, dx, nx; default = nothing)
    isnothing(idx) && throw(ArgumentError("One of '$index_key' or '$x_key' must be provided."))
    return idx
end

function monitor_from_dict(data::AbstractDict, dx::Real, nx::Int)
    return FieldMonitor(
        String(data["name"]),
        resolve_required_grid_index(data, "index", "x", dx, nx),
    )
end

function resolve_layer_bounds(data::AbstractDict, dx::Real, nx::Int)
    has_start = haskey(data, "start_index") || haskey(data, "start_x")
    has_end = haskey(data, "end_index") || haskey(data, "end_x")
    has_center = haskey(data, "center_x")
    has_thickness = haskey(data, "thickness")

    if has_start || has_end
        (has_start && has_end) ||
            throw(ArgumentError("Layer bounds require both start and end definitions."))
        !(has_center || has_thickness) ||
            throw(ArgumentError("Use either start/end or center_x/thickness for a layer, not both."))

        return (
            resolve_required_grid_index(data, "start_index", "start_x", dx, nx),
            resolve_required_grid_index(data, "end_index", "end_x", dx, nx),
        )
    elseif has_center || has_thickness
        (has_center && has_thickness) ||
            throw(ArgumentError("Layer center/thickness definition requires both 'center_x' and 'thickness'."))

        center_x = Float64(data["center_x"])
        thickness = Float64(data["thickness"])
        center_x >= 0 || throw(ArgumentError("'center_x' must be nonnegative."))
        thickness > 0 || throw(ArgumentError("'thickness' must be positive."))

        start_x = center_x - thickness / 2
        end_x = center_x + thickness / 2
        start_x >= 0 || throw(ArgumentError("Layer thickness extends below x = 0."))

        return (
            resolve_required_grid_index(Dict("start_x" => start_x), "start_index", "start_x", dx, nx),
            resolve_required_grid_index(Dict("end_x" => end_x), "end_index", "end_x", dx, nx),
        )
    end

    throw(ArgumentError("Each material layer must define either start/end bounds or center_x/thickness."))
end

function make_layered_material(nx::Int, dx::Real, layers_data)
    eps_r_e = ones(nx)
    mu_r_h = ones(nx - 1)
    sigma_e = zeros(nx)

    for entry in layers_data
        start_index, end_index = resolve_layer_bounds(entry, dx, nx)
        layer = MaterialLayer(
            start_index,
            end_index,
            Float64(get(entry, "eps_r", 1.0)),
            Float64(get(entry, "mu_r", 1.0)),
            Float64(get(entry, "sigma_e", 0.0)),
        )

        1 <= layer.start_index <= layer.end_index <= nx ||
            throw(ArgumentError("Material layer indices must satisfy 1 <= start_index <= end_index <= nx."))
        layer.eps_r > 0 || throw(ArgumentError("Material eps_r must be positive."))
        layer.mu_r > 0 || throw(ArgumentError("Material mu_r must be positive."))
        layer.sigma_e >= 0 || throw(ArgumentError("Material sigma_e must be nonnegative."))

        eps_r_e[layer.start_index:layer.end_index] .= layer.eps_r
        sigma_e[layer.start_index:layer.end_index] .= layer.sigma_e

        h_start = max(1, layer.start_index)
        h_end = min(nx - 1, layer.end_index - 1)
        if h_start <= h_end
            mu_r_h[h_start:h_end] .= layer.mu_r
        end
    end

    return GridMaterial(eps_r_e, mu_r_h, sigma_e)
end

function make_debye_material(nx::Int, dx::Real, layers_data)
    eps_inf_e = ones(nx)
    eps_static_e = ones(nx)
    tau_e = ones(nx)
    mu_r_h = ones(nx - 1)
    sigma_e = zeros(nx)

    for entry in layers_data
        start_index, end_index = resolve_layer_bounds(entry, dx, nx)
        eps_inf = Float64(get(entry, "eps_inf", 1.0))
        eps_static = Float64(get(entry, "eps_static", eps_inf))
        tau = Float64(get(entry, "tau", 1.0))
        mu_r = Float64(get(entry, "mu_r", 1.0))
        sigma = Float64(get(entry, "sigma_e", 0.0))

        1 <= start_index <= end_index <= nx ||
            throw(ArgumentError("Debye material layer indices must satisfy 1 <= start_index <= end_index <= nx."))
        eps_inf > 0 || throw(ArgumentError("Debye eps_inf must be positive."))
        eps_static >= eps_inf || throw(ArgumentError("Debye eps_static must be >= eps_inf."))
        tau > 0 || throw(ArgumentError("Debye tau must be positive."))
        mu_r > 0 || throw(ArgumentError("Debye mu_r must be positive."))
        sigma >= 0 || throw(ArgumentError("Debye sigma_e must be nonnegative."))

        eps_inf_e[start_index:end_index] .= eps_inf
        eps_static_e[start_index:end_index] .= eps_static
        tau_e[start_index:end_index] .= tau
        sigma_e[start_index:end_index] .= sigma

        h_start = max(1, start_index)
        h_end = min(nx - 1, end_index - 1)
        if h_start <= h_end
            mu_r_h[h_start:h_end] .= mu_r
        end
    end

    return DebyeMaterial(eps_inf_e, eps_static_e, tau_e, mu_r_h, sigma_e)
end

function make_drude_material(nx::Int, dx::Real, layers_data)
    eps_inf_e = ones(nx)
    omega_p_e = zeros(nx)
    gamma_e = zeros(nx)
    mu_r_h = ones(nx - 1)
    sigma_e = zeros(nx)

    for entry in layers_data
        start_index, end_index = resolve_layer_bounds(entry, dx, nx)
        eps_inf = Float64(get(entry, "eps_inf", 1.0))
        omega_p = Float64(get(entry, "omega_p", 0.0))
        gamma = Float64(get(entry, "gamma", 0.0))
        mu_r = Float64(get(entry, "mu_r", 1.0))
        sigma = Float64(get(entry, "sigma_e", 0.0))

        1 <= start_index <= end_index <= nx ||
            throw(ArgumentError("Drude material layer indices must satisfy 1 <= start_index <= end_index <= nx."))
        eps_inf > 0 || throw(ArgumentError("Drude eps_inf must be positive."))
        omega_p >= 0 || throw(ArgumentError("Drude omega_p must be nonnegative."))
        gamma >= 0 || throw(ArgumentError("Drude gamma must be nonnegative."))
        mu_r > 0 || throw(ArgumentError("Drude mu_r must be positive."))
        sigma >= 0 || throw(ArgumentError("Drude sigma_e must be nonnegative."))

        eps_inf_e[start_index:end_index] .= eps_inf
        omega_p_e[start_index:end_index] .= omega_p
        gamma_e[start_index:end_index] .= gamma
        sigma_e[start_index:end_index] .= sigma

        h_start = max(1, start_index)
        h_end = min(nx - 1, end_index - 1)
        if h_start <= h_end
            mu_r_h[h_start:h_end] .= mu_r
        end
    end

    return DrudeMaterial(eps_inf_e, omega_p_e, gamma_e, mu_r_h, sigma_e)
end

function make_lorentz_material(nx::Int, dx::Real, layers_data)
    eps_inf_e = ones(nx)
    delta_eps_e = zeros(nx)
    omega_0_e = zeros(nx)
    gamma_e = zeros(nx)
    mu_r_h = ones(nx - 1)
    sigma_e = zeros(nx)

    for entry in layers_data
        start_index, end_index = resolve_layer_bounds(entry, dx, nx)
        eps_inf = Float64(get(entry, "eps_inf", 1.0))
        delta_eps = Float64(get(entry, "delta_eps", 0.0))
        omega_0 = Float64(get(entry, "omega_0", 0.0))
        gamma = Float64(get(entry, "gamma", 0.0))
        mu_r = Float64(get(entry, "mu_r", 1.0))
        sigma = Float64(get(entry, "sigma_e", 0.0))

        1 <= start_index <= end_index <= nx ||
            throw(ArgumentError("Lorentz material layer indices must satisfy 1 <= start_index <= end_index <= nx."))
        eps_inf > 0 || throw(ArgumentError("Lorentz eps_inf must be positive."))
        delta_eps >= 0 || throw(ArgumentError("Lorentz delta_eps must be nonnegative."))
        omega_0 > 0 || throw(ArgumentError("Lorentz omega_0 must be positive."))
        gamma >= 0 || throw(ArgumentError("Lorentz gamma must be nonnegative."))
        mu_r > 0 || throw(ArgumentError("Lorentz mu_r must be positive."))
        sigma >= 0 || throw(ArgumentError("Lorentz sigma_e must be nonnegative."))

        eps_inf_e[start_index:end_index] .= eps_inf
        delta_eps_e[start_index:end_index] .= delta_eps
        omega_0_e[start_index:end_index] .= omega_0
        gamma_e[start_index:end_index] .= gamma
        sigma_e[start_index:end_index] .= sigma

        h_start = max(1, start_index)
        h_end = min(nx - 1, end_index - 1)
        if h_start <= h_end
            mu_r_h[h_start:h_end] .= mu_r
        end
    end

    return LorentzMaterial(eps_inf_e, delta_eps_e, omega_0_e, gamma_e, mu_r_h, sigma_e)
end

function material_from_dict(data::AbstractDict, nx::Int, dx::Real)
    material_name = lowercase(String(get(data, "type", "vacuum")))

    if material_name == "vacuum"
        return Vacuum()
    elseif material_name in ("layered", "piecewise_constant")
        layers = get(data, "layers", Any[])
        return make_layered_material(nx, dx, layers)
    elseif material_name in ("debye", "dispersive_debye")
        layers = get(data, "layers", Any[])
        return make_debye_material(nx, dx, layers)
    elseif material_name in ("drude", "dispersive_drude")
        layers = get(data, "layers", Any[])
        return make_drude_material(nx, dx, layers)
    elseif material_name in ("lorentz", "dispersive_lorentz")
        layers = get(data, "layers", Any[])
        return make_lorentz_material(nx, dx, layers)
    end

    throw(ArgumentError("Unknown material type '$material_name'. Supported values: vacuum, layered, debye, drude, lorentz."))
end

function available_source_templates()
    return Dict(
        :gaussian => GaussianCurrent(1.0, 40e-12, 12e-12),
        :cosine => CosineCurrent(0.2, 1.5e9, 0.0),
        :gaussian_modulated_cosine => GaussianModulatedCosineCurrent(1.0, 1.5e9, 60e-12, 18e-12, 0.0),
        :ricker => RickerCurrent(1.0, 1.0e9, 50e-12),
    )
end

function boundary_from_name(name::AbstractString)
    normalized = lowercase(strip(name))
    if normalized == "pec"
        return PECBoundary()
    elseif normalized in ("abc", "mur", "mur1", "mur_1")
        return MurABCBoundary()
    elseif normalized == "pml"
        return PMLBoundary()
    elseif normalized in ("none", "open")
        return NoBoundaryCondition()
    end
    throw(ArgumentError("Unknown boundary '$name'. Supported values: pec, abc, pml, none."))
end

function boundary_parameter(data::AbstractDict, side::AbstractString, key::AbstractString, default)
    side_key = "$(side)_$(key)"
    if haskey(data, side_key)
        return data[side_key]
    elseif haskey(data, key)
        return data[key]
    end
    return default
end

function boundary_region_cells_from_dict(
    data::AbstractDict,
    side::AbstractString,
    region_name::AbstractString,
    dx::Real,
    default::Int;
    cells_suffix::AbstractString = "cells",
    length_suffix::AbstractString = "thickness",
)
    global_cells_key = "$(region_name)_$(cells_suffix)"
    global_length_key = "$(region_name)_$(length_suffix)"
    cells_key = "$(side)_$(region_name)_$(cells_suffix)"
    length_key = "$(side)_$(region_name)_$(length_suffix)"
    has_side_cells = haskey(data, cells_key)
    has_side_length = haskey(data, length_key)
    has_global_cells = haskey(data, global_cells_key)
    has_global_length = haskey(data, global_length_key)

    if has_side_cells && has_side_length
        throw(ArgumentError("Specify only one of '$cells_key' or '$length_key'."))
    elseif has_side_cells
        return Int(data[cells_key])
    elseif has_side_length
        length_value = Float64(data[length_key])
        length_value > 0 || throw(ArgumentError("'$length_key' must be positive."))
        return max(1, round(Int, length_value / dx))
    end

    if has_global_cells && has_global_length
        throw(ArgumentError("Specify only one of '$global_cells_key' or '$global_length_key'."))
    elseif has_global_cells
        return Int(data[global_cells_key])
    elseif has_global_length
        length_value = Float64(data[global_length_key])
        length_value > 0 || throw(ArgumentError("'$global_length_key' must be positive."))
        return max(1, round(Int, length_value / dx))
    end

    return default
end

function boundary_from_dict(data::AbstractDict, side::AbstractString, dx::Real)
    boundary = boundary_from_name(String(get(data, side, "pec")))
    if boundary isa PMLBoundary
        return PMLBoundary(
            ncells = boundary_region_cells_from_dict(data, side, "pml", dx, boundary.ncells),
            grading_order = Int(boundary_parameter(data, side, "pml_order", boundary.grading_order)),
            sigma_max = Float64(boundary_parameter(data, side, "pml_sigma_max", boundary.sigma_max)),
            target_reflection = Float64(boundary_parameter(data, side, "pml_target_reflection", boundary.target_reflection)),
        )
    end
    return boundary
end

function source_from_dict(data::AbstractDict)
    kind = Symbol(lowercase(String(get(data, "type", "gaussian"))))

    if kind === :gaussian
        return GaussianCurrent(
            Float64(get(data, "amplitude", 1.0)),
            Float64(get(data, "t0", 40e-12)),
            Float64(get(data, "spread", 12e-12)),
        )
    elseif kind === :cosine
        return CosineCurrent(
            Float64(get(data, "amplitude", 0.2)),
            Float64(get(data, "frequency", 1.5e9)),
            Float64(get(data, "phase", 0.0)),
        )
    elseif kind in (:gaussian_modulated_cosine, :modulated_cosine)
        return GaussianModulatedCosineCurrent(
            Float64(get(data, "amplitude", 1.0)),
            Float64(get(data, "frequency", 1.5e9)),
            Float64(get(data, "t0", 60e-12)),
            Float64(get(data, "spread", 18e-12)),
            Float64(get(data, "phase", 0.0)),
        )
    elseif kind === :ricker
        return RickerCurrent(
            Float64(get(data, "amplitude", 1.0)),
            Float64(get(data, "frequency", 1.0e9)),
            Float64(get(data, "t0", 50e-12)),
        )
    end

    throw(ArgumentError("Unknown source type '$kind'."))
end

function excitation_from_dict(data::AbstractDict, nx::Int, dx::Real)
    mode = lowercase(String(get(data, "type", "soft")))

    if mode in ("soft", "current", "soft_current")
        return SoftCurrentExcitation()
    elseif mode in ("tfsf", "tfsf_plane_wave")
        start_index = resolve_grid_index(data, "start_index", "start_x", dx, nx; default = fld(nx, 4))
        end_index = resolve_grid_index(data, "end_index", "end_x", dx, nx; default = nx - fld(nx, 4))
        incident_eps_r = Float64(get(data, "incident_eps_r", 1.0))
        incident_mu_r = Float64(get(data, "incident_mu_r", 1.0))
        return TFSFExcitation(start_index, end_index, incident_eps_r, incident_mu_r)
    end

    throw(ArgumentError("Unknown excitation type '$mode'. Supported values: soft, tfsf."))
end

function load_simulation_parameters(path::AbstractString)
    data = TOML.parsefile(path)
    sim = get(data, "simulation", Dict{String,Any}())
    source_data = get(data, "source", Dict{String,Any}())
    excitation_data = get(data, "excitation", Dict{String,Any}())
    boundary_data = get(data, "boundary", Dict{String,Any}())
    material_data = get(data, "material", Dict{String,Any}())
    monitor_data = get(data, "monitors", Any[])

    nx = Int(get(sim, "nx", 401))
    dx = Float64(get(sim, "dx", 1e-3))
    source_position = resolve_grid_index(sim, "source_position", "source_x", dx, nx; default = nothing)

    return (
        config = SimulationConfig(
            nx = nx,
            dx = dx,
            courant_factor = Float64(get(sim, "courant_factor", 0.99)),
            nsteps = Int(get(sim, "nsteps", 1_200)),
            source_position = source_position,
            save_every = Int(get(sim, "save_every", 2)),
            source = source_from_dict(source_data),
            excitation = excitation_from_dict(excitation_data, nx, dx),
            left_boundary = boundary_from_dict(boundary_data, "left", dx),
            right_boundary = boundary_from_dict(boundary_data, "right", dx),
            material = material_from_dict(material_data, nx, dx),
            monitors = [monitor_from_dict(entry, dx, nx) for entry in monitor_data],
        ),
        output = get(data, "output", Dict{String,Any}()),
        diagnostics = get(data, "diagnostics", Dict{String,Any}()),
    )
end

struct OutputManager
    root_dir::String
    run_dir::String
    case_name::String
    config_path::Union{Nothing,String}
end

function sanitize_case_name(name::AbstractString)
    sanitized = replace(lowercase(strip(name)), r"[^a-z0-9_.-]+" => "_")
    sanitized = strip(sanitized, '_')
    return isempty(sanitized) ? "fdtd_run" : sanitized
end

timestamp_for_path() = Libc.strftime("%Y%m%d_%H%M%S", time())
timestamp_iso() = Libc.strftime("%Y-%m-%dT%H:%M:%S", time())

function create_output_manager(
    config_path::Union{Nothing,AbstractString} = nothing;
    root::AbstractString = "outputs",
    case_name::Union{Nothing,AbstractString} = nothing,
    timestamped::Bool = true,
)
    resolved_case = isnothing(case_name) || isempty(strip(case_name)) ?
        (isnothing(config_path) ? "fdtd_run" : splitext(basename(String(config_path)))[1]) :
        String(case_name)
    safe_case = sanitize_case_name(resolved_case)
    stamp = timestamp_for_path()
    run_name = timestamped ? "$(safe_case)_$(stamp)" : safe_case
    run_dir = joinpath(String(root), run_name)
    mkpath(run_dir)
    return OutputManager(String(root), run_dir, safe_case, isnothing(config_path) ? nothing : String(config_path))
end

function output_path(manager::OutputManager, requested::AbstractString; default_name::AbstractString = "output.dat")
    name = isempty(strip(requested)) ? default_name : basename(String(requested))
    return joinpath(manager.run_dir, name)
end

function copy_input_config!(manager::OutputManager; filename::AbstractString = "config.toml")
    isnothing(manager.config_path) && return nothing
    destination = joinpath(manager.run_dir, filename)
    cp(manager.config_path, destination; force = true)
    return destination
end

function write_csv_row(io, values)
    println(io, join(values, ","))
end

function write_monitor_traces_csv(path::AbstractString, result::SimulationResult)
    mkpath(dirname(abspath(path)))
    monitor_names = sort(collect(keys(result.monitor_traces)))
    open(path, "w") do io
        write_csv_row(io, ["time_s"; monitor_names])
        for i in eachindex(result.times)
            write_csv_row(io, [result.times[i]; [result.monitor_traces[name][i] for name in monitor_names]])
        end
    end
    return path
end

function write_material_profile_csv(path::AbstractString, config::SimulationConfig)
    mkpath(dirname(abspath(path)))
    ϵ = electric_permittivity_profile(config.material, config.nx) ./ eps0
    σe = electric_conductivity_profile(config.material, config.nx)
    μ = magnetic_permeability_profile(config.material, config.nx) ./ mu0
    open(path, "w") do io
        println(io, "index,x_m,eps_r_e,sigma_e,mu_r_left_h,mu_r_right_h")
        for i in 1:config.nx
            left_mu = i > 1 ? μ[i - 1] : ""
            right_mu = i <= length(μ) ? μ[i] : ""
            write_csv_row(io, (i, (i - 1) * config.dx, ϵ[i], σe[i], left_mu, right_mu))
        end
    end
    return path
end

function write_spectrum_csv(path::AbstractString, spectrum)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "frequency_hz,reflection,transmission,incident_real,incident_imag,reflected_real,reflected_imag,transmitted_real,transmitted_imag")
        for i in eachindex(spectrum.frequencies)
            write_csv_row(io, (
                spectrum.frequencies[i],
                spectrum.reflection[i],
                spectrum.transmission[i],
                real(spectrum.incident_spectrum[i]),
                imag(spectrum.incident_spectrum[i]),
                real(spectrum.reflected_spectrum[i]),
                imag(spectrum.reflected_spectrum[i]),
                real(spectrum.transmitted_spectrum[i]),
                imag(spectrum.transmitted_spectrum[i]),
            ))
        end
    end
    return path
end

function write_sparameters_csv(path::AbstractString, sparams)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "frequency_hz,s11_real,s11_imag,s21_real,s21_imag,s12_real,s12_imag,s22_real,s22_imag")
        for i in eachindex(sparams.frequencies)
            write_csv_row(io, (
                sparams.frequencies[i],
                real(sparams.s11[i]), imag(sparams.s11[i]),
                real(sparams.s21[i]), imag(sparams.s21[i]),
                real(sparams.s12[i]), imag(sparams.s12[i]),
                real(sparams.s22[i]), imag(sparams.s22[i]),
            ))
        end
    end
    return path
end

function toml_sanitize(value)
    if value isa AbstractDict
        return Dict{String,Any}(string(k) => toml_sanitize(v) for (k, v) in value if !isnothing(v))
    elseif value isa NamedTuple
        return Dict{String,Any}(string(k) => toml_sanitize(getproperty(value, k)) for k in keys(value) if !isnothing(getproperty(value, k)))
    elseif value isa Tuple
        return [toml_sanitize(v) for v in value if !isnothing(v)]
    elseif value isa AbstractVector
        return [toml_sanitize(v) for v in value if !isnothing(v)]
    elseif value isa Complex
        return string(value)
    elseif isnothing(value)
        return nothing
    end
    return value
end

function write_run_summary(
    manager::OutputManager,
    result::SimulationResult;
    diagnostics = nothing,
    spectrum = nothing,
    sparams = nothing,
    files::AbstractDict = Dict{String,Any}(),
)
    summary = Dict{String,Any}(
        "case_name" => manager.case_name,
        "run_dir" => manager.run_dir,
        "timestamp" => timestamp_iso(),
        "grid" => Dict{String,Any}(
            "nx" => result.config.nx,
            "dx" => result.config.dx,
            "nsteps" => result.config.nsteps,
            "save_every" => result.config.save_every,
            "courant_factor" => result.config.courant_factor,
        ),
        "types" => Dict{String,Any}(
            "source" => string(typeof(result.config.source)),
            "excitation" => string(typeof(result.config.excitation)),
            "material" => string(typeof(result.config.material)),
            "left_boundary" => string(typeof(result.config.left_boundary)),
            "right_boundary" => string(typeof(result.config.right_boundary)),
        ),
        "files" => Dict{String,Any}(files),
    )
    if !isnothing(diagnostics)
        summary["diagnostics"] = Dict{String,Any}(string(k) => v for (k, v) in diagnostics)
    end
    if !isnothing(spectrum)
        summary["spectrum"] = Dict{String,Any}(
            "frequency_min" => minimum(spectrum.frequencies),
            "frequency_max" => maximum(spectrum.frequencies),
            "frequency_count" => length(spectrum.frequencies),
            "window" => spectrum.window,
            "gate_start" => spectrum.gate_start,
            "gate_end" => spectrum.gate_end,
        )
    end
    if !isnothing(sparams)
        summary["sparameters"] = Dict{String,Any}(
            "frequency_min" => minimum(sparams.frequencies),
            "frequency_max" => maximum(sparams.frequencies),
            "frequency_count" => length(sparams.frequencies),
            "left_shift" => sparams.left_shift,
            "right_shift" => sparams.right_shift,
        )
        if hasproperty(sparams, :metadata)
            summary["sparameters"]["metadata"] = Dict{String,Any}(string(k) => getproperty(sparams.metadata, k) for k in keys(sparams.metadata))
        end
    end

    path = joinpath(manager.run_dir, "summary.toml")
    open(path, "w") do io
        TOML.print(io, toml_sanitize(summary))
    end
    return path
end

current_value(src::GaussianCurrent, t) =
    src.amplitude * exp(-((t - src.t0)^2) / (2 * src.spread^2))

current_value(src::CosineCurrent, t) =
    src.amplitude * cos(2 * pi * src.frequency * t + src.phase)

current_value(src::GaussianModulatedCosineCurrent, t) =
    src.amplitude *
    exp(-((t - src.t0)^2) / (2 * src.spread^2)) *
    cos(2 * pi * src.frequency * t + src.phase)

function current_value(src::RickerCurrent, t)
    a = pi * src.frequency * (t - src.t0)
    return src.amplitude * (1.0 - 2.0 * a^2) * exp(-a^2)
end

function validate(config::SimulationConfig)
    config.nx >= 3 || throw(ArgumentError("nx must be at least 3."))
    config.dx > 0 || throw(ArgumentError("dx must be positive."))
    config.nsteps >= 1 || throw(ArgumentError("nsteps must be at least 1."))
    config.save_every >= 1 || throw(ArgumentError("save_every must be at least 1."))
    config.courant_factor > 0 || throw(ArgumentError("courant_factor must be positive."))
    config.courant_factor <= 1 ||
        throw(ArgumentError("courant_factor must be <= 1 for the 1D vacuum CFL limit."))

    left_interior, right_interior = interior_domain_bounds(config)
    left_interior <= right_interior ||
        throw(ArgumentError("Boundary regions leave no interior cells. Reduce the PML thickness or increase nx."))

    for (side, boundary) in (("left", config.left_boundary), ("right", config.right_boundary))
        if boundary isa PMLBoundary
            boundary.ncells >= 1 ||
                throw(ArgumentError("$(side) PML requires pml_cells >= 1."))
            boundary.ncells <= config.nx - 2 ||
                throw(ArgumentError("$(side) PML requires pml_cells <= nx - 2."))
            boundary.grading_order >= 1 ||
                throw(ArgumentError("$(side) PML requires pml_order >= 1."))
            boundary.sigma_max >= 0 ||
                throw(ArgumentError("$(side) PML requires pml_sigma_max >= 0."))
            0 < boundary.target_reflection < 1 ||
                throw(ArgumentError("$(side) PML requires 0 < pml_target_reflection < 1."))
        end
    end

    src_idx = resolve_source_position(config)
    if config.excitation isa SoftCurrentExcitation
        left_interior <= src_idx <= right_interior ||
            throw(ArgumentError("source_position must lie outside the PML regions and inside the domain interior."))
    elseif config.excitation isa TFSFExcitation
        excitation = config.excitation
        2 <= excitation.start_index < excitation.end_index <= config.nx - 2 ||
            throw(ArgumentError("For TF/SF, require 2 <= start_index < end_index <= nx - 2."))
        left_interior <= excitation.start_index ||
            throw(ArgumentError("TF/SF start_index must lie outside the left PML region."))
        excitation.end_index + 1 <= right_interior ||
            throw(ArgumentError("TF/SF end_index must lie outside the right PML region."))
        excitation.incident_eps_r > 0 ||
            throw(ArgumentError("incident_eps_r must be positive for TF/SF excitation."))
        excitation.incident_mu_r > 0 ||
            throw(ArgumentError("incident_mu_r must be positive for TF/SF excitation."))
    end

    seen_monitors = Set{String}()
    for monitor in config.monitors
        1 <= monitor.index <= config.nx ||
            throw(ArgumentError("Monitor $(monitor.name) index must satisfy 1 <= index <= nx."))
        left_interior <= monitor.index <= right_interior ||
            throw(ArgumentError("Monitor $(monitor.name) must lie outside the PML regions."))
        monitor.name in seen_monitors &&
            throw(ArgumentError("Duplicate monitor name '$(monitor.name)'."))
        push!(seen_monitors, monitor.name)
    end

    electric_permittivity_profile(config.material, config.nx)
    magnetic_permeability_profile(config.material, config.nx)
    electric_conductivity_profile(config.material, config.nx)
    initialize_material_state(config.material, config.nx, config.courant_factor * config.dx / c0)
end

resolve_source_position(config::SimulationConfig) =
    isnothing(config.source_position) ? fld(config.nx + 1, 2) : config.source_position

pml_cells(::BoundaryCondition) = 0
pml_cells(boundary::PMLBoundary) = boundary.ncells

function interior_domain_bounds(config::SimulationConfig)
    left = pml_cells(config.left_boundary) + 1
    right = config.nx - pml_cells(config.right_boundary)
    return left, right
end

apply_boundary!(::NoBoundaryCondition, E, side, boundary_cache, cfl_profile) = E

function apply_boundary!(::PECBoundary, E, side::Symbol, boundary_cache, cfl_profile)
    if side === :left
        E[1] = 0.0
    elseif side === :right
        E[end] = 0.0
    else
        throw(ArgumentError("Unknown boundary side: $side"))
    end
    return E
end

function apply_boundary!(::MurABCBoundary, E, side::Symbol, boundary_cache, cfl_profile)
    if side === :left
        coeff = (cfl_profile.left - 1) / (cfl_profile.left + 1)
        E[1] = boundary_cache.left_inner_prev + coeff * (E[2] - boundary_cache.left_prev)
    elseif side === :right
        coeff = (cfl_profile.right - 1) / (cfl_profile.right + 1)
        E[end] = boundary_cache.right_inner_prev + coeff * (E[end - 1] - boundary_cache.right_prev)
    else
        throw(ArgumentError("Unknown boundary side: $side"))
    end
    return E
end

function apply_boundary!(::PMLBoundary, E, side::Symbol, boundary_cache, cfl_profile)
    if side === :left
        E[1] = 0.0
    elseif side === :right
        E[end] = 0.0
    else
        throw(ArgumentError("Unknown boundary side: $side"))
    end
    return E
end

Base.@kwdef mutable struct BoundaryCache
    left_prev::Float64 = 0.0
    left_inner_prev::Float64 = 0.0
    right_prev::Float64 = 0.0
    right_inner_prev::Float64 = 0.0
end

Base.@kwdef mutable struct DebyePolarizationState
    P::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
end

Base.@kwdef mutable struct DrudePolarizationState
    J::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
end

Base.@kwdef mutable struct LorentzPolarizationState
    P::Vector{Float64}
    J::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
end

function update_boundary_cache!(cache::BoundaryCache, E)
    cache.left_prev = E[1]
    cache.left_inner_prev = E[2]
    cache.right_prev = E[end]
    cache.right_inner_prev = E[end - 1]
    return cache
end

function boundary_sigma_max(boundary::PMLBoundary, dx::Real, η_edge::Real)
    if boundary.sigma_max > 0
        return boundary.sigma_max
    end
    thickness = boundary.ncells * dx
    return -((boundary.grading_order + 1) * log(boundary.target_reflection)) / (2 * η_edge * thickness)
end

graded_depth_from_outer(offset::Int, ncells::Int) =
    ncells <= 1 ? 1.0 : offset / (ncells - 1)

function apply_left_pml_profiles!(σe_pml, σm_pml, boundary::PMLBoundary, dx, ϵ, μ)
    η_edge = sqrt(μ[1] / ϵ[1])
    σmax = boundary_sigma_max(boundary, dx, η_edge)

    for i in 1:min(boundary.ncells, length(σe_pml))
        depth = graded_depth_from_outer(boundary.ncells - i, boundary.ncells)
        σe_pml[i] = max(σe_pml[i], σmax * depth^boundary.grading_order)
    end

    for i in 1:min(boundary.ncells, length(σm_pml))
        depth = graded_depth_from_outer(boundary.ncells - i, boundary.ncells)
        εh = 0.5 * (ϵ[i] + ϵ[min(i + 1, length(ϵ))])
        σm_pml[i] = max(σm_pml[i], σmax * depth^boundary.grading_order * μ[i] / εh)
    end
    return nothing
end

function apply_right_pml_profiles!(σe_pml, σm_pml, boundary::PMLBoundary, dx, ϵ, μ)
    η_edge = sqrt(μ[end] / ϵ[end])
    σmax = boundary_sigma_max(boundary, dx, η_edge)
    start_e = length(σe_pml) - boundary.ncells + 1
    start_h = length(σm_pml) - boundary.ncells + 1

    for i in max(1, start_e):length(σe_pml)
        depth = graded_depth_from_outer(i - start_e, boundary.ncells)
        σe_pml[i] = max(σe_pml[i], σmax * depth^boundary.grading_order)
    end

    for i in max(1, start_h):length(σm_pml)
        depth = graded_depth_from_outer(i - start_h, boundary.ncells)
        εh = 0.5 * (ϵ[i] + ϵ[i + 1])
        σm_pml[i] = max(σm_pml[i], σmax * depth^boundary.grading_order * μ[i] / εh)
    end
    return nothing
end

function pml_profiles(config::SimulationConfig, ϵ, μ)
    σe_pml = zeros(length(ϵ))
    σm_pml = zeros(length(μ))

    config.left_boundary isa PMLBoundary &&
        apply_left_pml_profiles!(σe_pml, σm_pml, config.left_boundary, config.dx, ϵ, μ)
    config.right_boundary isa PMLBoundary &&
        apply_right_pml_profiles!(σe_pml, σm_pml, config.right_boundary, config.dx, ϵ, μ)

    return σe_pml, σm_pml
end

initialize_material_state(::AbstractMaterial, nx::Int, dt::Real) = nothing

function initialize_material_state(material::DebyeMaterial, nx::Int, dt::Real)
    length(material.eps_static_e) == nx ||
        throw(ArgumentError("DebyeMaterial eps_static_e length must match nx."))
    length(material.tau_e) == nx ||
        throw(ArgumentError("DebyeMaterial tau_e length must match nx."))

    a = similar(material.tau_e)
    b = similar(material.tau_e)

    @inbounds for i in eachindex(material.tau_e)
        α = dt / (2 * material.tau_e[i])
        a[i] = (1 - α) / (1 + α)
        b[i] = eps0 * (material.eps_static_e[i] - material.eps_inf_e[i]) * α / (1 + α)
    end

    return DebyePolarizationState(
        P = zeros(nx),
        a = a,
        b = b,
    )
end

function initialize_material_state(material::DrudeMaterial, nx::Int, dt::Real)
    length(material.omega_p_e) == nx ||
        throw(ArgumentError("DrudeMaterial omega_p_e length must match nx."))
    length(material.gamma_e) == nx ||
        throw(ArgumentError("DrudeMaterial gamma_e length must match nx."))

    a = similar(material.omega_p_e)
    b = similar(material.omega_p_e)

    @inbounds for i in eachindex(material.omega_p_e)
        denom = 1 / dt + material.gamma_e[i] / 2
        a[i] = (1 / dt - material.gamma_e[i] / 2) / denom
        b[i] = eps0 * material.omega_p_e[i]^2 / (2 * denom)
    end

    return DrudePolarizationState(
        J = zeros(nx),
        a = a,
        b = b,
    )
end

function initialize_material_state(material::LorentzMaterial, nx::Int, dt::Real)
    length(material.delta_eps_e) == nx ||
        throw(ArgumentError("LorentzMaterial delta_eps_e length must match nx."))
    length(material.omega_0_e) == nx ||
        throw(ArgumentError("LorentzMaterial omega_0_e length must match nx."))
    length(material.gamma_e) == nx ||
        throw(ArgumentError("LorentzMaterial gamma_e length must match nx."))

    a = similar(material.delta_eps_e)
    b = similar(material.delta_eps_e)
    c = similar(material.delta_eps_e)

    @inbounds for i in eachindex(material.delta_eps_e)
        denom = 1 / dt + material.gamma_e[i] / 2 + material.omega_0_e[i]^2 * dt / 4
        a[i] = (1 / dt - material.gamma_e[i] / 2 - material.omega_0_e[i]^2 * dt / 4) / denom
        b[i] = -material.omega_0_e[i]^2 / denom
        c[i] = eps0 * material.delta_eps_e[i] * material.omega_0_e[i]^2 / (2 * denom)
    end

    return LorentzPolarizationState(
        P = zeros(nx),
        J = zeros(nx),
        a = a,
        b = b,
        c = c,
    )
end

function update_h!(H, E, dt, dx, μ, σm)
    @inbounds for i in eachindex(H)
        loss = σm[i] * dt / (2 * μ[i])
        da = (1 - loss) / (1 + loss)
        db = (dt / (μ[i] * dx)) / (1 + loss)
        H[i] = da * H[i] - db * (E[i + 1] - E[i])
    end
    return H
end

function update_e!(E, H, dt, dx, ϵ, σe)
    @inbounds for i in 2:(length(E) - 1)
        loss = σe[i] * dt / (2 * ϵ[i])
        ca = (1 - loss) / (1 + loss)
        cb = (dt / (ϵ[i] * dx)) / (1 + loss)
        E[i] = ca * E[i] - cb * (H[i] - H[i - 1])
    end
    return E
end

function update_e!(E, H, dt, dx, ϵ, σe, state::DebyePolarizationState)
    @inbounds for i in 2:(length(E) - 1)
        Eold = E[i]
        bdt = state.b[i] / dt
        loss = σe[i] / 2
        denom = ϵ[i] / dt + bdt + loss
        numer_e = ϵ[i] / dt - bdt - loss
        numer_p = (1 - state.a[i]) * state.P[i] / dt
        curl = (H[i] - H[i - 1]) / dx
        Enew = (numer_e * Eold + numer_p - curl) / denom
        E[i] = Enew
        state.P[i] = state.a[i] * state.P[i] + state.b[i] * (Enew + Eold)
    end
    return E
end

function update_e!(E, H, dt, dx, ϵ, σe, state::DrudePolarizationState)
    @inbounds for i in 2:(length(E) - 1)
        Eold = E[i]
        jscale = state.b[i] / 2
        loss = σe[i] / 2
        denom = ϵ[i] / dt + jscale + loss
        numer_e = ϵ[i] / dt - jscale - loss
        curl = (H[i] - H[i - 1]) / dx
        numer_j = -0.5 * (1 + state.a[i]) * state.J[i]
        Enew = (numer_e * Eold + numer_j - curl) / denom
        E[i] = Enew
        state.J[i] = state.a[i] * state.J[i] + state.b[i] * (Enew + Eold)
    end
    return E
end

function update_e!(E, H, dt, dx, ϵ, σe, state::LorentzPolarizationState)
    @inbounds for i in 2:(length(E) - 1)
        Eold = E[i]
        jscale = state.c[i] / 2
        loss = σe[i] / 2
        denom = ϵ[i] / dt + jscale + loss
        numer_e = ϵ[i] / dt - jscale - loss
        curl = (H[i] - H[i - 1]) / dx
        numer_aux = -0.5 * (1 + state.a[i]) * state.J[i] - 0.5 * state.b[i] * state.P[i]
        Enew = (numer_e * Eold + numer_aux - curl) / denom
        Jnew = state.a[i] * state.J[i] + state.b[i] * state.P[i] + state.c[i] * (Enew + Eold)
        E[i] = Enew
        state.P[i] += dt / 2 * (Jnew + state.J[i])
        state.J[i] = Jnew
    end
    return E
end

function inject_current_source!(E, source::CurrentSource, time, dt, ϵ, idx)
    E[idx] -= (dt / ϵ) * current_value(source, time)
    return E
end

incident_permittivity(excitation::TFSFExcitation) = eps0 * excitation.incident_eps_r
incident_permeability(excitation::TFSFExcitation) = mu0 * excitation.incident_mu_r
incident_impedance(excitation::TFSFExcitation) =
    sqrt(incident_permeability(excitation) / incident_permittivity(excitation))

incident_phase_velocity(excitation::TFSFExcitation) =
    1 / sqrt(incident_permeability(excitation) * incident_permittivity(excitation))

function initialize_monitor_traces(config::SimulationConfig, nsaved::Int)
    traces = Dict{String,Vector{Float64}}()
    for monitor in config.monitors
        traces[monitor.name] = zeros(nsaved)
    end
    return traces
end

function save_monitor_values!(traces::Dict{String,Vector{Float64}}, monitors, E, save_idx)
    for monitor in monitors
        traces[monitor.name][save_idx] = E[monitor.index]
    end
    return traces
end

function config_with_material(config::SimulationConfig, material::AbstractMaterial)
    return SimulationConfig(
        nx = config.nx,
        dx = config.dx,
        courant_factor = config.courant_factor,
        nsteps = config.nsteps,
        source_position = config.source_position,
        source = config.source,
        excitation = config.excitation,
        left_boundary = config.left_boundary,
        right_boundary = config.right_boundary,
        material = material,
        monitors = config.monitors,
        save_every = config.save_every,
    )
end

mirror_index(nx::Int, idx::Int) = nx + 1 - idx

function mirror_material(material::Vacuum)
    return material
end

function mirror_material(material::GridMaterial)
    return GridMaterial(
        reverse(material.eps_r_e),
        reverse(material.mu_r_h),
        reverse(material.sigma_e),
    )
end

function mirror_material(material::DebyeMaterial)
    return DebyeMaterial(
        reverse(material.eps_inf_e),
        reverse(material.eps_static_e),
        reverse(material.tau_e),
        reverse(material.mu_r_h),
        reverse(material.sigma_e),
    )
end

function mirror_material(material::DrudeMaterial)
    return DrudeMaterial(
        reverse(material.eps_inf_e),
        reverse(material.omega_p_e),
        reverse(material.gamma_e),
        reverse(material.mu_r_h),
        reverse(material.sigma_e),
    )
end

function mirror_material(material::LorentzMaterial)
    return LorentzMaterial(
        reverse(material.eps_inf_e),
        reverse(material.delta_eps_e),
        reverse(material.omega_0_e),
        reverse(material.gamma_e),
        reverse(material.mu_r_h),
        reverse(material.sigma_e),
    )
end

function mirror_excitation(excitation::SoftCurrentExcitation, nx::Int)
    return excitation
end

function mirror_excitation(excitation::TFSFExcitation, nx::Int)
    return TFSFExcitation(
        nx - excitation.end_index,
        nx - excitation.start_index,
        excitation.incident_eps_r,
        excitation.incident_mu_r,
    )
end

function mirror_monitors(monitors::Vector{FieldMonitor}, nx::Int)
    return [FieldMonitor(monitor.name, mirror_index(nx, monitor.index)) for monitor in monitors]
end

function mirror_config(config::SimulationConfig)
    mirrored_source_position = isnothing(config.source_position) ? nothing :
        mirror_index(config.nx, config.source_position)
    return SimulationConfig(
        nx = config.nx,
        dx = config.dx,
        courant_factor = config.courant_factor,
        nsteps = config.nsteps,
        source_position = mirrored_source_position,
        source = config.source,
        excitation = mirror_excitation(config.excitation, config.nx),
        left_boundary = config.right_boundary,
        right_boundary = config.left_boundary,
        material = mirror_material(config.material),
        monitors = mirror_monitors(config.monitors, config.nx),
        save_every = config.save_every,
    )
end

function monitor_index(config::SimulationConfig, name::AbstractString)
    for monitor in config.monitors
        if monitor.name == name
            return monitor.index
        end
    end
    throw(ArgumentError("Monitor '$name' not found in config."))
end

monitor_position(config::SimulationConfig, name::AbstractString) =
    (monitor_index(config, name) - 1) * config.dx

function incident_field_trace(config::SimulationConfig, times, monitor_name::AbstractString)
    config.excitation isa TFSFExcitation ||
        throw(ArgumentError("Excitation-derived incident fields are currently supported for TF/SF excitation."))

    excitation = config.excitation
    idx = monitor_index(config, monitor_name)
    if idx < excitation.start_index || idx > excitation.end_index + 1
        return zeros(length(times))
    end
    delay = (idx - excitation.start_index) * config.dx / incident_phase_velocity(excitation)
    return [current_value(config.source, time - delay) for time in times]
end

function compute_scattering_diagnostics(
    result::SimulationResult;
    incident_monitor::AbstractString = "incident",
    reflected_monitor::AbstractString = "incident",
    transmitted_monitor::AbstractString = "transmitted",
)
    config = result.config

    haskey(result.monitor_traces, reflected_monitor) ||
        throw(ArgumentError("Monitor '$reflected_monitor' not found in simulation result."))
    haskey(result.monitor_traces, transmitted_monitor) ||
        throw(ArgumentError("Monitor '$transmitted_monitor' not found in simulation result."))

    incident_normalization = incident_field_trace(config, result.times, incident_monitor)
    incident_reflected = incident_field_trace(config, result.times, reflected_monitor)
    incident_transmitted = incident_field_trace(config, result.times, transmitted_monitor)
    reflected_sig = result.monitor_traces[reflected_monitor] .- incident_reflected
    transmitted_sig = result.monitor_traces[transmitted_monitor]

    incident_energy = sum(abs2, incident_normalization)
    transmitted_incident_energy = sum(abs2, incident_transmitted)
    reflected_energy = sum(abs2, reflected_sig)
    transmitted_energy = sum(abs2, transmitted_sig)

    return Dict(
        "incident_energy" => incident_energy,
        "reflected_energy" => reflected_energy,
        "transmitted_energy" => transmitted_energy,
        "reflection_coefficient" => incident_energy > 0 ? reflected_energy / incident_energy : 0.0,
        "transmission_coefficient" => transmitted_incident_energy > 0 ? transmitted_energy / transmitted_incident_energy : 0.0,
    )
end

function rfft_spectrum(signal::AbstractVector{<:Real}, dt::Real)
    n = length(signal)
    spectrum = FFTW.rfft(signal) .* dt
    frequencies = collect(0:length(spectrum)-1) ./ (n * dt)
    return frequencies, spectrum
end

function interpolate_complex_spectrum(frequencies, spectrum, frequency)
    if frequency < first(frequencies) || frequency > last(frequencies)
        return 0.0 + 0.0im
    end

    idx = searchsortedlast(frequencies, frequency)
    if idx <= 0
        return spectrum[1]
    elseif idx >= length(frequencies)
        return spectrum[end]
    end

    f0 = frequencies[idx]
    f1 = frequencies[idx + 1]
    weight = (frequency - f0) / (f1 - f0)
    return (1 - weight) * spectrum[idx] + weight * spectrum[idx + 1]
end

function sample_spectrum(signal, dt, target_frequencies)
    frequencies, spectrum = rfft_spectrum(signal, dt)
    return [interpolate_complex_spectrum(frequencies, spectrum, frequency) for frequency in target_frequencies]
end

function fft_window(name::AbstractString, n::Int)
    normalized = lowercase(strip(name))
    if normalized in ("none", "rect", "rectangular")
        return ones(n)
    elseif normalized == "hann"
        return n == 1 ? ones(1) : 0.5 .- 0.5 .* cos.(2π .* (0:n-1) ./ (n - 1))
    elseif normalized == "hamming"
        return n == 1 ? ones(1) : 0.54 .- 0.46 .* cos.(2π .* (0:n-1) ./ (n - 1))
    elseif normalized == "blackman"
        if n == 1
            return ones(1)
        end
        phase = 2π .* (0:n-1) ./ (n - 1)
        return 0.42 .- 0.5 .* cos.(phase) .+ 0.08 .* cos.(2 .* phase)
    end
    throw(ArgumentError("Unknown FFT window '$name'. Supported values: none, hann, hamming, blackman."))
end

function apply_fft_window(signal::AbstractVector{<:Real}, window::AbstractString)
    weights = fft_window(window, length(signal))
    normalization = sum(weights) / length(weights)
    return collect(signal .* weights ./ normalization)
end

function gated_signal_segment(
    signal::AbstractVector{<:Real},
    times::AbstractVector{<:Real};
    gate_start::Union{Nothing,Real} = nothing,
    gate_end::Union{Nothing,Real} = nothing,
)
    length(signal) == length(times) ||
        throw(ArgumentError("Signal and time arrays must have the same length for gating."))

    in_gate(t) = (isnothing(gate_start) || t >= gate_start) &&
                 (isnothing(gate_end) || t <= gate_end)
    first_idx = findfirst(in_gate, times)
    isnothing(first_idx) && return Float64[]

    last_idx = findlast(in_gate, times)
    return collect(signal[first_idx:last_idx])
end

function sample_windowed_spectrum(
    signal::AbstractVector{<:Real},
    times::AbstractVector{<:Real},
    dt::Real,
    frequencies::AbstractVector{<:Real};
    window::AbstractString,
    gate_start::Union{Nothing,Real} = nothing,
    gate_end::Union{Nothing,Real} = nothing,
)
    segment = gated_signal_segment(signal, times; gate_start = gate_start, gate_end = gate_end)
    isempty(segment) && return zeros(ComplexF64, length(frequencies))
    return sample_spectrum(apply_fft_window(segment, window), dt, frequencies)
end

function compute_frequency_scattering_diagnostics(
    result::SimulationResult;
    frequencies::AbstractVector{<:Real},
    incident_monitor::AbstractString = "incident",
    reflected_monitor::AbstractString = "incident",
    transmitted_monitor::AbstractString = "transmitted",
    window::AbstractString = "hann",
    gate_start::Union{Nothing,Real} = nothing,
    gate_end::Union{Nothing,Real} = nothing,
)
    config = result.config

    haskey(result.monitor_traces, reflected_monitor) ||
        throw(ArgumentError("Monitor '$reflected_monitor' not found in simulation result."))
    haskey(result.monitor_traces, transmitted_monitor) ||
        throw(ArgumentError("Monitor '$transmitted_monitor' not found in simulation result."))

    dt_sample = length(result.times) > 1 ? result.times[2] - result.times[1] : config.save_every * config.courant_factor * config.dx / c0
    incident_normalization = incident_field_trace(config, result.times, incident_monitor)
    incident_reflected = incident_field_trace(config, result.times, reflected_monitor)
    incident_transmitted = incident_field_trace(config, result.times, transmitted_monitor)
    reflected_sig = result.monitor_traces[reflected_monitor] .- incident_reflected
    transmitted_sig = result.monitor_traces[transmitted_monitor]

    incident_spectrum = sample_windowed_spectrum(
        incident_normalization,
        result.times,
        dt_sample,
        frequencies;
        window = window,
        gate_start = gate_start,
        gate_end = gate_end,
    )
    reflected_spectrum = sample_windowed_spectrum(
        reflected_sig,
        result.times,
        dt_sample,
        frequencies;
        window = window,
        gate_start = gate_start,
        gate_end = gate_end,
    )
    transmitted_spectrum = sample_windowed_spectrum(
        transmitted_sig,
        result.times,
        dt_sample,
        frequencies;
        window = window,
        gate_start = gate_start,
        gate_end = gate_end,
    )
    transmitted_incident_spectrum = sample_windowed_spectrum(
        incident_transmitted,
        result.times,
        dt_sample,
        frequencies;
        window = window,
        gate_start = gate_start,
        gate_end = gate_end,
    )

    reflection = similar(collect(frequencies), Float64)
    transmission = similar(collect(frequencies), Float64)

    for i in eachindex(frequencies)
        reflection[i] = abs2(incident_spectrum[i]) > 0 ? abs2(reflected_spectrum[i] / incident_spectrum[i]) : 0.0
        transmission[i] = abs2(transmitted_incident_spectrum[i]) > 0 ?
            abs2(transmitted_spectrum[i] / transmitted_incident_spectrum[i]) : 0.0
    end

    return (
        frequencies = collect(frequencies),
        reflection = reflection,
        transmission = transmission,
        incident_spectrum = incident_spectrum,
        reflected_spectrum = reflected_spectrum,
        transmitted_spectrum = transmitted_spectrum,
        transmitted_reference_spectrum = transmitted_incident_spectrum,
        transmitted_incident_spectrum = transmitted_incident_spectrum,
        window = window,
        gate_start = gate_start,
        gate_end = gate_end,
    )
end

function safe_complex_ratio(numerator::Complex, denominator::Complex)
    return abs2(denominator) > 0 ? numerator / denominator : 0.0 + 0.0im
end

function port_shift_factor(frequencies, distance::Real, phase_velocity::Real)
    return exp.(1im .* (2π .* collect(frequencies) ./ phase_velocity) .* distance)
end

function deembed_sparameters(
    sparams;
    left_shift::Real = 0.0,
    right_shift::Real = 0.0,
    phase_velocity::Real = c0,
    metadata = NamedTuple(),
)
    left_reflection_factor = port_shift_factor(sparams.frequencies, 2 * left_shift, phase_velocity)
    right_reflection_factor = port_shift_factor(sparams.frequencies, 2 * right_shift, phase_velocity)
    transmission_factor = port_shift_factor(sparams.frequencies, left_shift + right_shift, phase_velocity)

    return (
        frequencies = sparams.frequencies,
        s11 = sparams.s11 .* left_reflection_factor,
        s21 = sparams.s21 .* transmission_factor,
        s12 = sparams.s12 .* transmission_factor,
        s22 = sparams.s22 .* right_reflection_factor,
        forward = sparams.forward,
        reverse = sparams.reverse,
        mirrored_result = sparams.mirrored_result,
        window = sparams.window,
        gate_start = sparams.gate_start,
        gate_end = sparams.gate_end,
        left_shift = left_shift,
        right_shift = right_shift,
        phase_velocity = phase_velocity,
        metadata = metadata,
    )
end

function compute_sparameters(
    result::SimulationResult;
    frequencies::AbstractVector{<:Real},
    incident_monitor::AbstractString = "incident",
    reflected_monitor::AbstractString = "incident",
    transmitted_monitor::AbstractString = "transmitted",
    port1_monitor::Union{Nothing,AbstractString} = nothing,
    port1_reflected_monitor::Union{Nothing,AbstractString} = nothing,
    port2_monitor::Union{Nothing,AbstractString} = nothing,
    window::AbstractString = "hann",
    gate_start::Union{Nothing,Real} = nothing,
    gate_end::Union{Nothing,Real} = nothing,
    left_reference_plane::Union{Nothing,Real} = nothing,
    right_reference_plane::Union{Nothing,Real} = nothing,
    port1_reference_plane::Union{Nothing,Real} = nothing,
    port2_reference_plane::Union{Nothing,Real} = nothing,
)
    port1_monitor_name = String(isnothing(port1_monitor) ? incident_monitor : port1_monitor)
    port1_reflected_monitor_name = String(isnothing(port1_reflected_monitor) ? reflected_monitor : port1_reflected_monitor)
    port2_monitor_name = String(isnothing(port2_monitor) ? transmitted_monitor : port2_monitor)

    forward = compute_frequency_scattering_diagnostics(
        result;
        frequencies = frequencies,
        incident_monitor = port1_monitor_name,
        reflected_monitor = port1_reflected_monitor_name,
        transmitted_monitor = port2_monitor_name,
        window = window,
        gate_start = gate_start,
        gate_end = gate_end,
    )

    mirrored_result = run_fdtd(mirror_config(result.config))
    reverse = compute_frequency_scattering_diagnostics(
        mirrored_result;
        frequencies = frequencies,
        incident_monitor = port1_monitor_name,
        reflected_monitor = port1_reflected_monitor_name,
        transmitted_monitor = port2_monitor_name,
        window = window,
        gate_start = gate_start,
        gate_end = gate_end,
    )

    s11 = [safe_complex_ratio(forward.reflected_spectrum[i], forward.incident_spectrum[i]) for i in eachindex(frequencies)]
    s21 = [safe_complex_ratio(forward.transmitted_spectrum[i], forward.transmitted_reference_spectrum[i]) for i in eachindex(frequencies)]
    s22 = [safe_complex_ratio(reverse.reflected_spectrum[i], reverse.incident_spectrum[i]) for i in eachindex(frequencies)]
    s12 = [safe_complex_ratio(reverse.transmitted_spectrum[i], reverse.transmitted_reference_spectrum[i]) for i in eachindex(frequencies)]

    raw = (
        frequencies = collect(frequencies),
        s11 = s11,
        s21 = s21,
        s12 = s12,
        s22 = s22,
        forward = forward,
        reverse = reverse,
        mirrored_result = mirrored_result,
        window = window,
        gate_start = gate_start,
        gate_end = gate_end,
    )

    port1_reference_value = if !isnothing(port1_reference_plane)
        Float64(port1_reference_plane)
    elseif !isnothing(left_reference_plane)
        Float64(left_reference_plane)
    else
        monitor_position(result.config, port1_monitor_name)
    end
    port2_reference_value = if !isnothing(port2_reference_plane)
        Float64(port2_reference_plane)
    elseif !isnothing(right_reference_plane)
        Float64(right_reference_plane)
    else
        monitor_position(result.config, port2_monitor_name)
    end

    port1_monitor_position = monitor_position(result.config, port1_monitor_name)
    port1_reflected_monitor_position = monitor_position(result.config, port1_reflected_monitor_name)
    port2_monitor_position = monitor_position(result.config, port2_monitor_name)
    left_shift = port1_reference_value - port1_monitor_position
    right_shift = port2_monitor_position - port2_reference_value
    phase_velocity = result.config.excitation isa TFSFExcitation ?
        incident_phase_velocity(result.config.excitation) : c0
    metadata = (
        port1_monitor = port1_monitor_name,
        port1_reflected_monitor = port1_reflected_monitor_name,
        port2_monitor = port2_monitor_name,
        port1_monitor_position = port1_monitor_position,
        port1_reflected_monitor_position = port1_reflected_monitor_position,
        port2_monitor_position = port2_monitor_position,
        port1_reference_plane = port1_reference_value,
        port2_reference_plane = port2_reference_value,
        port1_shift = left_shift,
        port2_shift = right_shift,
        phase_velocity = phase_velocity,
        window = window,
        gate_start = gate_start,
        gate_end = gate_end,
        nx = result.config.nx,
        dx = result.config.dx,
        nsteps = result.config.nsteps,
        save_every = result.config.save_every,
    )

    return deembed_sparameters(
        raw;
        left_shift = left_shift,
        right_shift = right_shift,
        phase_velocity = phase_velocity,
        metadata = metadata,
    )
end

function db20(values)
    return [20 * log10(max(abs(value), eps(Float64))) for value in values]
end

function phase_deg(values)
    return rad2deg.(angle.(values))
end

function unwrap_phase(phases::AbstractVector{<:Real})
    unwrapped = collect(Float64, phases)
    for i in 2:length(unwrapped)
        delta = unwrapped[i] - unwrapped[i - 1]
        if delta > π
            unwrapped[i:end] .-= 2π
        elseif delta < -π
            unwrapped[i:end] .+= 2π
        end
    end
    return unwrapped
end

function compute_sparameter_delays(sparams; parameter::Symbol = :s21)
    hasproperty(sparams, parameter) ||
        throw(ArgumentError("S-parameter set does not contain '$parameter'."))
    values = getproperty(sparams, parameter)
    frequencies = collect(Float64, sparams.frequencies)
    ω = 2π .* frequencies
    phase = unwrap_phase(angle.(values))
    phase_delay = similar(frequencies)
    group_delay = similar(frequencies)

    for i in eachindex(frequencies)
        phase_delay[i] = iszero(ω[i]) ? 0.0 : -phase[i] / ω[i]
    end

    if length(frequencies) == 1
        group_delay[1] = 0.0
    else
        group_delay[1] = -(phase[2] - phase[1]) / (ω[2] - ω[1])
        for i in 2:(length(frequencies) - 1)
            group_delay[i] = -(phase[i + 1] - phase[i - 1]) / (ω[i + 1] - ω[i - 1])
        end
        group_delay[end] = -(phase[end] - phase[end - 1]) / (ω[end] - ω[end - 1])
    end

    return (
        frequencies = frequencies,
        parameter = parameter,
        phase_rad = phase,
        phase_deg = rad2deg.(phase),
        phase_delay = phase_delay,
        group_delay = group_delay,
    )
end

function plot_sparameters(
    sparams;
    output::AbstractString = "sparameters.png",
    title_prefix::AbstractString = "S-parameters",
)
    mkpath(dirname(abspath(output)))
    frequencies_ghz = sparams.frequencies ./ 1e9

    p_mag = plot(
        frequencies_ghz,
        db20(sparams.s11);
        label = "S11",
        xlabel = "Frequency (GHz)",
        ylabel = "Magnitude (dB)",
        linewidth = 2,
        title = "$(title_prefix) magnitude",
    )
    plot!(p_mag, frequencies_ghz, db20(sparams.s21); label = "S21", linewidth = 2)
    plot!(p_mag, frequencies_ghz, db20(sparams.s12); label = "S12", linewidth = 2, linestyle = :dash)
    plot!(p_mag, frequencies_ghz, db20(sparams.s22); label = "S22", linewidth = 2, linestyle = :dashdot)

    p_phase = plot(
        frequencies_ghz,
        phase_deg(sparams.s11);
        label = "S11",
        xlabel = "Frequency (GHz)",
        ylabel = "Phase (deg)",
        linewidth = 2,
        title = "$(title_prefix) phase",
    )
    plot!(p_phase, frequencies_ghz, phase_deg(sparams.s21); label = "S21", linewidth = 2)
    plot!(p_phase, frequencies_ghz, phase_deg(sparams.s12); label = "S12", linewidth = 2, linestyle = :dash)
    plot!(p_phase, frequencies_ghz, phase_deg(sparams.s22); label = "S22", linewidth = 2, linestyle = :dashdot)

    combined = plot(p_mag, p_phase; layout = (2, 1), size = (950, 850))
    savefig(combined, output)
    return combined
end

function plot_sparameter_delays(
    delays;
    output::AbstractString = "sparameter_delays.png",
    title_prefix::AbstractString = "S-parameter delays",
)
    mkpath(dirname(abspath(output)))
    frequencies_ghz = delays.frequencies ./ 1e9
    phase_delay_ns = delays.phase_delay .* 1e9
    group_delay_ns = delays.group_delay .* 1e9
    parameter_label = uppercase(String(delays.parameter))

    p1 = plot(
        frequencies_ghz,
        phase_delay_ns;
        label = "$(parameter_label) phase delay",
        xlabel = "Frequency (GHz)",
        ylabel = "Phase delay (ns)",
        linewidth = 2,
        title = "$(title_prefix)",
    )

    p2 = plot(
        frequencies_ghz,
        group_delay_ns;
        label = "$(parameter_label) group delay",
        xlabel = "Frequency (GHz)",
        ylabel = "Group delay (ns)",
        linewidth = 2,
        title = "$(title_prefix)",
    )

    combined = plot(p1, p2; layout = (2, 1), size = (950, 850))
    savefig(combined, output)
    return combined
end

function touchstone_pair(value, format::AbstractString)
    normalized = lowercase(strip(format))
    if normalized == "ri"
        return real(value), imag(value)
    elseif normalized == "ma"
        return abs(value), rad2deg(angle(value))
    elseif normalized == "db"
        return 20 * log10(max(abs(value), eps(Float64))), rad2deg(angle(value))
    end
    throw(ArgumentError("Unknown Touchstone format '$format'. Supported values: RI, MA, DB."))
end

function touchstone_header(format::AbstractString, reference_impedance::Real)
    normalized = uppercase(strip(format))
    normalized in ("RI", "MA", "DB") ||
        throw(ArgumentError("Unknown Touchstone format '$format'. Supported values: RI, MA, DB."))
    return "# Hz S $(normalized) R $(Float64(reference_impedance))"
end

function write_touchstone_metadata(io, sparams)
    hasproperty(sparams, :metadata) || return nothing
    metadata = sparams.metadata
    isempty(keys(metadata)) && return nothing

    println(io, "! Port/reference metadata")
    for key in keys(metadata)
        value = getproperty(metadata, key)
        if !isnothing(value)
            println(io, "! $(key): $(value)")
        end
    end
    return nothing
end

function write_touchstone_s1p(
    sparams,
    path::AbstractString;
    parameter::Symbol = :s11,
    format::AbstractString = "RI",
    reference_impedance::Real = 50.0,
)
    hasproperty(sparams, parameter) ||
        throw(ArgumentError("S-parameter set does not contain '$parameter'."))
    values = getproperty(sparams, parameter)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "! Generated by FDTD1D.jl")
        write_touchstone_metadata(io, sparams)
        println(io, touchstone_header(format, reference_impedance))
        for i in eachindex(sparams.frequencies)
            v1, v2 = touchstone_pair(values[i], format)
            println(io, join((sparams.frequencies[i], v1, v2), " "))
        end
    end
    return path
end

function write_touchstone_s2p(
    sparams,
    path::AbstractString;
    format::AbstractString = "RI",
    reference_impedance::Real = 50.0,
)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "! Generated by FDTD1D.jl")
        write_touchstone_metadata(io, sparams)
        println(io, touchstone_header(format, reference_impedance))
        for i in eachindex(sparams.frequencies)
            s11_1, s11_2 = touchstone_pair(sparams.s11[i], format)
            s21_1, s21_2 = touchstone_pair(sparams.s21[i], format)
            s12_1, s12_2 = touchstone_pair(sparams.s12[i], format)
            s22_1, s22_2 = touchstone_pair(sparams.s22[i], format)
            println(
                io,
                join((
                    sparams.frequencies[i],
                    s11_1, s11_2,
                    s21_1, s21_2,
                    s12_1, s12_2,
                    s22_1, s22_2,
                ), " "),
            )
        end
    end
    return path
end

function analytical_slab_rt(
    frequencies::AbstractVector{<:Real};
    slab_eps_r::Real,
    slab_mu_r::Real = 1.0,
    thickness::Real,
    incident_eps_r::Real = 1.0,
    incident_mu_r::Real = 1.0,
    exit_eps_r::Real = incident_eps_r,
    exit_mu_r::Real = incident_mu_r,
)
    η1 = sqrt(mu0 * incident_mu_r / (eps0 * incident_eps_r))
    η2 = sqrt(mu0 * slab_mu_r / (eps0 * slab_eps_r))
    η3 = sqrt(mu0 * exit_mu_r / (eps0 * exit_eps_r))
    n2 = sqrt(slab_eps_r * slab_mu_r)

    r12 = (η2 - η1) / (η2 + η1)
    r23 = (η3 - η2) / (η3 + η2)
    #t12 = 2η2 / (η2 + η1)
    t12 = 1 - r12
    #t23 = 2η3 / (η3 + η2)
    t23 = 1 - r23

    reflection = zeros(length(frequencies))
    transmission = zeros(length(frequencies))
    reflection_amplitude = zeros(ComplexF64, length(frequencies))
    transmission_amplitude = zeros(ComplexF64, length(frequencies))

    for i in eachindex(frequencies)
        βd = 2π * frequencies[i] * n2 * thickness / c0
        phase = exp(-2im * βd)
        denominator = 1 + r12 * r23 * phase
        r = (r12 + r23 * phase) / denominator
        t = (t12 * t23 * exp(-1im * βd)) / denominator

        reflection_amplitude[i] = r
        transmission_amplitude[i] = t
        reflection[i] = abs2(r)
        transmission[i] = real(η1 / η3) * abs2(t)
    end

    return (
        frequencies = collect(frequencies),
        reflection = reflection,
        transmission = transmission,
        reflection_amplitude = reflection_amplitude,
        transmission_amplitude = transmission_amplitude,
    )
end

function analytical_dispersive_slab_rt(
    frequencies::AbstractVector{<:Real};
    model::Symbol,
    thickness::Real,
    slab_mu_r::Real = 1.0,
    incident_eps_r::Real = 1.0,
    incident_mu_r::Real = 1.0,
    exit_eps_r::Real = incident_eps_r,
    exit_mu_r::Real = incident_mu_r,
    eps_inf::Real = 1.0,
    eps_static::Union{Nothing,Real} = nothing,
    tau::Union{Nothing,Real} = nothing,
    omega_p::Union{Nothing,Real} = nothing,
    gamma::Union{Nothing,Real} = nothing,
    delta_eps::Union{Nothing,Real} = nothing,
    omega_0::Union{Nothing,Real} = nothing,
)
    η1 = sqrt(mu0 * incident_mu_r / (eps0 * incident_eps_r))
    η3 = sqrt(mu0 * exit_mu_r / (eps0 * exit_eps_r))
    reflection = zeros(length(frequencies))
    transmission = zeros(length(frequencies))
    reflection_amplitude = zeros(ComplexF64, length(frequencies))
    transmission_amplitude = zeros(ComplexF64, length(frequencies))
    eps_complex = zeros(ComplexF64, length(frequencies))

    for i in eachindex(frequencies)
        ω = 2π * frequencies[i]
        εr = if model === :debye
            isnothing(eps_static) && throw(ArgumentError("Debye analytical model requires eps_static."))
            isnothing(tau) && throw(ArgumentError("Debye analytical model requires tau."))
            eps_inf + (eps_static - eps_inf) / (1 + im * ω * tau)
        elseif model === :drude
            isnothing(omega_p) && throw(ArgumentError("Drude analytical model requires omega_p."))
            isnothing(gamma) && throw(ArgumentError("Drude analytical model requires gamma."))
            if iszero(ω)
                eps_inf - omega_p^2 / complex(0.0, -(gamma == 0 ? eps() : gamma * eps()))
            else
                eps_inf - omega_p^2 / (ω^2 - im * gamma * ω)
            end
        elseif model === :lorentz
            isnothing(delta_eps) && throw(ArgumentError("Lorentz analytical model requires delta_eps."))
            isnothing(omega_0) && throw(ArgumentError("Lorentz analytical model requires omega_0."))
            isnothing(gamma) && throw(ArgumentError("Lorentz analytical model requires gamma."))
            eps_inf + delta_eps * omega_0^2 / (omega_0^2 - ω^2 + im * gamma * ω)
        else
            throw(ArgumentError("Unsupported dispersive analytical model '$model'."))
        end

        eps_complex[i] = εr
        η2 = sqrt(mu0 * slab_mu_r / (eps0 * εr))
        n2 = sqrt(εr * slab_mu_r)
        r12 = (η2 - η1) / (η2 + η1)
        r23 = (η3 - η2) / (η3 + η2)
        #t12 = 2η2 / (η2 + η1)
        t12 = 1 - r12
        #t23 = 2η3 / (η3 + η2)
        t23 = 1 - r23
        βd = 2π * frequencies[i] * n2 * thickness / c0
        phase = exp(-2im * βd)
        denominator = 1 + r12 * r23 * phase
        r = (r12 + r23 * phase) / denominator
        t = (t12 * t23 * exp(-1im * βd)) / denominator

        reflection_amplitude[i] = r
        transmission_amplitude[i] = t
        reflection[i] = abs2(r)
        transmission[i] = real(η1 / η3) * abs2(t)
    end

    return (
        frequencies = collect(frequencies),
        reflection = reflection,
        transmission = transmission,
        reflection_amplitude = reflection_amplitude,
        transmission_amplitude = transmission_amplitude,
        eps_complex = eps_complex,
        model = model,
        thickness = thickness,
    )
end

function apply_excitation_h!(::SoftCurrentExcitation, H, config, time, dt, dx, μ)
    return H
end

function apply_excitation_h!(excitation::TFSFExcitation, H, config, time, dt, dx, μ)
    incident_dx = dx / incident_phase_velocity(excitation)
    e_inc_left = current_value(config.source, time)
    e_inc_right = current_value(config.source, time - (excitation.end_index - excitation.start_index + 1) * incident_dx)
    left_h = excitation.start_index - 1
    right_h = excitation.end_index

    H[left_h] += (dt / (μ[left_h] * dx)) * e_inc_left
    H[right_h] -= (dt / (μ[right_h] * dx)) * e_inc_right
    return H
end

function apply_excitation_e!(::SoftCurrentExcitation, E, config, time, dt, dx, ϵ)
    src_idx = resolve_source_position(config)
    inject_current_source!(E, config.source, time, dt, ϵ[src_idx], src_idx)
    return E
end

function apply_excitation_e!(excitation::TFSFExcitation, E, config, time, dt, dx, ϵ)
    incident_dx = dx / incident_phase_velocity(excitation)
    η_inc = incident_impedance(excitation)
    h_inc_left = current_value(config.source, time + 0.5 * dt + 0.5 * incident_dx) / η_inc
    h_inc_right = current_value(
        config.source,
        time + 0.5 * dt - (excitation.end_index - excitation.start_index + 0.5) * incident_dx,
    ) / η_inc
    left_e = excitation.start_index
    right_e = excitation.end_index + 1

    E[left_e] += (dt / (ϵ[left_e] * dx)) * h_inc_left
    E[right_e] -= (dt / (ϵ[right_e] * dx)) * h_inc_right
    return E
end

function run_fdtd(config::SimulationConfig)
    validate(config)

    dt = config.courant_factor * config.dx / c0
    ϵ = electric_permittivity_profile(config.material, config.nx)
    μ = magnetic_permeability_profile(config.material, config.nx)
    σe_material = electric_conductivity_profile(config.material, config.nx)
    σe_pml, σm_pml = pml_profiles(config, ϵ, μ)
    σe = σe_material .+ σe_pml
    material_state = initialize_material_state(config.material, config.nx, dt)

    x = collect(0:config.dx:(config.nx - 1) * config.dx)
    E = zeros(config.nx)
    H = zeros(config.nx - 1)
    boundary_cache = BoundaryCache()
    update_boundary_cache!(boundary_cache, E)
    cfl_profile = (
        left = dt / (config.dx * sqrt(μ[1] * ϵ[2])),
        right = dt / (config.dx * sqrt(μ[end] * ϵ[end - 1])),
    )

    nsaved = cld(config.nsteps, config.save_every)
    e_history = zeros(config.nx, nsaved)
    h_history = zeros(config.nx - 1, nsaved)
    times = zeros(nsaved)
    monitor_traces = initialize_monitor_traces(config, nsaved)

    save_idx = 1
    for n in 1:config.nsteps
        update_h!(H, E, dt, config.dx, μ, σm_pml)

        t = (n - 1) * dt
        apply_excitation_h!(config.excitation, H, config, t, dt, config.dx, μ)
        if isnothing(material_state)
            update_e!(E, H, dt, config.dx, ϵ, σe)
        else
            update_e!(E, H, dt, config.dx, ϵ, σe, material_state)
        end
        apply_excitation_e!(config.excitation, E, config, t, dt, config.dx, ϵ)
        apply_boundary!(config.left_boundary, E, :left, boundary_cache, cfl_profile)
        apply_boundary!(config.right_boundary, E, :right, boundary_cache, cfl_profile)
        update_boundary_cache!(boundary_cache, E)

        if n % config.save_every == 0 || n == config.nsteps
            e_history[:, save_idx] .= E
            h_history[:, save_idx] .= H
            times[save_idx] = n * dt
            save_monitor_values!(monitor_traces, config.monitors, E, save_idx)
            save_idx += 1
        end
    end

    last_idx = save_idx - 1
    return SimulationResult(
        x,
        times[1:last_idx],
        e_history[:, 1:last_idx],
        h_history[:, 1:last_idx],
        Dict(name => trace[1:last_idx] for (name, trace) in monitor_traces),
        config,
    )
end

function animate_field(
    result::SimulationResult;
    field::Symbol = :E,
    output::AbstractString = "fdtd_1d.gif",
    fps::Int = 20,
    title_prefix::AbstractString = "1D FDTD",
    show_material::Bool = true,
    label_materials::Bool = true,
    label_boundaries::Bool = true,
)
    mkpath(dirname(abspath(output)))
    history = field === :E ? result.e_history :
              field === :H ? result.h_history :
              throw(ArgumentError("field must be :E or :H"))

    x = field === :E ? result.x : @view(result.x[1:end-1]) .+ result.config.dx / 2
    ylabel = field === :E ? "Electric field (V/m)" : "Magnetic field (A/m)"
    ylimit = maximum(abs, history)
    ylimit = iszero(ylimit) ? 1.0 : 1.1 * ylimit
    source_markers = excitation_markers(result.config, result.x)
    regions = show_material ? material_regions(result.config.material, result.x, result.config.dx) : NamedTuple[]
    eps_values = isempty(regions) ? [1.0] : [region.eps_r for region in regions]
    eps_min = minimum(eps_values)
    eps_max = maximum(eps_values)
    left_boundary_text = boundary_label(result.config.left_boundary)
    right_boundary_text = boundary_label(result.config.right_boundary)
    excitation_text = excitation_label(result.config.excitation)

    anim = Plots.Animation()
    for frame in axes(history, 2)
        Plots.plot(
            x,
            history[:, frame];
            xlabel = "x (m)",
            ylabel = ylabel,
            title = string(title_prefix, " - ", field, " at t = ", round(result.times[frame] * 1e9; digits = 3), " ns"),
            linewidth = 2,
            legend = show_material && !isempty(regions),
            ylim = (-ylimit, ylimit),
        )
        for region in regions
            color = region_color(region, eps_min, eps_max)
            alpha = region.sigma_e > 0 ? 0.18 : 0.12
            Plots.vspan!([region.x0, region.x1]; color = color, alpha = alpha, label = region_label(region))
            if label_materials
                Plots.annotate!(
                    region.x_center,
                    0.84 * ylimit,
                    Plots.text(region_label(region), 8, color, :center),
                )
            end
        end
        if !isempty(source_markers)
            Plots.vline!(source_markers; linestyle = :dash, linewidth = 1.5, color = :black)
        end
        if label_boundaries
            Plots.annotate!(
                result.x[1] + 0.03 * (result.x[end] - result.x[1]),
                0.94 * ylimit,
                Plots.text("Left: $left_boundary_text", 9, :black, :left),
            )
            Plots.annotate!(
                result.x[end] - 0.03 * (result.x[end] - result.x[1]),
                0.94 * ylimit,
                Plots.text("Right: $right_boundary_text", 9, :black, :right),
            )
            Plots.annotate!(
                result.x[1] + 0.5 * (result.x[end] - result.x[1]),
                0.94 * ylimit,
                Plots.text(excitation_text, 9, :black, :center),
            )
        end
        Plots.frame(anim)
    end

    Plots.gif(anim, output; fps = fps)
    return output
end

end
