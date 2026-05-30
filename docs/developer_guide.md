# Developer Guide

This guide describes how the solver is organized and how to extend it safely.

## Module Structure

The implementation is contained in `src/FDTD1D.jl`. The main groups are:

- type definitions for materials, boundaries, sources, excitations, monitors, and results
- TOML parsing utilities
- output manager utilities
- validation and profile generation
- field update kernels
- boundary and PML profiles
- TF/SF source corrections
- diagnostics, FFTs, S-parameters, and Touchstone export
- analytical comparison helpers
- plotting and animation

## Core Types

`SimulationConfig` is the central input object. It contains:

- grid size and spacing
- Courant factor and number of time steps
- source waveform and excitation type
- boundary types
- material model
- monitor list
- save interval

`SimulationResult` contains:

- grid coordinates
- saved times
- electric and magnetic field histories
- monitor traces
- a copy of the config

## Time-Stepping Flow

`run_fdtd(config)` performs:

1. validate the configuration
2. compute `Δt`, material profiles, PML profiles, and dispersive state
3. allocate `E`, `H`, histories, and monitor traces
4. update `H`
5. apply TF/SF magnetic corrections
6. update `E`
7. apply TF/SF electric corrections or soft current source
8. apply boundary conditions
9. save fields and monitor values

When adding features, preserve this order unless there is a clear numerical reason to change it.

## Adding a Source

Add a subtype of `CurrentSource`:

```julia
struct MySource <: CurrentSource
    amplitude::Float64
end
```

Define:

```julia
current_value(src::MySource, t) = ...
```

Add TOML support in `source_from_dict`.

## Adding a Boundary Condition

Add a subtype of `BoundaryCondition`:

```julia
struct MyBoundary <: BoundaryCondition end
```

Define:

```julia
apply_boundary!(::MyBoundary, E, side::Symbol, boundary_cache, cfl_profile)
```

Add TOML support in `boundary_from_name` or `boundary_from_dict`.

For advanced absorbing boundaries, consider whether additional state must be initialized before the time loop.

## Adding a Material Model

Add a subtype of `AbstractMaterial`, then implement:

```julia
electric_permittivity_profile(material, nx)
magnetic_permeability_profile(material, nx)
electric_conductivity_profile(material, nx)
```

For dispersive media, also implement:

```julia
initialize_material_state(material, nx, dt)
update_e!(E, H, dt, dx, ϵ, σe, state)
```

Add TOML parsing and material-region plotting support if needed.

## Adding Diagnostics

Diagnostics should usually operate on `SimulationResult` without modifying the solver state.

Recommended pattern:

1. add a pure function that accepts `SimulationResult`
2. add optional TOML parsing in `examples/run_from_config.jl`
3. export CSV or plot data through the output manager
4. add a small regression test

## Output Manager

The output manager groups all artifacts from a run in one directory:

- copied config
- monitor traces
- material profile
- animation
- spectra
- S-parameters
- summary file

New export functions should follow the existing pattern:

```julia
write_new_artifact(path, data)
```

and should create parent directories with `mkpath(dirname(abspath(path)))`.

## Testing Strategy

The tests in `test/runtests.jl` include:

- nondispersive analytical validations
- boundary absorption checks
- S-parameter metadata checks
- output manager checks
- dispersive slab regressions

When adding features:

- add a focused test for the new behavior
- prefer analytical comparisons where possible
- use short simulations for unit-style tests
- keep generated files under `outputs/` or temporary directories

Run:

```bash
julia --project=. test/runtests.jl
```

## Numerical Cautions

- Use `courant_factor <= 1` unless you have verified stability.
- Broadband FFT comparisons are sensitive to time windows and source bandwidth.
- For analytical slab comparisons, prefer `window = "none"` unless the gate fully contains the relevant event.
- Keep PML regions outside the device under test.
- Put reflected monitors in the scattered-field region for TF/SF reflection studies.
- Dispersive media can require smaller time steps for stability and accuracy.

## Suggested Future Extensions

- CPML with auxiliary convolution variables
- higher-order Mur or other absorbing boundaries
- multi-pole Debye/Drude/Lorentz materials
- HDF5/JLD2 export for large field histories
- a command-line interface with subcommands
- documentation generated with Documenter.jl

