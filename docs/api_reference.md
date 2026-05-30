# API Reference

This is a high-level reference for the public API exported by `FDTD1D`.

## Configuration and Results

### `SimulationConfig`

Main simulation input object.

Important fields:

- `nx`: number of electric-field grid nodes
- `dx`: grid spacing in meters
- `courant_factor`: `Δt c0 / Δx`
- `nsteps`: number of time steps
- `source_position`: optional source index for soft excitation
- `source`: current-source waveform
- `excitation`: `SoftCurrentExcitation` or `TFSFExcitation`
- `left_boundary`, `right_boundary`: boundary condition objects
- `material`: material model
- `monitors`: vector of `FieldMonitor`
- `save_every`: field/history save interval

### `SimulationResult`

Returned by `run_fdtd`.

Fields:

- `x`
- `times`
- `e_history`
- `h_history`
- `monitor_traces`
- `config`

## Materials

- `Vacuum()`
- `GridMaterial(eps_r_e, mu_r_h, sigma_e)`
- `DebyeMaterial(eps_inf_e, eps_static_e, tau_e, mu_r_h, sigma_e)`
- `DrudeMaterial(eps_inf_e, omega_p_e, gamma_e, mu_r_h, sigma_e)`
- `LorentzMaterial(eps_inf_e, delta_eps_e, omega_0_e, gamma_e, mu_r_h, sigma_e)`

Profile functions:

- `electric_permittivity_profile(material, nx)`
- `magnetic_permeability_profile(material, nx)`
- `electric_conductivity_profile(material, nx)`

## Boundaries

- `PECBoundary()`
- `MurABCBoundary()`
- `PMLBoundary(ncells, grading_order, sigma_max, target_reflection)`
- `NoBoundaryCondition()`

## Sources

- `GaussianCurrent(amplitude, t0, spread)`
- `CosineCurrent(amplitude, frequency, phase)`
- `GaussianModulatedCosineCurrent(amplitude, frequency, t0, spread, phase)`
- `RickerCurrent(amplitude, frequency, t0)`

Use:

```julia
current_value(source, t)
```

to evaluate a source waveform.

## Excitations

- `SoftCurrentExcitation()`
- `TFSFExcitation(start_index, end_index, incident_eps_r, incident_mu_r)`

TF/SF helper functions:

- `incident_impedance(excitation)`
- `incident_phase_velocity(excitation)`

## Monitors

```julia
FieldMonitor(name, index)
```

Monitor traces are stored in:

```julia
result.monitor_traces[name]
```

## Running and Plotting

### `run_fdtd(config)`

Runs the solver and returns a `SimulationResult`.

### `animate_field(result; kwargs...)`

Creates a GIF animation of `E` or `H`.

Common keyword arguments:

- `field = :E`
- `output = "fdtd_1d.gif"`
- `fps = 20`
- `title_prefix = "1D FDTD"`
- `show_material = true`
- `label_materials = true`
- `label_boundaries = true`

## TOML Parsing

### `load_simulation_parameters(path)`

Returns a named tuple:

```julia
(
    config = SimulationConfig(...),
    output = Dict(...),
    diagnostics = Dict(...),
)
```

Other parsing helpers:

- `source_from_dict`
- `boundary_from_name`

## Diagnostics

### `compute_scattering_diagnostics(result; kwargs...)`

Computes time-domain energy ratios from monitor traces.

### `compute_frequency_scattering_diagnostics(result; frequencies, kwargs...)`

Computes frequency-domain reflection and transmission spectra.

Important keyword arguments:

- `incident_monitor`
- `reflected_monitor`
- `transmitted_monitor`
- `window`
- `gate_start`
- `gate_end`

### `compute_sparameters(result; frequencies, kwargs...)`

Computes complex `S11`, `S21`, `S12`, and `S22`.

Important keyword arguments:

- `port1_monitor`
- `port1_reflected_monitor`
- `port2_monitor`
- `port1_reference_plane`
- `port2_reference_plane`
- `window`
- `gate_start`
- `gate_end`

### `compute_sparameter_delays(sparams; parameter = :s21)`

Computes phase delay and group delay.

## Export

Output-manager helpers:

- `create_output_manager`
- `output_path`
- `copy_input_config!`
- `write_run_summary`
- `write_monitor_traces_csv`
- `write_material_profile_csv`
- `write_spectrum_csv`
- `write_sparameters_csv`

Touchstone helpers:

- `write_touchstone_s1p`
- `write_touchstone_s2p`

Plot helpers:

- `plot_sparameters`
- `plot_sparameter_delays`

## Analytical Helpers

### `analytical_slab_rt(frequencies; kwargs...)`

Computes normal-incidence reflection and transmission for a homogeneous slab.

### `analytical_dispersive_slab_rt(frequencies; kwargs...)`

Computes normal-incidence reflection and transmission for Debye, Drude, or Lorentz slabs.

