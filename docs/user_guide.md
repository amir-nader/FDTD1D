# User Guide

## Installation

From the project root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the regression tests:

```bash
julia --project=. test/runtests.jl
```

## Running a Simulation

The recommended workflow is TOML based:

```bash
julia --project=. examples/run_from_config.jl config/default.toml
```

By default, results are written to a managed run directory:

```text
outputs/<case_name>_<timestamp>/
├── config.toml
├── <animation>.gif
├── material_profile.csv
├── monitor_traces.csv
└── summary.toml
```

If frequency diagnostics are enabled, additional files such as `spectrum.csv`, S-parameter plots, delay plots, and Touchstone files are also written.

## Direct Julia Usage

```julia
using Pkg
Pkg.activate(".")
using FDTD1D

config = SimulationConfig(
    nx = 401,
    dx = 1e-3,
    courant_factor = 0.99,
    nsteps = 1200,
    source_position = 201,
    source = GaussianCurrent(1.0, 40e-12, 12e-12),
    left_boundary = PECBoundary(),
    right_boundary = PECBoundary(),
    save_every = 4,
)

result = run_fdtd(config)
animate_field(result; output = "outputs/manual_run/e_field.gif")
```

## TOML Configuration

Most examples are in `config/`. A configuration can contain:

- `[simulation]`: grid, time-step factor, run length, save rate, soft-source location
- `[source]`: source waveform and parameters
- `[excitation]`: `soft` or `tfsf`
- `[boundary]`: left/right boundary settings
- `[material]`: material model and layer definitions
- `[[monitors]]`: named field probes
- `[diagnostics]`: reflection/transmission, spectra, S-parameters
- `[output]`: animation and run-directory settings

## Spatial Coordinates

Most spatial settings support either an index or a physical coordinate in meters:

```toml
source_position = 201
```

or

```toml
source_x = 0.20
```

Use one form for a given item, not both.

## Sources

Supported source types:

```toml
[source]
type = "gaussian"
amplitude = 1.0
t0 = 4.0e-11
spread = 1.2e-11
```

Other types:

- `cosine`
- `gaussian_modulated_cosine`
- `modulated_cosine`
- `ricker`

## Excitation Types

### Soft Source

```toml
[excitation]
type = "soft"
```

The source is injected at `source_position` or `source_x`.

### TF/SF

```toml
[excitation]
type = "tfsf"
start_x = 0.08
end_x = 0.33
incident_eps_r = 1.0
incident_mu_r = 1.0
```

TF/SF is recommended for reflection, transmission, slab, half-space, and S-parameter studies.

## Boundaries

Supported boundary names:

- `pec`
- `abc`
- `pml`
- `none`

Example PML:

```toml
[boundary]
left = "pml"
right = "pml"
pml_thickness = 0.04
pml_order = 3
pml_target_reflection = 1.0e-6
```

Side-specific settings are also supported:

```toml
left_pml_thickness = 0.04
right_pml_cells = 30
```

## Materials

### Vacuum

```toml
[material]
type = "vacuum"
```

### Layered Nondispersive Material

```toml
[material]
type = "layered"

[[material.layers]]
start_x = 0.18
end_x = 0.24
eps_r = 4.0
mu_r = 1.0
sigma_e = 0.0
```

Layer geometry can use:

- `start_index` / `end_index`
- `start_x` / `end_x`
- `center_x` with `thickness`

### Debye

```toml
[material]
type = "debye"

[[material.layers]]
center_x = 0.25
thickness = 0.06
eps_inf = 2.0
eps_static = 4.5
tau = 6.0e-11
```

### Drude

```toml
[material]
type = "drude"

[[material.layers]]
center_x = 0.25
thickness = 0.06
eps_inf = 2.0
omega_p = 1.8849555922e10
gamma = 1.2566370614e9
```

### Lorentz

```toml
[material]
type = "lorentz"

[[material.layers]]
center_x = 0.25
thickness = 0.06
eps_inf = 1.5
delta_eps = 2.5
omega_0 = 1.1309733553e10
gamma = 7.5398223686e8
```

## Monitors

```toml
[[monitors]]
name = "incident"
x = 0.13

[[monitors]]
name = "reflected"
x = 0.06

[[monitors]]
name = "transmitted"
x = 0.29
```

For TF/SF reflection studies, place `reflected` in the left scattered-field region.

## Frequency Diagnostics

```toml
[diagnostics]
enabled = true
incident_monitor = "incident"
reflected_monitor = "reflected"
transmitted_monitor = "transmitted"
window = "none"
frequency_min = 2.0e8
frequency_max = 3.0e9
frequency_count = 121
```

For broadband analytical comparisons, start with `window = "none"` and no time gate. Use `gate_start` and `gate_end` only when the selected segment fully contains the event being transformed.

## S-Parameters

```toml
[diagnostics]
enabled = true
sparameters = true
port1_monitor = "incident"
port1_reflected_monitor = "reflected"
port2_monitor = "transmitted"
port1_reference_plane_x = 0.18
port2_reference_plane_x = 0.24
sparameter_plot = "slab_sparameters.png"
s2p_file = "slab_sparameters.s2p"
touchstone_format = "RI"
```

The runner writes:

- S-parameter plot
- `sparameters.csv`
- Touchstone `.s2p` or `.s1p`
- port metadata in the Touchstone comments

## Output Management

Useful `[output]` settings:

```toml
[output]
case_name = "my_case"
output_root = "outputs"
timestamped = true
use_run_directory = true
file = "field_animation.gif"
fps = 20
show_material = true
label_materials = true
label_boundaries = true
```

If `use_run_directory = false`, paths in the config are used directly.

