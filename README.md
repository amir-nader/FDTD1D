# FDTD1D: 1D Electromagnetic FDTD in Julia

FDTD1D is a modular one-dimensional finite-difference time-domain solver for Maxwell's equations. It is designed for educational electromagnetic simulations, analytical validation studies, and further development.

## Features

- 1D Yee-grid FDTD solver for transverse electromagnetic plane-wave propagation
- soft electric-current excitation and total-field/scattered-field excitation
- PEC, first-order Mur ABC, and conductivity-graded PML boundaries
- vacuum, layered dielectric/lossy media, and Debye/Drude/Lorentz dispersive media
- TOML configuration files with positions specified by grid index or meters
- field animations, monitor traces, material-profile export, spectra, S-parameters, and Touchstone export
- run-output manager that stores configs, CSV data, plots, animations, and summaries together
- analytical validation examples for slabs, half-spaces, dispersive materials, and absorbing boundaries

## Documentation

Start here:

- `docs/index.md`: documentation overview
- `docs/theory.md`: equations, FDTD updates, materials, boundaries, TF/SF, spectra, and S-parameters
- `docs/user_guide.md`: installation, running simulations, TOML configuration, and outputs
- `docs/examples.md`: example scripts and configuration files
- `docs/developer_guide.md`: code organization, extension points, and testing strategy
- `docs/api_reference.md`: public types and functions

## Repository Layout

```text
.
├── Project.toml
├── Manifest.toml
├── src/FDTD1D.jl
├── config/
├── examples/
├── test/
├── docs/
└── outputs/              # generated results, ignored by git
```

## Quick Start

Install/instantiate dependencies:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run a simple configuration:

```bash
julia --project=. examples/run_from_config.jl config/default.toml
```

Run all validation tests:

```bash
julia --project=. test/runtests.jl
```

## Basic Julia Usage

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
animate_field(result; field = :E, output = "outputs/manual/e_field.gif")
```

## TOML Workflow

The recommended workflow is to edit a TOML file in `config/` and run:

```bash
julia --project=. examples/run_from_config.jl config/tfsf_slab_sparameters.toml
```

The runner creates a timestamped output directory such as:

```text
outputs/tfsf_slab_sparameters_YYYYMMDD_HHMMSS/
├── config.toml
├── monitor_traces.csv
├── material_profile.csv
├── summary.toml
├── spectrum.csv
├── sparameters.csv
├── slab_sparameters.png
├── slab_sparameters.s2p
└── tfsf_slab_sparameters.gif
```

Set `use_run_directory = false` in `[output]` if exact file paths from the config should be used instead of a managed run directory.

## Common Runs

Boundary examples:

```bash
julia --project=. examples/run_from_config.jl config/abc_gaussian.toml
julia --project=. examples/run_from_config.jl config/pml_gaussian.toml
julia --project=. examples/compare_boundary_absorption.jl
```

Layered and dispersive media:

```bash
julia --project=. examples/run_from_config.jl config/layered_dielectric.toml
julia --project=. examples/run_from_config.jl config/debye_slab.toml
julia --project=. examples/run_from_config.jl config/drude_slab.toml
julia --project=. examples/run_from_config.jl config/lorentz_slab.toml
```

Analytical comparisons:

```bash
julia --project=. examples/compare_slab_analytical.jl config/tfsf_slab_analytical_compare.toml
julia --project=. examples/compare_halfspace_reflection.jl config/tfsf_halfspace_reflection.toml
julia --project=. examples/compare_dispersive_slab_analytical.jl config/tfsf_debye_slab_analytical_compare.toml
```

S-parameters:

```bash
julia --project=. examples/run_from_config.jl config/tfsf_slab_sparameters.toml
julia --project=. examples/run_from_config.jl config/tfsf_slab_sparameters_advanced.toml
julia --project=. examples/run_from_config.jl config/tfsf_slab_s11_only.toml
```

## Development Notes

- Keep generated files under `outputs/`.
- Add new validation tests in `test/runtests.jl` when changing numerical behavior.
- Prefer extending through Julia dispatch: new source, boundary, material, or diagnostic types should be implemented with dedicated methods.
- See `docs/developer_guide.md` before adding new physics.

