# FDTD1D Documentation

FDTD1D is a Julia implementation of a one-dimensional finite-difference time-domain solver for Maxwell's equations. It is intended for education, rapid experimentation, and extension toward more advanced electromagnetic modeling.

The code supports:

- soft electric-current excitation and total-field/scattered-field excitation
- PEC, first-order Mur ABC, and conductivity-graded PML boundaries
- vacuum, nondispersive layered media, lossy media, and Debye/Drude/Lorentz dispersive media
- field monitors, animations, spectra, S-parameters, and Touchstone export
- TOML-based simulation setup with distances specified by index or meters
- managed output directories containing configs, plots, CSV exports, and summaries

## Documentation Map

- [Theory](theory.md): Maxwell equations, Yee grid, update equations, boundary conditions, TF/SF, materials, spectra, and S-parameters.
- [User Guide](user_guide.md): installation, running examples, TOML configuration, output files, and recommended workflows.
- [Examples](examples.md): overview of the provided configurations and scripts.
- [Developer Guide](developer_guide.md): module organization, extension points, validation tests, and coding conventions.
- [API Reference](api_reference.md): important public types and functions.

## Repository Layout

```text
.
├── Project.toml
├── Manifest.toml
├── README.md
├── src/FDTD1D.jl
├── config/
├── examples/
├── test/
├── docs/
└── outputs/              # generated, ignored by git
```

## Quick Start

Instantiate dependencies:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run a configuration:

```bash
julia --project=. examples/run_from_config.jl config/default.toml
```

Run all tests:

```bash
julia --project=. test/runtests.jl
```

The standard runner writes results into a timestamped directory under `outputs/`, for example:

```text
outputs/default_YYYYMMDD_HHMMSS/
├── config.toml
├── fdtd_from_config.gif
├── material_profile.csv
├── monitor_traces.csv
└── summary.toml
```

