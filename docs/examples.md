# Examples

All examples assume the project root as the working directory.

## General Runner

Use the TOML runner for most simulations:

```bash
julia --project=. examples/run_from_config.jl config/default.toml
```

This creates a timestamped run directory under `outputs/`.

## Basic Sources

Run predefined source demonstrations:

```bash
julia --project=. examples/run_demo.jl
```

This runs Gaussian, cosine, Gaussian-modulated cosine, and Ricker sources.

## Boundary Examples

PEC default:

```bash
julia --project=. examples/run_from_config.jl config/default.toml
```

Mur ABC:

```bash
julia --project=. examples/run_from_config.jl config/abc_gaussian.toml
```

PML:

```bash
julia --project=. examples/run_from_config.jl config/pml_gaussian.toml
```

Boundary comparison:

```bash
julia --project=. examples/compare_boundary_absorption.jl
```

## Material Examples

Layered dielectric:

```bash
julia --project=. examples/run_from_config.jl config/layered_dielectric.toml
```

Lossy dielectric:

```bash
julia --project=. examples/run_from_config.jl config/lossy_layered_dielectric.toml
```

Position-based layer definition:

```bash
julia --project=. examples/run_from_config.jl config/center_thickness_slab.toml
```

## TF/SF Examples

Gaussian TF/SF injection:

```bash
julia --project=. examples/run_from_config.jl config/tfsf_gaussian.toml
```

Dielectric slab:

```bash
julia --project=. examples/run_from_config.jl config/tfsf_dielectric_slab.toml
```

Position-based slab:

```bash
julia --project=. examples/run_from_config.jl config/position_based_tfsf_slab.toml
```

## Analytical Comparisons

Lossless slab:

```bash
julia --project=. examples/compare_slab_analytical.jl config/tfsf_slab_analytical_compare.toml
```

Half-space reflection:

```bash
julia --project=. examples/compare_halfspace_reflection.jl config/tfsf_halfspace_reflection.toml
```

Dispersive slabs:

```bash
julia --project=. examples/compare_dispersive_slab_analytical.jl config/tfsf_debye_slab_analytical_compare.toml
julia --project=. examples/compare_dispersive_slab_analytical.jl config/tfsf_drude_slab_analytical_compare.toml
julia --project=. examples/compare_dispersive_slab_analytical.jl config/tfsf_lorentz_slab_analytical_compare.toml
```

## S-Parameters

Two-port S-parameters:

```bash
julia --project=. examples/run_from_config.jl config/tfsf_slab_sparameters.toml
```

Advanced S-parameters with reference-plane shifts and delay plots:

```bash
julia --project=. examples/run_from_config.jl config/tfsf_slab_sparameters_advanced.toml
```

One-port S11 export:

```bash
julia --project=. examples/run_from_config.jl config/tfsf_slab_s11_only.toml
```

## Tests

Run all validation tests:

```bash
julia --project=. test/runtests.jl
```

The tests cover:

- vacuum TF/SF propagation
- dielectric half-space reflection
- lossless dielectric slab spectra
- ABC/PML absorption
- S-parameter metadata and Touchstone export
- output manager artifacts
- Debye, Drude, and Lorentz slab comparisons

