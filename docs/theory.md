# Theory

This page summarizes the equations and numerical methods implemented in `src/FDTD1D.jl`.

## 1D Maxwell System

The solver models a normally incident plane wave propagating along `x`. The electric field has one transverse component and the magnetic field has the orthogonal transverse component:

```text
E = E_z(x,t) ẑ
H = H_y(x,t) ŷ
```

With the sign convention used by the implementation, the 1D curl equations are

```text
∂H/∂t = -(1/μ) ∂E/∂x
∂D/∂t = -∂H/∂x - J_s
```

For a nondispersive conducting medium,

```text
D = εE
J_loss = σE
```

so the electric-field equation is written in update form as

```text
ε ∂E/∂t + σE = -∂H/∂x - J_s .
```

The wave impedance and phase velocity in a homogeneous nondispersive medium are

```text
η = sqrt(μ/ε)
v = 1/sqrt(με)
```

where `ε = ε0 εr` and `μ = μ0 μr`.

## Yee Grid

FDTD1D uses the standard staggered Yee layout:

- electric fields `E[i]` live on integer grid points `x_i = (i-1)Δx`
- magnetic fields `H[i]` live halfway between electric nodes
- `E` and `H` are staggered in time by half a time step

The explicit time step is

```text
Δt = S Δx / c0
```

where `S` is `courant_factor`. For vacuum 1D propagation, stability requires approximately `S ≤ 1`. In inhomogeneous media, the user should choose a Courant factor that is safe for the fastest wave speed in the model.

## Nondispersive Update Equations

The magnetic-field update is

```text
H_i^{n+1/2} = d_a H_i^{n-1/2}
              - d_b (E_{i+1}^n - E_i^n)
```

with

```text
loss_m = σm_i Δt / (2 μ_i)
d_a = (1 - loss_m) / (1 + loss_m)
d_b = (Δt / (μ_i Δx)) / (1 + loss_m)
```

The electric-field update is

```text
E_i^{n+1} = c_a E_i^n
            - c_b (H_i^{n+1/2} - H_{i-1}^{n+1/2})
```

with

```text
loss_e = σe_i Δt / (2 ε_i)
c_a = (1 - loss_e) / (1 + loss_e)
c_b = (Δt / (ε_i Δx)) / (1 + loss_e)
```

This trapezoidal treatment of conductivity gives better behavior than a purely explicit loss term and is also used for the conductivity-graded PML profiles.

## Sources

### Soft Current Source

For `SoftCurrentExcitation`, the source is injected into the electric update at one node:

```text
E_i ← E_i - (Δt/ε_i) J_s(t)
```

Available source waveforms include Gaussian, cosine, Gaussian-modulated cosine, and Ricker pulses.

### TF/SF Plane-Wave Source

For `TFSFExcitation`, the incident plane wave is introduced at the total-field/scattered-field boundaries. The incident electric waveform is generated from the configured source, and the corresponding magnetic waveform is

```text
H_inc = E_inc / η_inc
```

where

```text
η_inc = sqrt(μ_inc/ε_inc)
```

The left and right TF/SF corrections are applied to the adjacent `H` and `E` nodes so the total-field region contains the incident plus scattered fields, while the surrounding scattered-field regions contain only scattered fields.

The diagnostics exploit this by deriving the incident signal analytically from the configured TF/SF waveform rather than from a separate reference simulation.

## Boundary Conditions

### PEC

Perfect electric conductor boundaries enforce

```text
E = 0
```

at the boundary node.

### First-Order Mur ABC

The Mur absorbing boundary approximates a one-way outgoing wave equation at the edge. In the implementation the update is

```text
E_boundary^{n+1} = E_inner^n + κ (E_inner^{n+1} - E_boundary^n)
```

where

```text
κ = (S_edge - 1)/(S_edge + 1)
S_edge = Δt / (Δx sqrt(με))
```

Mur ABC is simple and inexpensive, but it is less broadband and less robust than PML.

### Conductivity-Graded PML

The current PML is implemented as matched electric and magnetic conductivity profiles. The electric conductivity increases gradually into the boundary region:

```text
σ_e(x) = σ_max (depth)^m
```

where `m` is `pml_order`. The magnetic conductivity is scaled to preserve impedance matching:

```text
σ_m = σ_e μ/ε
```

If `pml_sigma_max` is not specified, the code estimates it from the target reflection:

```text
σ_max = - (m+1) ln(R_target) / (2 η d)
```

where `d` is the PML thickness and `η` is the edge impedance.

## Material Models

### Nondispersive Layered Media

The `GridMaterial` model stores per-cell profiles:

```text
ε_i = ε0 εr_i
μ_i = μ0 μr_i
σ_i = σe_i
```

Layer geometry can be specified with indices or meters.

### Debye Dispersion

The Debye model is represented by

```text
εr(ω) = ε∞ + (εs - ε∞)/(1 + jωτ)
```

The time-domain implementation uses an auxiliary polarization `P` with trapezoidal coefficients. This allows the electric update to include dispersive polarization without storing the full time history.

### Drude Dispersion

The Drude model is represented by

```text
εr(ω) = ε∞ - ωp²/(ω² - jγω)
```

The implementation evolves an auxiliary polarization current `J`, damped by `γ` and driven by `ωp² E`.

### Lorentz Dispersion

The Lorentz model is represented by

```text
εr(ω) = ε∞ + Δε ω0²/(ω0² - ω² + jγω)
```

The implementation evolves both polarization `P` and polarization current `J`.

## Scattering Diagnostics

For TF/SF simulations, the incident field at a monitor inside the total-field region is computed as a delayed copy of the configured source waveform:

```text
E_inc(x,t) = source(t - (x - x_TFSF_start)/v_inc)
```

If the reflected monitor is in the left scattered-field region, its incident field is zero, so the measured trace is already the scattered/reflected field.

Frequency-domain reflection and transmission are computed from FFT samples:

```text
R(f) = |E_ref(f)/E_inc,port1(f)|²
T(f) = |E_trn(f)/E_inc,port2(f)|²
```

Optional time gating extracts a signal segment before applying the FFT window.

## S-Parameters

The S-parameter workflow computes complex ratios:

```text
S11 = E_ref,1 / E_inc,1
S21 = E_trn,2 / E_inc,2
```

A mirrored reverse-incidence simulation is run internally to obtain `S22` and `S12`.

Reference-plane de-embedding applies phase shifts:

```text
S11' = S11 exp(j 2β Δx1)
S22' = S22 exp(j 2β Δx2)
S21' = S21 exp(j β (Δx1 + Δx2))
S12' = S12 exp(j β (Δx1 + Δx2))
```

where `β = 2πf/v_inc`.

Touchstone export supports `RI`, `MA`, and `DB` formats and includes metadata comments describing monitors, reference planes, grid settings, and de-embedding shifts.

## Analytical Comparisons

The code includes analytical normal-incidence comparisons for:

- dielectric slabs using standard transfer-matrix/Fabry-Perot formulas
- Debye, Drude, and Lorentz dispersive slabs
- dielectric half-space reflection through the Fresnel coefficient

For a single slab between identical incident and exit media:

```text
r = (r12 + r23 exp(-2jβd)) / (1 + r12 r23 exp(-2jβd))
t = (t12 t23 exp(-jβd)) / (1 + r12 r23 exp(-2jβd))
R = |r|²
T = Re(η1/η3) |t|²
```

