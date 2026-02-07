# PFLOTRAN MPHASE CO2 Injection Simulation (mphase case)

Date: 2026-02-03

## Purpose (what this model demonstrates)
This case simulates **large-scale CO2 injection into a brine-filled formation** with an overlying **low-permeability caprock**, using PFLOTRAN’s **multiphase (MPHASE) CO2–H2O flow** model. The goal is to illustrate pressure response, plume migration, phase appearance/disappearance, and buoyant segregation in a simplified layered system.

## Toolchain / files
- Input deck: [co2/mphase/pflotran.in](co2/mphase/pflotran.in)
- Primary results: [co2/mphase/pflotran.h5](co2/mphase/pflotran.h5)
- Text output: [co2/mphase/pflotran.out](co2/mphase/pflotran.out)
- Reference slides/derivation (provided with the case): [co2/mphase/doc/co2.tex](co2/mphase/doc/co2.tex) and [co2/mphase/doc/co2-schematic.pdf](co2/mphase/doc/co2-schematic.pdf)

## Physics and governing equations (high level)
**Process model:** `SUBSURFACE_FLOW` with `MODE MPHASE`.

**Phases and components** (as documented in [co2/mphase/doc/co2.tex](co2/mphase/doc/co2.tex)):
- **Phases:** aqueous water (liquid) and supercritical CO2 (gas/sc phase)
- **Components:** H2O and CO2

The model solves coupled, non-isothermal multiphase flow with component transport:
- **Mass conservation** for each component (H2O, CO2), including advection in each phase and diffusive fluxes.
- **Darcy flow** for each phase using relative permeability and viscosity.
- **Energy conservation** including sensible heat in fluids + rock and conduction.

**Constitutive/EOS pieces** (case documentation):
- Supercritical CO2 EOS (Span–Wagner)
- Brine/mixture density and CO2 solubility (Duan-style correlations)
- Capillary pressure / relative permeability via Van Genuchten–Mualem

PFLOTRAN confirms the MPHASE unknown set and DOFs:
- “number of dofs = 3, number of phases = 2; mode = MPH: p, T, s/X” (see [co2/mphase/pflotran.out](co2/mphase/pflotran.out)).

## Model domain and discretization
**Domain size (structured grid):**
- Extent: **2500 m × 2500 m × 1000 m** (X × Y × Z)
- Origin: (0, 0, 0)
- Bounds: (0, 0, 0) to (2500, 2500, 1000)

**Grid resolution:**
- `NXYZ 21 21 20` → **8820 cells** total
- Uniform cell volumes (PFLOTRAN reports min=max volume)

## Stratigraphy / materials
The model uses four vertical layers (all laterally continuous):

| Unit | Z range [m] | Material ID | Porosity | Permeability (isotropic) |
|---|---:|---:|---:|---:|
| Basement | 0–400 | 4 | 0.10 | 1e-16 m² |
| Formation (injection interval) | 400–500 | 1 | 0.15 | 2e-14 m² |
| Caprock | 500–600 | 2 | 0.01 | 1e-17 m² |
| Overburden | 600–1000 | 3 | 0.15 | 2e-14 m² |

All materials use the same saturation function `sf2` (Van Genuchten + Mualem) with:
- Residual liquid saturation: 0.1
- Residual gas saturation: 0.0
- LAMBDA ($\lambda$): 0.762
- $\alpha$: 7.5e-4
- Max capillary pressure: 1e6 Pa

## Initial and boundary conditions
### Initial condition (entire domain)
Defined by `FLOW_CONDITION initial` applied to `REGION all`.

Key settings:
- **Hydrostatic liquid pressure** with a **datum at z = 1000 m**
- Datum pressure: **2.0e7 Pa (20 MPa)**
- **Temperature gradient:** 0.025 °C/m (negative sign in input indicates increasing T with depth)
- Reference temperature at the datum: **75 °C**
- Initial CO2 “concentration” is set extremely small (1e-8; used by PFLOTRAN for initialization)

### Boundary conditions
All six outer faces (top, bottom, north, south, east, west) are assigned the **same `initial` flow condition** as a Dirichlet-type boundary. In practice, this behaves like a large surrounding aquifer with hydrostatic pressure/temperature conditions.

## Sources/sinks (CO2 injection)
**Injection region:** a vertical line at the domain center:
- X = 1250 m, Y = 1250 m
- Z from 400 to 500 m (formation layer)

**Injection schedule (mass rate):**
- 0–20 years: **10 kg/s** (CO2 mass injection)
- after 20 years: **0 kg/s**

This is specified by `FLOW_CONDITION src` with a `RATE LIST` and applied via `SOURCE_SINK src`.

## Simulation time control and outputs
- Final time: **100 years**
- Initial timestep: **1e-8 years**
- Max timestep: **25 years**

**Output cadence:** 40 output times, with dense output around 20 years and then coarser to 100 years (see `OUTPUT TIMES` in [co2/mphase/pflotran.in](co2/mphase/pflotran.in)). The HDF5 actually contains **40 time groups** from **0 to 100 years**.

## Numerical methods / solver settings
From [co2/mphase/pflotran.in](co2/mphase/pflotran.in) and [co2/mphase/pflotran.out](co2/mphase/pflotran.out):
- Newton: ATOL 1e-12, RTOL 1e-8, max iters 25
- Time stepper acceleration enabled (`TS_ACCELERATION 8`)
- Linear solver: BCGS with block Jacobi preconditioner

The run completed in ~5 seconds wall-clock for 4 MPI ranks (as reported at the end of [co2/mphase/pflotran.out](co2/mphase/pflotran.out)).

## Key results (quantitative summary)
Values below are global statistics extracted from [co2/mphase/pflotran.h5](co2/mphase/pflotran.h5).

### Gas saturation (CO2 phase appearance / plume)
- At **t=0 y**: max Gas_Saturation = 0
- Peak max Gas_Saturation over time: **~0.7848** at **~50 y**
- At **t=100 y**: max Gas_Saturation = **~0.7826**; domain-average Gas_Saturation = **~0.031**

Interpretation: CO2 accumulates as a buoyant gas/sc phase in parts of the domain and persists after injection stops (20 y), with redistribution continuing toward later times.

### Pressure response
- Gas pressure (t=0 y): mean ~2.480e7 Pa
- Gas pressure (t=100 y): mean ~2.481e7 Pa
- Pressure change (100y - 0y):
  - min ~ -9.4e3 Pa
  - max ~ +1.17e5 Pa
  - mean ~ +4.35e3 Pa

Interpretation: injection generates a localized overpressure that relaxes/spreads after injection stops.

### Temperature
- Temperature is primarily controlled by the initial geothermal gradient and shows very small changes over time.

### CO2 composition (gas phase)
- Gas mole fraction of CO2 evolves from ~4.5e-7 initially to near 1.0 in CO2-rich regions by 100 years.

## Figures and animations (for slides)
We generated a **curated set of 6 variables** across **three sections** with:
- one animation per variable (GIF)
- five evenly spaced snapshots per variable (PNG)
- plus a last-time slice PNG

**Variables plotted (the “important 6”):**
- Gas_Saturation
- Liquid_Saturation
- Gas_Pressure [Pa]
- Liquid_Pressure [Pa]
- Gas_Mole_Fraction_CO2
- Temperature [C]

**Section outputs:**
- XY slices: [co2/mphase/figures_important_xy](co2/mphase/figures_important_xy)
- XZ slices: [co2/mphase/figures_important_xz](co2/mphase/figures_important_xz)
- YZ slices: [co2/mphase/figures_important_yz](co2/mphase/figures_important_yz)

Notes for presentation:
- Use **Gas_Saturation** animations to show plume growth (0–20 y) and post-injection migration (20–100 y).
- Use **Gas/Liquid pressure** to show overpressure near the injector and pressure diffusion.
- Use **Gas_Mole_Fraction_CO2** to highlight the CO2-rich gas region.

## Slide-ready outline (copy/paste)
### Slide 1 — Problem statement
- 3D CO2 injection into brine-filled formation with caprock
- PFLOTRAN `MPHASE` (non-isothermal, 2-phase, 2-component)
- Goal: plume migration + pressure response during/after injection

### Slide 2 — Domain + stratigraphy
- 2500×2500×1000 m; structured 21×21×20 (8820 cells)
- Layers: basement (0–400), formation (400–500), caprock (500–600), overburden (600–1000)
- Caprock permeability 1e-17 m² (low-K seal)

### Slide 3 — Initial & boundary conditions
- Hydrostatic initial pressure, datum at z=1000 m with 20 MPa
- Geothermal gradient 0.025 °C/m, T=75 °C at datum
- All outer faces held to the same hydrostatic/thermal condition

### Slide 4 — Injection schedule
- Vertical well at (1250,1250), z=400–500 m
- CO2 mass rate: 10 kg/s from 0–20 y, then shut-in
- Total simulation time 100 y

### Slide 5 — Results highlights
- Gas saturation plume appears and persists; peak max Sg ~0.785 (≈50 y)
- Overpressure localized; max ΔP ~1.17e5 Pa at 100 y relative to t0
- Temperature nearly unchanged; mainly geothermal structure

### Slide 6 — Visuals
- XY/XZ/YZ animations + 5 snapshots per animation
- Point audience to the three output folders and pick 1–2 key variables per slide

## Caveats / notes
- The current run uses `CO2_DATABASE /home/lal/software/pflotran/database/co2data0.dat` (absolute path) in the input deck; older relative paths in logs may reflect an earlier attempt.
- This case is flow-only (no geochemistry, no mineral trapping) and is intended as a conceptual/benchmark-style demonstration.
