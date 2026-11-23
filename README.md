# MPCD — C++ Simulation Library with Python Research Interface

A research implementation of **Multiparticle Collision Dynamics (MPCD)** and coarse-grained polymer dynamics.

The project consists of two layers:

- a **C++ simulation library**, developed from scratch, containing the numerical core, physical models, observables, thermostatting, solvent–polymer coupling, and transport calculations;
- a **Python research layer**, used to configure simulations, run validation experiments, compare results with theory, and visualize the resulting physics.

The aim is to keep the simulation itself in C++ while using Python as a convenient interface for research and analysis.

---

## Reference model

The physical formulation follows:

**R. G. Winkler, “Flow Simulations with Multiparticle Collision Dynamics,” in _Hierarchical Methods for Dynamics in Complex Molecular Systems_, J. Grotendorst, G. Sutmann, G. Gompper, and D. Marx (eds.), Forschungszentrum Jülich, 2012, pp. 417–443.**

In particular:

- MPCD/SRD solvent dynamics and random grid shifting follow **Section 2**;
- solvent coupling to embedded particles / polymer monomers follows **Section 3**, especially Eq. (9);
- the cell-level canonical thermostat follows **Section 4**, especially Eqs. (14)–(17);
- VACF and self-diffusion analysis follow **Section 5.1**, Eqs. (18)–(26);
- the mean-free-path comparison follows **Figure 2** and the accompanying discussion.

This README focuses on the implementation and on project-specific choices rather than reproducing the derivations from the reference.

---

## Architecture

```text
                         MPCD project
                              │
              ┌───────────────┴───────────────┐
              │                               │
        C++ simulation library          Python research layer
              │                               │
        MPCD / polymer dynamics         simulation setup
        thermostatting                  validation
        observables                     parameter studies
        transport analysis              plotting / notebooks
              │                               │
              └──────── pybind11 ─────────────┘
```

The C++ core performs the numerical evolution and sampling. Python provides the high-level interface.

---

## C++ library

### MPCD solvent

The solvent implementation follows Section 2 of Winkler.

Implemented functionality includes:

- three-dimensional solvent particles;
- uniform position initialization;
- Maxwell-Boltzmann / Gaussian velocity initialization;
- center-of-mass drift removal;
- ballistic streaming;
- periodic boundary conditions;
- collision-cell assignment;
- random grid shifting;
- stochastic rotation collisions around randomly oriented axes;
- cell center-of-mass velocities;
- momentum-conserving collision dynamics.

### Polymer dynamics

Solvent–polymer collision coupling follows the embedded-particle approach in Section 3 of Winkler.

Project-specific choices:

- ideal bead-spring polymer;
- harmonic nearest-neighbour bonds;
- velocity-Verlet integration;
- separate polymer sampling and diagnostics.

### Canonical thermostat

The cell-level thermostat follows Section 4 of Winkler.

Project-specific implementation details:

- thermostatting is selectable at run time;
- the same interface is used for solvent and coupled solvent–polymer simulations;
- temperature, momentum, energy, and relative-velocity diagnostics are sampled directly in C++.

### Observables

The C++ library exposes observables for the solvent, polymer, and coupled system, including:

**Solvent**
- kinetic energy;
- temperature;
- total momentum;
- relative-speed diagnostics;
- velocity autocorrelation function;
- averaged VACF sample counts;
- self-diffusion coefficient.

**Polymer**
- kinetic, potential, and total energy;
- momentum;
- bond lengths and maximum bond length;
- end-to-end distance;
- bond vectors;
- sampled polymer configurations.

**Coupled system**
- solvent and polymer temperatures;
- total momentum;
- energy diagnostics.

---

## Transport analysis

Transport analysis is based on Section 5.1 of Winkler.

### VACF

The numerical VACF is compared with the analytical molecular-chaos result in Eqs. (23)–(24) and with the behaviour shown in Figure 2.

A project-specific addition is **multi-time-origin averaging**.

During one equilibrium run, the C++ `VacfAccumulator` periodically stores new velocity reference states and accumulates correlations for each lag. It also tracks the number of samples contributing to every lag and removes expired origins once they exceed the requested maximum lag.

This provides a statistically cleaner VACF without requiring separate simulation runs.

### Mean free path

The Figure 2 comparison is reproduced for representative values of

- \(\lambda/a = 0.1\);
- \(\lambda/a = 1.0\).

As in Winkler, \(\lambda/a\) is derived from the simulation parameters rather than stored as an independent system variable.

### Self-diffusion

The numerical diffusion coefficient is obtained from the Green-Kubo relation in Eq. (18) and compared with the molecular-chaos prediction in Eq. (26).

The implementation evaluates the Green-Kubo sum from the **averaged VACF** and uses a configurable maximum lag so that the numerical result can be checked for cutoff convergence.

---

## Validation

The implementation is tested against conservation laws, equilibrium statistics, and analytical predictions.

Current checks include:

- solvent momentum conservation;
- total-energy conservation in unthermostatted equilibrium runs;
- momentum conservation with the canonical thermostat;
- solvent and polymer temperature behaviour;
- Maxwell-Boltzmann velocity statistics;
- relaxation of anisotropic initial velocity distributions;
- isolated-polymer energy and bond stability;
- \(C_v(0)\) at equilibrium;
- VACF versus the molecular-chaos prediction;
- comparison of the \(\lambda/a=0.1\) and \(\lambda/a=1\) regimes;
- Green-Kubo self-diffusion versus the analytical result.

---

## Python interface

The library is exposed through **pybind11**.

A typical pure-solvent run looks like:

```python
import numpy as np
import mpcd_cpp as mpcd

box = np.array([10.0, 10.0, 10.0])

system = mpcd.System(
    5000,
    box,
    1.0,
    0.1,
    1.0,
    0.75,
    np.pi / 2,
    42
)

system.initPositionsUniform()
system.initVelocitiesNormal()

samples = mpcd.runSolvent(
    system=system,
    steps=5000,
    sample_period=100,
    vacf_max_lag=100,
    vacf_origin_stride=10
)
```

Sampled quantities are then available directly in Python:

```python
lags = np.asarray(samples.averagedVacfLags)
vacf = np.asarray(samples.averagedVelocityAutocorrelation)
counts = np.asarray(samples.averagedVacfCounts)

vacf_normalized = vacf / vacf[0]
```

The simulation and physical observables remain in C++; Python is used for configuration, comparison, and visualization.

---

## Building

### Requirements

- C++20 compiler
- CMake
- Eigen3
- pybind11
- Python 3
- NumPy

For analysis and testing, the project also uses Matplotlib and pytest.

### Build

```bash
git clone https://github.com/olegoniy/mpcd.git
cd mpcd

cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
```

The build produces the C++ core and the `mpcd_cpp` Python extension.

---

## Project structure

```text
cpp/
├── include/mpcd/       # public C++ interfaces
├── src/                # implementations
└── bindings/           # pybind11 bindings
```

The Python side contains notebooks and analysis code for simulation setup, validation, parameter studies, and plotting.

---

## Current directions

The next research steps are:

- compare the self-diffusion coefficient over several collision times / mean free paths;
- add self-avoidance to the polymer model;
- implement the gyration tensor and radius of gyration;
- compare rigidity-related observables for ideal and self-avoiding polymers;
- extend the polymer interaction model;
- study heat transport using spatially separated hot and cold thermostat regions and the resulting temperature gradient.

---

## Testing

Python-facing tests can be run with

```bash
pytest
```

Physics-level validation is additionally performed in notebooks through conservation checks, distributions, time series, and comparison with analytical predictions.

---

## Technologies

C++20 · Eigen · CMake · pybind11 · Python · NumPy · Matplotlib · pytest

---

## Reference

[1] **R. G. Winkler**, “Flow Simulations with Multiparticle Collision Dynamics,” in _Hierarchical Methods for Dynamics in Complex Molecular Systems_, J. Grotendorst, G. Sutmann, G. Gompper, and D. Marx (eds.), Forschungszentrum Jülich, 2012, pp. 417–443.

---

## Author

Developed by **Oleg Onyshchenko** as a computational / soft-matter physics research project.

Repository:  
https://github.com/olegoniy/mpcd
