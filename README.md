# MPCD Polymer Simulation

## Overview

MPCD implementation using combination of classical MD for polymer behavior and coarse-grained simulation for a solvent

## Motivation

MPCD is cutting corners on elements of simulation, where we have no need for exact atomic interactions and replaces them with stochastic collision proccess instead. This allows for better performance, while still giving ability to see atomic impact on polymer beads. 

## Current Features

Coupled system with:
- periodic boundary conditions
- Galilean invariance
- Stochastic rotation for collision
- harmonic polymer bonds
- validation plots and tests

## Validation

| Check | What it verifies | Expected behavior |
|---|---|---|
| Rotation matrix orthogonality | Collision is a true rotation | `R.T @ R ≈ I` |
| Rotation matrix determinant | No reflection or scaling is introduced | `det(R) ≈ 1` |
| Solvent momentum conservation | MPCD collision treats cell center-of-mass motion correctly | Total solvent momentum remains constant |
| Solvent kinetic energy conservation | Relative velocities are rotated without changing their magnitude | Solvent kinetic energy remains constant |
| Polymer force balance | Bond forces are internal forces | Sum of polymer forces is zero |
| Polymer momentum conservation | Velocity Verlet with internal forces does not create net drift | Polymer momentum remains constant |
| Polymer energy conservation | Polymer MD is numerically stable | Small bounded energy error |
| Coupled cell momentum conservation | Mass-weighted solvent-polymer COM velocity is correct | Cell momentum remains constant |
| Coupled COM velocity conservation | Coupled collision preserves cell center-of-mass motion | COM velocity remains constant |
| Coupled relative speed conservation | Rotation preserves velocities relative to the COM | Relative speeds remain unchanged |
| Full coupled momentum conservation | Solvent-polymer exchange does not create net momentum | Total momentum remains constant |
| Full coupled kinetic energy conservation | Coupled collision does not heat or cool the system | Total kinetic energy remains constant during collision |
| NVE total energy plot | Full simulation remains stable without thermostat | Total energy deviation remains small |


## Repository Structure
```text
mpcd_project/
├── notebooks/
│   ├── plots.ipynb
│   └── SimulationDiagnostics.ipynb
│
├── src/
│   ├── system.py          # solvent state and parameters
│   ├── polymer.py         # polymer state and initialization
│   ├── forces.py          # polymer force calculations
│   ├── md.py              # velocity Verlet integrator
│   ├── mpcd.py            # MPCD streaming, cells, collision
│   └── observables.py     # diagnostics and measured quantities
│
├── scripts/
│  └─── run_sim.py
│
├── tests/
│   ├── test_mpcd.py
│   ├── test_md.py
│   └── test_physics.py
│
└── README.md
```

## Usage

Run all tests:

```bash
pytest
```

Run one test file:

```bash
pytest tests/test_polymer_forces.py
```

Run one test:

```bash
pytest tests/test_polymer_forces.py::test_total_bond_force_zero
```

## Recommended test parameters

```python
    system = System(
        N=10_000,
        box=[10.0, 10.0, 10.0],
        a=1.0,
        h=0.1,
        m=1.0,
        kBT=1.0,
        alpha=2/3*np.pi,
        seed=12345,
    )

    polymer = Polymer(
        nMonomers=100,
        box=[10.0, 10.0, 10.0],
        bondLength=0.25,
        m=2.0,
        k=100.0,
        kBT=1.0,
        dt=0.0001,
        seed=54321,
    )
```

---

# Progress checklist by book structure

## Chapter 2 — Basic MPCD algorithm

### Solvent state

- [x] Create `System` class
- [x] Store positions `r`, velocities `v`, box size, particle mass, cell size, time step, temperature, and rotation angle
- [x] Initialize solvent positions uniformly in the box
- [x] Initialize solvent velocities from a Gaussian distribution
- [x] Remove solvent center-of-mass velocity

### Streaming step

- [x] Implement ballistic streaming
- [x] Apply periodic boundary conditions after streaming
- [x] Test that all particles remain inside the simulation box

### Collision cells

- [x] Implement cell assignment
- [x] Implement random grid shift
- [x] Use shifted positions only for cell assignment
- [x] Test cell assignment with hand-picked positions
- [x] Test that every particle is assigned exactly once

### Collision step

- [x] Implement random rotation matrix
- [x] Test rotation matrix orthogonality
- [x] Test rotation matrix determinant
- [x] Compute cell center-of-mass velocity
- [x] Rotate relative velocities in each cell
- [x] Test cell momentum conservation
- [x] Test relative kinetic energy conservation in one cell
- [x] Test full-system momentum conservation

### Solvent diagnostics

- [x] Check velocity distribution
- [x] Check cell occupancy distribution
- [ ] Save solvent observables
- [ ] Plot solvent observables from saved data

---

## Chapter 3 — Embedded particles, polymers, and boundaries

### Polymer object

- [x] Create `Polymer` class
- [x] Store monomer positions `r`, velocities `v`, forces `f`, mass, bond length, and bond stiffness
- [x] Generate straight-chain initial positions
- [x] Generate random-walk initial positions
- [x] Initialize polymer velocities
- [ ] Remove polymer center-of-mass velocity
- [x] Plot polymer in 3D
- [x] Plot polymer with monomer indices

### Periodic distances

- [x] Implement minimum image convention
- [x] Use minimum image for bond vectors
- [x] Test 1D boundary case
- [x] Test 3D boundary case

### Polymer forces

- [x] Implement harmonic bond forces
- [x] Test zero force at equilibrium bond length
- [x] Test total bond force equals zero
- [x] Test equal and opposite pair forces

### Polymer MD

- [x] Implement velocity Verlet integrator
- [x] Compute forces before the first Verlet step
- [x] Test that equilibrium polymer with zero velocity stays fixed
- [x] Test polymer momentum conservation
- [x] Test polymer energy conservation
- [x] Test that bond lengths remain stable
- [x] Plot average bond length over time
- [x] Plot total energy over time

### MPCD-polymer coupling

- [x] Assign polymer beads to MPCD collision cells
- [x] Use the same shifted grid for solvent particles and polymer beads
- [x] Compute mass-weighted cell center-of-mass velocity
- [x] Rotate solvent relative velocities
- [x] Rotate polymer bead relative velocities
- [x] Test coupled-cell momentum conservation
- [x] Test coupled-cell kinetic energy conservation
- [x] Run solvent + polymer simulation
- [x] Track polymer bond lengths during coupled simulation
- [ ] Track polymer center-of-mass motion

---

## Chapter 4 — Thermostats and temperature control

### Temperature observables

- [ ] Implement solvent kinetic temperature
- [ ] Implement polymer kinetic temperature
- [ ] Implement total kinetic temperature for coupled systems
- [ ] Account for removed center-of-mass degrees of freedom where needed

### Solvent thermostat

- [ ] Implement cell-level velocity scaling thermostat
- [ ] Implement stochastic cell-level thermostat if needed
- [ ] Test that target temperature is maintained
- [ ] Check velocity distribution after thermostatting

### Polymer thermal behavior

- [ ] Check polymer temperature during standalone MD
- [ ] Check polymer temperature during MPCD coupling
- [ ] Decide whether polymer needs an explicit thermostat or only MPCD coupling

### Validation

- [ ] Plot temperature versus simulation step
- [ ] Test that thermostat does not create net momentum drift
- [ ] Test stability for long runs

---

## Chapter 5 — Transport and dynamical observables

To be written
