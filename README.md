# MPCD Polymer Simulation

Python/NumPy implementation of a basic MPCD solvent and a coarse-grained polymer model, developed step by step following the structure of the MPCD lecture notes/book.

## Scope

Current focus:

- basic MPCD solvent
- periodic boundary conditions
- collision cells and random grid shift
- stochastic rotation collision step
- coarse-grained polymer beads
- harmonic bond forces
- velocity Verlet integration
- validation tests for conservation laws and numerical stability

## Project structure

```text
mpcd_project/
├── src/mpcd/
│   ├── system.py          # solvent state and parameters
│   ├── polymer.py         # polymer state and initialization
│   ├── forces.py          # polymer force calculations
│   ├── md.py              # velocity Verlet integrator
│   ├── mpcd.py            # MPCD streaming, cells, collision
│   ├── coupling.py        # polymer-solvent coupling
│   ├── observables.py     # diagnostics and measured quantities
│   └── visualization.py   # debugging plots
│
├── scripts/
│   ├── run_solvent.py
│   ├── run_polymer_md.py
│   └── run_polymer_mpcd.py
│
├── tests/
│   ├── test_rotation.py
│   ├── test_cells.py
│   ├── test_polymer_forces.py
│   ├── test_md.py
│   └── test_coupling.py
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

Run scripts:

```bash
python scripts/run_solvent.py
python scripts/run_polymer_md.py
python scripts/run_polymer_mpcd.py
```

## Recommended test parameters

```python
nMonomers = 20
box = np.array([10.0, 10.0, 10.0])
m = 1.0
kBT = 1.0
bondLength = 0.5
k = 10.0
t = 0.001
Number of steps = 10000
SEED = <Any positive integerß>
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

- [ ] Assign polymer beads to MPCD collision cells
- [ ] Use the same shifted grid for solvent particles and polymer beads
- [ ] Compute mass-weighted cell center-of-mass velocity
- [ ] Rotate solvent relative velocities
- [ ] Rotate polymer bead relative velocities
- [ ] Test coupled-cell momentum conservation
- [ ] Test coupled-cell kinetic energy conservation
- [ ] Run solvent + polymer simulation
- [ ] Track polymer bond lengths during coupled simulation
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
