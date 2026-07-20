# MPCD Polymer Simulation

A small Python/NumPy project for building and testing a basic MPCD solvent and a coarse-grained polymer coupled to it.

## Current scope

Implemented / in progress:

- MPCD solvent particles
- periodic boundary conditions
- cell assignment
- streaming step
- collision step with random rotations
- basic polymer object
- harmonic bond forces
- velocity Verlet integration for polymer MD
- basic tests for conservation laws and stability

## Project structure

```text
mpcd_project/
├── src/mpcd/
│   ├── system.py          # MPCD solvent system
│   ├── polymer.py         # polymer data and initialization
│   ├── forces.py          # bond, bending, nonbonded forces
│   ├── md.py              # velocity Verlet integrator
│   ├── mpcd.py            # streaming, cells, collision
│   ├── coupling.py        # solvent-polymer MPCD coupling
│   ├── observables.py     # energy, momentum, bond lengths, etc.
│   └── visualization.py   # plots for debugging
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

## Basic usage

Run tests:

```bash
pytest
```

Run a specific test file:

```bash
pytest tests/test_polymer_forces.py
```

Run one test:

```bash
pytest tests/test_polymer_forces.py::test_total_bond_force_zero
```


## Recommended test parameters

```python
N_MONOMERS = 20
BOX = np.array([10, 10, 10])
MASS = 1.0
KBT = 1.0
BOND_LENGTH = 0.5
K_BOND = 10.0
DT = 0.001
N_STEPS = 10000
SEED = 12345
```

## Progress checklist up to Chapter 5

### 1. Basic MPCD solvent

- [*] Create `System` class
- [*] Initialize solvent positions uniformly in the box
- [*] Initialize solvent velocities from a Gaussian distribution
- [*] Implement periodic boundary conditions
- [*] Implement streaming step
- [*] Implement random grid shift
- [*] Implement cell assignment
- [*] Implement random rotation matrix
- [*] Implement MPCD collision step
- [*] Test rotation matrix determinant and orthogonality
- [*] Test momentum conservation in one cell
- [*] Test kinetic energy conservation in one cell
- [*] Test full solvent momentum conservation
- [*] Check velocity distribution 
- [*] Check cell occupancy distribution

### 2. Polymer MD without solvent

- [*] Create `Polymer` class
- [*] Generate straight-chain initial positions
- [ ] Generate random-walk initial positions
- [ ] Initialize polymer velocities
- [ ] Remove polymer center-of-mass velocity
- [*] Implement minimum image convention
- [*] Implement harmonic bond forces
- [*] Implement velocity Verlet integrator
- [*] Test initial bond lengths
- [ ] Test zero force at equilibrium bond length
- [ ] Test total bond force equals zero
- [ ] Test polymer momentum conservation
- [ ] Test polymer energy conservation
- [*] Test bond lengths remain stable
- [*] Plot polymer in 3D
- [*] Plot bond length evolution
- [ ] Plot total energy evolution

### 3. MPCD-polymer coupling

- [ ] Assign polymer beads to MPCD collision cells
- [ ] Assign solvent particles and polymer beads using the same shifted grid
- [ ] Compute mass-weighted cell center-of-mass velocity
- [ ] Rotate solvent relative velocities
- [ ] Rotate polymer bead relative velocities
- [ ] Test total momentum conservation in coupled collision
- [ ] Test kinetic energy conservation in coupled collision
- [ ] Run solvent + polymer simulation
- [ ] Track polymer bond lengths during coupled simulation
- [ ] Track polymer center-of-mass motion

### 4. Polymer observables

- [ ] Implement end-to-end distance
- [ ] Implement radius of gyration
- [ ] Track polymer center of mass
- [ ] Track bond length statistics
- [ ] Track kinetic, potential, and total energy
- [ ] Save observables to file
- [ ] Plot observables after simulation

### 5. Extensions toward Chapter 5

- [ ] Add bending stiffness
- [ ] Test bending energy calculation
- [ ] Add bending forces
- [ ] Test stability with bending stiffness
- [ ] Add excluded-volume interaction
- [ ] Add neighbor list or cell list for nonbonded forces
- [ ] Add multiple polymers
- [ ] Add inter-polymer interactions
- [ ] Add temperature gradient if needed
- [ ] Compare flexible and semiflexible polymer behavior
- [ ] Document final simulation parameters

## Main validation goals

Before adding new physics, confirm:

- [*] solvent collision conserves momentum
- [ ] polymer bond forces conserve total momentum
- [ ] polymer MD approximately conserves total energy
- [ ] coupled MPCD-polymer collision conserves momentum
- [ ] bond lengths stay stable
- [ ] simulation does not produce NaN or infinite values