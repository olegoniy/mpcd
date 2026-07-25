import numpy as np
from observables import coupledCMVelocity

def stream(system):
    system.r += system.h * system.v
    system.r %= system.box

def distributeToCellsSolvent(system):
    particle_cells = np.floor(system.r/system.a).astype(int)    #What cells particle are in
    n_cells = np.ceil(system.box/system.a).astype(int)          #Create "boxes"

    cells = np.empty((n_cells[0], n_cells[1], n_cells[2]), dtype=object)

    for index in np.ndindex(cells.shape):
        cells[index] = []

    for i, cell in enumerate(particle_cells):
        ix, iy, iz = cell
        cells[ix, iy, iz].append(i)

    return cells
    
def rotateInCell(system, cell, rotation):
    cell = np.array(cell)

    v_com = system.v[cell].mean(axis=0)

    dv = system.v[cell] - v_com

    system.v[cell] = v_com + dv @ rotation.T

def generateRotation(system):
    phi = system.rng.uniform(0, 2*np.pi)
    theta = system.rng.uniform(-1, 1)
    Rx = np.sqrt(1-theta**2)*np.cos(phi)
    Ry = np.sqrt(1-theta**2)*np.sin(phi)
    Rz = theta
    c = np.cos(system.alpha)
    s = np.sin(system.alpha)
    return np.array([
        [
            Rx**2 + (1 - Rx**2) * c,
            Rx * Ry * (1 - c) - Rz * s,
            Rx * Rz * (1 - c) + Ry * s
        ],
        [
            Rx * Ry * (1 - c) + Rz * s,
            Ry**2 + (1 - Ry**2) * c,
            Ry * Rz * (1 - c) - Rx * s
        ],
        [
            Rx * Rz * (1 - c) - Ry * s,
            Ry * Rz * (1 - c) + Rx * s,
            Rz**2 + (1 - Rz**2) * c
        ]
    ])

def collide(system):
    
    shift = system.rng.uniform(-system.a / 2, system.a / 2, size=3)
    shifted_positions = (system.r + shift) % system.box
    cells = distributeToCells(shifted_positions, system.box, system.a)
    for ix, iy, iz in np.ndindex(cells.shape):
        cell = cells[ix, iy, iz]
        rotationMatrix = generateRotation(system)
        if len(cell) > 1:
            rotateInCell(system, cell, rotationMatrix)

def rotateCoupledCell(rotation, solventIndices, system, polymerIndices = None, polymer = None ):
    v_com = coupledCMVelocity(solventIndices, system, polymerIndices, polymer)
    dv_solvent = system.v[solventIndices] - v_com

    system.v[solventIndices] = v_com + dv_solvent @ rotation.T

    if polymer is not None:
        dv_monomers = polymer.v[polymerIndices] - v_com
        polymer.v[polymerIndices] = v_com + dv_monomers @ rotation.T

def distributeToCells(positions, box, a):
    n_cells = np.ceil(box / a).astype(int)
    particle_cells = np.floor(positions / a).astype(int)

    cells = np.empty((n_cells[0], n_cells[1], n_cells[2]), dtype=object)

    for index in np.ndindex(cells.shape):
        cells[index] = []
    
    for i, cell in enumerate(particle_cells):
        ix, iy, iz = cell
        cells[ix, iy, iz].append(i)

    return cells

def collideCoupled(system, polymer):
    shift = system.rng.uniform(-system.a / 2, system.a / 2, size=3)

    shifted_solvent = (system.r + shift) % system.box
    shifted_polymer = (polymer.r + shift) % polymer.box

    cells_solvent = distributeToCells(shifted_solvent, system.box, system.a)
    cells_polymer = distributeToCells(shifted_polymer, polymer.box, system.a)

    for index in np.ndindex(cells_solvent.shape):
        solvent_indices = cells_solvent[index]
        polymer_indices = cells_polymer[index]

        if len(solvent_indices) + len(polymer_indices)> 1:
            rotation = generateRotation(system)

            rotateCoupledCell(
                rotation,
                solvent_indices,
                system,
                polymer_indices,
                polymer
            )
    