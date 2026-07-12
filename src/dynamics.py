import numpy as np

def streaming(system):
    system.r += system.h * system.v
    system.r %= system.box

def distributeToCells(system):
    particle_cells = np.floor(system.r/system.a).astype(int)    #What cells particle are in
    n_cells = np.ceil(system.box/system.a).astype(int)          #Create "boxes"

    cells = np.empty((n_cells[0], n_cells[1], n_cells[2]), dtype=object)

    for index in np.ndindex(cells.shape):
        cells[index] = []

    for i, cell in enumerate(particle_cells):
        ix, iy, iz = cell
        cells[ix, iy, iz].append(i)

    return cells
    
def rotationInCell(system, cell, rotation):
    v_com = system.r[cell].mean(axis=0)
    for velocity in system.v[cell]:
        velocity += (rotation-np.eye(3))@(velocity - v_com)

def generateRotation(system):
    alpha = system.rng.uniform(0, 2*np.pi)
    phi = system.rng.uniform(0, 2*np.pi)
    theta = system.rng.uniform(-1, 1)
    Rx = np.sqrt(1-theta**2)*np.cos(phi)
    Ry = np.sqrt(1-theta**2)*np.sin(phi)
    Rz = theta
    c = np.cos(alpha)
    s = np.sin(alpha)
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

def collision(system):
    cells = distributeToCells(system)
    rotationMatrix = generateRotation(system)
    for ix, iy, iz in np.ndindex(cells.shape):
        cell = cells[ix, iy, iz]
        if(system.r[cell].any()):
            rotationInCell(system, cell, rotationMatrix)

    

#[np.ceil(system.box[0]/system.a)][np.ceil(system.box[0]/system.a)][np.ceil(system.box[0]/system.a)]