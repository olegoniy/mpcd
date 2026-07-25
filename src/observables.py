import numpy as np

def total_momentum(system):
    return system.m * system.v.sum(axis=0)

def system_kinetic(system):
    return system.m * np.sum(system.v**2)*0.5

def coupled_cell_momentum(solventIndicies, system, polymerIndicies=None, polymer=None):
    monomersMomentum = np.zeros(shape=3)
    solventMomentum = np.sum(system.v[solventIndicies], axis=0) * system.m
    if polymer is not None:
        monomersMomentum = np.sum(polymer.v[polymerIndicies], axis=0) * polymer.m
    return solventMomentum + monomersMomentum

def coupledCMVelocity(solventIndicies, system, polymerIndicies=None, polymer=None):
    cellMomentum = coupled_cell_momentum(solventIndicies, system, polymerIndicies, polymer)
    if polymer is not None:
        massInCell = len(solventIndicies)*system.m + len(polymerIndicies)*polymer.m
        return cellMomentum/massInCell
    return system.v[solventIndicies].mean(axis=0)

def coupled_cm_kinetic_energy(solventIndicies, system, polymerIndicies=None, polymer=None):
    if polymer is not None:
        massInCell = len(solventIndicies)*system.m + len(polymerIndicies)*polymer.m
    else: 
        massInCell = len(solventIndicies)*system.m

    return massInCell * np.sum(coupledCMVelocity(solventIndicies, system, polymerIndicies, polymer)**2) * 0.5

def bond_lengths(polymer):
    res = []
    for i in range(polymer.nMonomers - 1):
        res.append(polymer.distInBC(i, i+1))
    return np.array(res)

def polymer_kinetic_energy(polymer):
    return polymer.kineticEnergy()

def polymer_bond_energy(polymer):
    return polymer.potentialEnergy()

def polymer_total_energy(polymer):
    return polymer.totalEnergy()

