import numpy as np

def system_momentum(system):
    return system.m * system.v.sum(axis=0)

def system_kinetic(system):
    return system.m * np.sum(system.v**2)*0.5

def coupledCellMomentum(solventIndicies, system, polymerIndicies=None, polymer=None):
    monomersMomentum = np.zeros(shape=3)
    solventMomentum = np.sum(system.v[solventIndicies], axis=0) * system.m
    if polymer is not None:
        monomersMomentum = np.sum(polymer.v[polymerIndicies], axis=0) * polymer.m
    return solventMomentum + monomersMomentum

def coupledCMVelocity(solventIndicies, system, polymerIndicies=None, polymer=None):
    cellMomentum = coupledCellMomentum(solventIndicies, system, polymerIndicies, polymer)
    if polymer is not None:
        massInCell = len(solventIndicies)*system.m + len(polymerIndicies)*polymer.m
        return cellMomentum/massInCell
    return system.v[solventIndicies].mean(axis=0)