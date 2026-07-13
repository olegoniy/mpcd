import numpy as np

def system_momentum(system):
    return system.m * system.v.sum(axis=0)

def system_kinetic(system):
    return system.m * np.sum(system.v**2)*0.5
