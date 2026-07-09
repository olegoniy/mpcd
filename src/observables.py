import numpy as np

def system_impuls(system):
    return system.m * system.v.sum(axis=0)

def system_kinetic(system):
    return system.m * sum(system.v**2)*0.5
