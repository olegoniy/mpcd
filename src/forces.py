import numpy as np

def harmonicSpringPotential(k, len, len_0):
    return k*0.5*(len - len_0)**2

def bondForce(k, len, len_0):
    return k*(len - len_0)

def addBondforces(polymer):
    for i in range(polymer.nMonomers-1):
        len = polymer.distInBC(i, i+1)
        direction = polymer.vecDiffInBC(i, i+1)/len
        f = direction * bondForce(polymer.k, len, polymer.bondLength)
        polymer.f[i] += f   
        polymer.f[i+1] -= f

def forcesOnPolymer(polymer):
    polymer.f[:] = 0
    addBondforces(polymer)