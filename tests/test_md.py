import numpy as np
import sys
import numpy as np
sys.path.append("../src")

from polymer import Polymer
from md import verlet_step
from forces import forcesOnPolymer

polymer = Polymer(
    nMonomers = 1000,
    bondLength = 0.5,
    m = 10,
    k = 10,
    box = [10, 10, 10],
    kBT = 1,
    seed = 988
)

def average_length():
    res = 0
    for i in range(polymer.nMonomers-1):
        res += polymer.distInBC(i, i+1)
    return res/polymer.nMonomers

def test_initial_length():
    for i in range(polymer.nMonomers-1):  
        np.testing.assert_allclose(
            polymer.distInBC(i, i+1),
            polymer.bondLength,
            rtol=0.05
        )

def test_zero_force():
    forcesOnPolymer(polymer)
    np.testing.assert_allclose(
        polymer.f.sum(axis=0),
        np.zeros(3),
        rtol=0.0,
        atol=1e-14                  # actually better, than I thaught, I was expecting 1e-12
    )


def test_momentum_conservation():
    initial_momentum = polymer.m*polymer.v.sum(axis=0)
    for step in range(10000):
        verlet_step(polymer, 0.0001, forcesOnPolymer)
    np.testing.assert_allclose(
        polymer.m*polymer.v.sum(axis=0),
        initial_momentum,
        rtol=0.0,
        atol=1e-12
    )

def test_energy_conservation():
    initial_energy = polymer.totalEnergy()
    for step in range(10000):
        verlet_step(polymer, 0.0001, forcesOnPolymer)
    final_energy = polymer.totalEnergy()
    np.testing.assert_allclose(
        final_energy,
        initial_energy,
        rtol=1e-4
)