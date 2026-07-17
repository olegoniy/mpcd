import numpy as np
import sys
import numpy as np
sys.path.append("../src")

from polymer import Polymer
from md import verlet_step
from forces import forcesOnPolymer

polymer = Polymer(
    nMonomers = 100,
    bondLength = 0.5,
    m = 10,
    k = 100,
    box = [10, 10, 10],
    kBT = 1,
    seed = 988
)

def average_length():
    res = 0
    for i in range(polymer.nMonomers-1):
        res += polymer.distInBC(i, i+1)
    return res/polymer.nMonomers

def test_av_lentgh_coservation():
    initialAverage = average_length()
    for step in range(10000):
        verlet_step(polymer, 0.001, forcesOnPolymer)
        new_av_length = average_length()
        np.testing.assert_allclose(
            initialAverage, 
            new_av_length,
            rtol=0.)

def test_momentum_conservation():
    initial_momentum = polymer.m*polymer.v.sum(axis=0)
    for step in range():
        pass