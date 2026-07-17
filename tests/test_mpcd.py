import sys
import numpy as np
sys.path.append("../src")

from system import System

from dynamics import generateRotation, distributeToCells, rotationInCell
from observables import system_momentum, system_kinetic

system = System(
    N = 10000, 
    box = [3, 3, 3], 
    a = 1.0, 
    h = 0.1, 
    m = 1.0,  
    kBT = 1.0, 
    alpha_deg=30, 
    seed=12345
    )


def test_is_orthogonal_matrix():
    
    EyeMrx = np.eye(3)
    for i in range(1000):
        R = generateRotation(system)
        np.testing.assert_allclose(
            R.T @ R,
            EyeMrx,
            rtol=0.0,
            atol=1e-15,
        )

def test_determinant():

    for i in range(100000):
        R = generateRotation(system)
        print(np.linalg.det(R))
        np.testing.assert_allclose(
            np.linalg.det(R),
            1,
            rtol=0.0, 
            atol=1e-15   
        )

def test_rotation():
    cells = distributeToCells(system)
    cell = cells[1,1,1]
    for i in range(10):
        R = generateRotation(system)
        start = system.v[cell]
        rotationInCell(system, cell, R)
        end = system.v[cell]
        assert not np.allclose(start, end)


