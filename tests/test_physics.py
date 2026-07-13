import sys
import numpy as np
sys.path.append("../src")

from system import System

from dynamics import generateRotation, distributeToCells, rotationInCell, collision, streaming

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


def test_lokal_momentum_after_collision():
    cells = distributeToCells(system)
    cell = cells[1,1,1]
    startMom = np.sum(system.v[cell]) * system.m
    for i in range(100):
        R = generateRotation(system)
        rotationInCell(system, cell, R)
        endMom = np.sum(system.v[cell]) * system.m
        np.testing.assert_allclose(
            startMom,
            endMom,
            rtol=0.0,
            atol=1e-13
        )
        startMom = endMom


def test_global_momentum_after_collision():
    start = system_momentum(system)
    for i in range(1000):
        collision(system)
        end = system_momentum(system)
        np.testing.assert_allclose(
                end,
                start,
                rtol=0.0,
                atol=1e-10
            )


def test_global_momentum_after_dynamic():
    start = system_momentum(system)
    for i in range(1000):
        streaming(system)
        collision(system)
        end = system_momentum(system)
        np.testing.assert_allclose(
                end,
                start,
                rtol=0.0,
                atol=1e-10
            )
        
def test_kinetic_energy_after_dynamic():
    start = system_kinetic(system)
    for i in range(1000):
        streaming(system)
        collision(system)
        end = system_kinetic(system)
        np.testing.assert_allclose(
                end,
                start,
                rtol=0.0,
                atol=1e-10
            )