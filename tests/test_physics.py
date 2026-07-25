import sys
import numpy as np
sys.path.append("../src")

from system import System

from mpcd import generateRotation, distributeToCells, rotateInCell, collide, stream

from observables import total_momentum, system_kinetic

system = System(
    N = 10000, 
    box = [3, 3, 3], 
    a = 1.0, 
    h = 0.1, 
    m = 1.0,
    kBT = 1.0, 
    alpha=np.pi/6, 
    seed=12345
    )


def test_local_momentum_after_collide():
    cells = distributeToCells(system.r, system.box, system.a)
    cell = cells[1,1,1]
    startMom = np.sum(system.v[cell]) * system.m
    for i in range(100):
        R = generateRotation(system)
        rotateInCell(system, cell, R)
        endMom = np.sum(system.v[cell]) * system.m
        np.testing.assert_allclose(
            startMom,
            endMom,
            rtol=0.0,
            atol=1e-13
        )
        startMom = endMom


def test_global_momentum_after_collide():
    start = total_momentum(system)
    for i in range(1000):
        collide(system)
        end = total_momentum(system)
        np.testing.assert_allclose(
                end,
                start,
                rtol=0.0,
                atol=1e-10
            )


def test_global_momentum_after_dynamic():
    start = total_momentum(system)
    for i in range(1000):
        stream(system)
        collide(system)
        end = total_momentum(system)
        np.testing.assert_allclose(
                end,
                start,
                rtol=0.0,
                atol=1e-10
            )
        
def test_kinetic_energy_after_dynamic():
    start = system_kinetic(system)
    for i in range(1000):
        stream(system)
        collide(system)
        end = system_kinetic(system)
        np.testing.assert_allclose(
                end,
                start,
                rtol=0.0,
                atol=1e-9
            )