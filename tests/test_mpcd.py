import sys
import numpy as np
import pytest
sys.path.append("../src")

from system import System
from polymer import Polymer

from mpcd import generateRotation, distributeToCells, rotateInCell, rotateCoupledCell, distributeToCellsSolvent, collideCoupled
from observables import * 

@pytest.fixture
def coupled_system():
    system = System(
        N=10_000,
        box=[10.0, 10.0, 10.0],
        a=1.0,
        h=0.1,
        m=1.0,
        kBT=1.0,
        alpha=2/3*np.pi,
        seed=12345,
    )

    polymer = Polymer(
        nMonomers=100,
        box=[10.0, 10.0, 10.0],
        bondLength=0.25,
        m=2.0,
        k=100.0,
        kBT=1.0,
        dt=0.0001,
        seed=54321,
    )

    return system, polymer

def test_is_orthogonal_matrix(coupled_system):
    system, _ = coupled_system
    
    EyeMrx = np.eye(3)
    for i in range(1000):
        R = generateRotation(system)
        np.testing.assert_allclose(
            R.T @ R,
            EyeMrx,
            rtol=0.0,
            atol=1e-15,
        )

def test_determinant(coupled_system):
    system, _ = coupled_system
    for i in range(1000):
        R = generateRotation(system)
        np.testing.assert_allclose(
            np.linalg.det(R),
            1,
            rtol=0.0, 
            atol=1e-14   
        )

def test_rotation(coupled_system):
    system, _ = coupled_system
    cells = distributeToCellsSolvent(system)
    cell = cells[1,1,1]
    for i in range(10):
        R = generateRotation(system)
        start = system.v[cell]
        rotateInCell(system, cell, R)
        end = system.v[cell]
        assert not np.allclose(start, end)

def test_rotate_coupled_cell_exact_result(coupled_system):

    system, polymer = coupled_system

    solventIndices = np.array([0, 1, 2, 3])
    polymerIndices = np.array([0, 1, 2])

    system.v[solventIndices] = np.array([
        [ 1.0,  0.0,  0.0],
        [-1.0,  2.0,  0.0],
        [ 0.5, -1.0,  1.0],
        [ 2.0,  0.5, -1.0],
    ])

    polymer.v[polymerIndices] = np.array([
        [ 0.0,  2.0,  0.5],
        [-2.0,  0.0,  1.0],
        [ 1.0, -1.0, -0.5],
    ])

    # 90° rotation around the z-axis
    rotation = np.array([
        [0.0, -1.0, 0.0],
        [1.0,  0.0, 0.0],
        [0.0,  0.0, 1.0],
    ])

    solvent_before = system.v[solventIndices].copy()
    polymer_before = polymer.v[polymerIndices].copy()

    total_momentum = (
        system.m * solvent_before.sum(axis=0)
        + polymer.m * polymer_before.sum(axis=0)
    )

    total_mass = (
        system.m * len(solventIndices)
        + polymer.m * len(polymerIndices)
    )

    expected_v_com = total_momentum / total_mass

    expected_solvent = (
        expected_v_com
        + (solvent_before - expected_v_com) @ rotation.T
    )

    expected_polymer = (
        expected_v_com
        + (polymer_before - expected_v_com) @ rotation.T
    )

    rotateCoupledCell(
        rotation,
        solventIndices,
        system,
        polymerIndices,
        polymer,
    )

    np.testing.assert_allclose(
        system.v[solventIndices],
        expected_solvent,
        rtol=0.0,
        atol=1e-13,
    )

    np.testing.assert_allclose(
        polymer.v[polymerIndices],
        expected_polymer,
        rtol=0.0,
        atol=1e-13,
    )

def test_rotate_coupled_cell_conserves_momentum(coupled_system):
    system, polymer = coupled_system

    solventIndices = np.array([0, 1, 2, 3, 4])
    polymerIndices = np.array([0, 1, 2])

    rotation = np.array([
        [0.0, -1.0, 0.0],
        [1.0,  0.0, 0.0],
        [0.0,  0.0, 1.0],
    ])

    momentum_before = coupled_cell_momentum(
        solventIndices,
        system,
        polymerIndices,
        polymer,
    ).copy()

    rotateCoupledCell(
        rotation,
        solventIndices,
        system,
        polymerIndices,
        polymer,
    )

    momentum_after = coupled_cell_momentum(
        solventIndices,
        system,
        polymerIndices,
        polymer,
    )

    np.testing.assert_allclose(
        momentum_after,
        momentum_before,
        rtol=1e-13,
        atol=1e-12,
    )

def test_rotate_coupled_cell_preserves_com_velocity(coupled_system):

    system, polymer = coupled_system

    solventIndices = np.array([5, 7, 11, 20])
    polymerIndices = np.array([3, 8, 12])

    rotation = np.array([
        [0.36, -0.48, 0.80],
        [0.80,  0.60, 0.00],
        [-0.48, 0.64, 0.60],
    ])

    v_com_before = coupledCMVelocity(
        solventIndices,
        system,
        polymerIndices,
        polymer,
    ).copy()

    rotateCoupledCell(
        rotation,
        solventIndices,
        system,
        polymerIndices,
        polymer,
    )

    v_com_after = coupledCMVelocity(
        solventIndices,
        system,
        polymerIndices,
        polymer,
    )

    np.testing.assert_allclose(
        v_com_after,
        v_com_before,
        rtol=1e-13,
        atol=1e-12,
    )

def test_rotate_coupled_cell_expected_velocities(
    coupled_system,
):
    system, polymer = coupled_system

    solventIndices = np.array([0, 1])
    polymerIndices = np.array([0, 1])

    system.v[solventIndices] = np.array([
        [1.0, 0.0, 0.0],
        [-1.0, 2.0, 0.0],
    ])

    polymer.v[polymerIndices] = np.array([
        [0.0, 1.0, 1.0],
        [2.0, -1.0, 0.0],
    ])

    rotation = np.array([
        [0.0, -1.0, 0.0],
        [1.0,  0.0, 0.0],
        [0.0,  0.0, 1.0],
    ])

    v_com = coupledCMVelocity(
        solventIndices,
        system,
        polymerIndices,
        polymer,
    ).copy()

    solvent_before = system.v[solventIndices].copy()
    polymer_before = polymer.v[polymerIndices].copy()

    expected_solvent = (
        v_com
        + (solvent_before - v_com) @ rotation.T
    )

    expected_polymer = (
        v_com
        + (polymer_before - v_com) @ rotation.T
    )

    rotateCoupledCell(
        rotation,
        solventIndices,
        system,
        polymerIndices,
        polymer,
    )

    np.testing.assert_allclose(
        system.v[solventIndices],
        expected_solvent,
        rtol=0.0,
        atol=1e-13,
    )

    np.testing.assert_allclose(
        polymer.v[polymerIndices],
        expected_polymer,
        rtol=0.0,
        atol=1e-13,
    )

def test_rotate_coupled_cell_preserves_relative_speeds(
    coupled_system,
):
    system, polymer = coupled_system

    solventIndices = np.array([0, 2, 4, 6])
    polymerIndices = np.array([1, 3, 5])

    rotation = np.array([
        [0.0, -1.0, 0.0],
        [1.0,  0.0, 0.0],
        [0.0,  0.0, 1.0],
    ])

    v_com = coupledCMVelocity(
        solventIndices,
        system,
        polymerIndices,
        polymer,
    ).copy()

    solvent_relative_before = (
        system.v[solventIndices] - v_com
    )
    polymer_relative_before = (
        polymer.v[polymerIndices] - v_com
    )

    solvent_speeds_before = np.linalg.norm(
        solvent_relative_before,
        axis=1,
    )
    polymer_speeds_before = np.linalg.norm(
        polymer_relative_before,
        axis=1,
    )

    rotateCoupledCell(
        rotation,
        solventIndices,
        system,
        polymerIndices,
        polymer,
    )

    solvent_speeds_after = np.linalg.norm(
        system.v[solventIndices] - v_com,
        axis=1,
    )
    polymer_speeds_after = np.linalg.norm(
        polymer.v[polymerIndices] - v_com,
        axis=1,
    )

    np.testing.assert_allclose(
        solvent_speeds_after,
        solvent_speeds_before,
        rtol=1e-13,
        atol=1e-12,
    )

    np.testing.assert_allclose(
        polymer_speeds_after,
        polymer_speeds_before,
        rtol=1e-13,
        atol=1e-12,
    )

def test_all_particles_assigned_once():
    rng = np.random.default_rng(12345)

    n_particles = 1000
    box = np.array([10.0, 10.0, 10.0])
    a = 1.0

    positions = rng.uniform(0, box, size=(n_particles, 3))

    cells = distributeToCells(positions, box, a)

    assigned = []

    for index in np.ndindex(cells.shape):
        assigned.extend(cells[index])

    assigned = np.array(assigned)

    assert len(assigned) == n_particles
    assert len(np.unique(assigned)) == n_particles
    assert set(assigned) == set(range(n_particles))

def test_coupled_cell_momentum_conservation(coupled_system):

    system, polymer = coupled_system

    solventIndices = np.array([0, 1, 2, 3])
    polymerIndices = np.array([0, 1])

    rotation = generateRotation(system)

    P_before = coupled_cell_momentum(
        solventIndices,
        system,
        polymerIndices,
        polymer
    )

    rotateCoupledCell(
        rotation,
        solventIndices,
        system,
        polymerIndices,
        polymer
    )

    P_after = coupled_cell_momentum(
        solventIndices,
        system,
        polymerIndices,
        polymer
    )

    np.testing.assert_allclose(P_after, P_before, rtol=0.0, atol=1e-10)

def test_coupled_momentum_conservation(coupled_system):
    system, polymer = coupled_system
    
    solventIndices = np.array([0, 1, 2, 3])
    polymerIndices = np.array([0, 1])

    rotation = generateRotation(system)

    P_before = coupled_cell_momentum(
        solventIndices,
        system,
        polymerIndices,
        polymer
    )
    for i in range(100):
        rotateCoupledCell(
            rotation,
            solventIndices,
            system,
            polymerIndices,
            polymer
        )

    P_after = coupled_cell_momentum(
        solventIndices,
        system,
        polymerIndices,
        polymer
    )

    np.testing.assert_allclose(P_after, P_before, rtol=0.0, atol=1e-14)

def test_coupled_cm_kinetic_energy_conservation(coupled_system):
    system, polymer = coupled_system

    solventIndicies = np.array([0, 1, 2, 3])
    polymerIndicies = np.array([0, 1])

    rotation = generateRotation(system)


    E_before = coupled_cm_kinetic_energy(
        solventIndicies, 
        system, 
        polymerIndicies, 
        polymer
    )

    for i in range(100):
        rotateCoupledCell(
            rotation,
            solventIndicies,
            system,
            polymerIndicies,
            polymer
        )

    E_after = coupled_cm_kinetic_energy(
        solventIndicies, 
        system, 
        polymerIndicies, 
        polymer
    )
    np.testing.assert_allclose(E_after, E_before, rtol=0.0, atol=1e-14)

def test_coupled_full_momentum_conservation(coupled_system):
    system, polymer = coupled_system

    p_full_before = total_momentum(system) + total_momentum(polymer)

    for i in range(100):
        collideCoupled(system, polymer)

    p_full_after = total_momentum(system) + total_momentum(polymer)

    np.testing.assert_allclose(p_full_after, p_full_before, rtol=0.0, atol=1e-12)

def test_coupled_full_kinetic_energy_conservation(coupled_system):
    system, polymer = coupled_system

    kinetic_before = system_kinetic(system) + polymer_kinetic_energy(polymer)

    for i in range(100):
            collideCoupled(system, polymer)

    kinetic_after = system_kinetic(system) + polymer_kinetic_energy(polymer)
    
    np.testing.assert_allclose(kinetic_after, kinetic_before, rtol=0.0, atol=1e-10)

