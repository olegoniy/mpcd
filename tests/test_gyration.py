import numpy as np
import mpcd_cpp as mpcd


def make_polymer(n):
    box = np.array([10.0, 10.0, 10.0])

    return mpcd.Polymer(
        n,
        box,
        0.0001,
        0.25,
        1.0,
        100.0,
        0.75,
        42
    )

def test_center_of_mass_simple_configuration():
    polymer = make_polymer(3)

    polymer.r[:] = np.array([
        [4.0, 5.0, 5.0],
        [5.0, 5.0, 5.0],
        [6.0, 5.0, 5.0],
    ])

    center = np.asarray(
        mpcd.polymerCenterOfMass(polymer)
    )

    assert np.allclose(
        center,
        [5.0, 5.0, 5.0]
    )


def test_relative_positions_have_zero_mean(polymer):
    polymer.initPositionsRandomWalk()

    relative = np.asarray(
        mpcd.polymerPositionsRelativeToCenterOfMass(polymer)
    )

    assert np.allclose(
        relative.mean(axis=0),
        0.0,
        atol=1e-12
    )


def test_gyration_tensor_two_particle_line():
    polymer = make_polymer(2)

    polymer.r[:] = np.array([
        [4.0, 5.0, 5.0],
        [6.0, 5.0, 5.0]
    ])

    tensor = np.asarray(
        mpcd.polymerGyrationTensor(polymer)
    )

    expected = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0]
    ])

    assert np.allclose(tensor, expected)


def test_radius_of_gyration_two_particle_line():
    polymer = make_polymer(2)

    polymer.r[:] = np.array([
        [4.0, 5.0, 5.0],
        [6.0, 5.0, 5.0]
    ])

    radius = mpcd.polymerRadiusOfGyration(polymer)

    assert np.isclose(radius, 1.0)

def test_gyration_tensor_is_symmetric():
    polymer = make_polymer(20)

    polymer.initPositionsRandomWalk()

    tensor = np.asarray(
        mpcd.polymerGyrationTensor(polymer)
    )

    assert np.allclose(
        tensor,
        tensor.T,
        atol=1e-12
    )

def test_gyration_tensor_has_nonnegative_eigenvalues():
    polymer = make_polymer(20)

    polymer.initPositionsRandomWalk()

    tensor = np.asarray(
        mpcd.polymerGyrationTensor(polymer)
    )

    eigenvalues = np.linalg.eigvalsh(tensor)

    assert np.all(eigenvalues >= -1e-12)


def test_radius_of_gyration_matches_tensor_trace():
    polymer = make_polymer(20)

    polymer.initPositionsRandomWalk()

    tensor = np.asarray(
        mpcd.polymerGyrationTensor(polymer)
    )

    radius = mpcd.polymerRadiusOfGyration(polymer)

    assert np.isclose(
        radius**2,
        np.trace(tensor),
        rtol=1e-12,
        atol=1e-12
    )


def test_gyration_tensor_is_translation_invariant():
    polymer1 = make_polymer(4)
    polymer2 = make_polymer(4)

    positions = np.array([
        [1.0, 1.0, 1.0],
        [2.0, 1.0, 1.0],
        [2.0, 2.0, 1.0],
        [2.0, 2.0, 2.0],
    ])

    translation = np.array([3.0, 2.0, 1.0])

    polymer1.r[:] = positions
    polymer2.r[:] = positions + translation

    tensor1 = np.asarray(
        mpcd.polymerGyrationTensor(polymer1)
    )

    tensor2 = np.asarray(
        mpcd.polymerGyrationTensor(polymer2)
    )

    assert np.allclose(
        tensor1,
        tensor2,
        atol=1e-12
    )


def test_gyration_tensor_correct_across_periodic_boundary():
    polymer = make_polymer(2)

    # Physical coordinates:
    # 9.8 and 10.2
    polymer.r[:] = np.array([
        [9.8, 5.0, 5.0],
        [10.2, 5.0, 5.0]
    ])

    polymer.wrapPositions()

    assert np.allclose(
        polymer.r[:, 0],
        [9.8, 0.2]
    )

    assert np.array_equal(
        polymer.images[:, 0],
        [0, 1]
    )

    tensor = np.asarray(
        mpcd.polymerGyrationTensor(polymer)
    )

    expected = np.array([
        [0.04, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0]
    ])

    assert np.allclose(
        tensor,
        expected,
        atol=1e-12
    )

    assert np.isclose(
        mpcd.polymerRadiusOfGyration(polymer),
        0.2
    )


def polymer_momentum(polymer):
    return polymer.m * np.sum(
        np.asarray(polymer.v),
        axis=0
    )


def test_isolated_polymer_conserves_momentum(polymer):
    polymer.initPositionsRandomWalk()
    polymer.initVelocitiesNormal()
    polymer.removeDrift()

    momentum_before = polymer_momentum(polymer)

    mpcd.runPolymer(
        polymer,
        1000,
        100,
        0,
        0
    )

    momentum_after = polymer_momentum(polymer)

    assert np.allclose(
        momentum_after,
        momentum_before,
        atol=1e-10
    )

def test_isolated_polymer_energy_is_stable(polymer):
    polymer.initPositionsRandomWalk()
    polymer.initVelocitiesNormal()
    polymer.removeDrift()

    samples = mpcd.runPolymer(
        polymer,
        5000,
        50,
        0,
        0
    )

    energies = np.asarray(samples.totalEnergy)

    relative_spread = (
        energies.max() - energies.min()
    ) / abs(energies.mean())

    assert relative_spread < 1e-3