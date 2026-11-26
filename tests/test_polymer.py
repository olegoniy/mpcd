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

def test_bond_lengths_known_configuration():
    polymer = make_polymer(3)

    polymer.r[:] = np.array([
        [4.0, 5.0, 5.0],
        [4.25, 5.0, 5.0],
        [4.50, 5.0, 5.0],
    ])

    lengths = np.asarray(
        mpcd.polymerBondLengths(polymer)
    )

    assert np.allclose(
        lengths,
        [0.25, 0.25],
        atol=1e-12
    )

def test_polymer_has_zero_spring_energy_at_rest_length():
    polymer = make_polymer(3)

    polymer.r[:] = np.array([
        [4.0, 5.0, 5.0],
        [4.25, 5.0, 5.0],
        [4.50, 5.0, 5.0],
    ])

    energy = mpcd.polymerPotentialEnergy(polymer)

    assert np.isclose(
        energy,
        0.0,
        atol=1e-12
    )

def test_bond_vectors_match_bond_lengths(polymer):
    polymer.initPositionsRandomWalk()

    vectors = np.asarray(
        mpcd.polymerBondVectors(polymer)
    )

    lengths = np.asarray(
        mpcd.polymerBondLengths(polymer)
    )

    vector_lengths = np.linalg.norm(vectors, axis=1)

    assert np.allclose(vector_lengths, lengths)


def test_distance_pbc_across_boundary(polymer):
    polymer.r[0] = [9.9, 5.0, 5.0]
    polymer.r[1] = [0.1, 5.0, 5.0]

    distance = polymer.distancePBC(0, 1)

    assert np.isclose(distance, 0.2)


def test_vector_pbc_across_boundary(polymer):
    polymer.r[0] = [9.9, 5.0, 5.0]
    polymer.r[1] = [0.1, 5.0, 5.0]

    vector = np.asarray(
        polymer.vecDiffPBC(0, 1)
    )

    assert np.isclose(np.linalg.norm(vector), 0.2)
    assert np.isclose(abs(vector[0]), 0.2)
    assert np.isclose(vector[1], 0.0)
    assert np.isclose(vector[2], 0.0)