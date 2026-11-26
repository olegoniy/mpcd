import numpy as np
import mpcd_cpp as mpcd


def kinetic_energy(system):
    v = np.asarray(system.v)
    return 0.5 * system.m * np.sum(v * v)


def momentum(system):
    return system.m * np.sum(np.asarray(system.v), axis=0)


def test_uniform_positions_inside_box(system):
    r = np.asarray(system.r)

    assert np.all(r >= 0.0)
    assert np.all(r < system.box)


def test_remove_drift(system):
    p = momentum(system)

    assert np.allclose(p, 0.0, atol=1e-12)


def test_stream_keeps_particles_inside_box(system):
    mpcd.stream(system)

    r = np.asarray(system.r)

    assert np.all(r >= 0.0)
    assert np.all(r < system.box)


def test_stream_does_not_change_velocities(system):
    velocities_before = np.asarray(system.v).copy()

    mpcd.stream(system)

    assert np.allclose(system.v, velocities_before)


def test_collision_conserves_total_momentum(system):
    momentum_before = momentum(system)

    mpcd.collide(system)

    momentum_after = momentum(system)

    assert np.allclose(momentum_after, momentum_before, atol=1e-10)


def test_collision_conserves_kinetic_energy(system):
    energy_before = kinetic_energy(system)

    mpcd.collide(system)

    energy_after = kinetic_energy(system)

    assert np.isclose(
        energy_after,
        energy_before,
        rtol=1e-12,
        atol=1e-10
    )


def test_cell_list_contains_every_particle_once():
    box = np.array([4.0, 4.0, 4.0])

    positions = np.array([
        [0.1, 0.1, 0.1],
        [1.2, 0.3, 0.4],
        [2.8, 3.1, 0.2],
        [3.9, 3.9, 3.9],
    ])

    cells = mpcd.distributeToCells(positions, box, 1.0)

    indices = np.asarray(cells.indices)

    assert len(indices) == len(positions)
    assert sorted(indices.tolist()) == list(range(len(positions)))
    assert cells.offsets[0] == 0
    assert cells.offsets[-1] == len(positions)