import numpy as np


def test_wrap_positions_positive_crossing(polymer):
    polymer.r[0] = [10.2, 5.0, 5.0]

    polymer.wrapPositions()

    assert np.allclose(
        polymer.r[0],
        [0.2, 5.0, 5.0]
    )

    assert np.array_equal(
        polymer.images[0],
        [1, 0, 0]
    )


def test_wrap_positions_negative_crossing(polymer):
    polymer.r[0] = [-0.3, 5.0, 5.0]

    polymer.wrapPositions()

    assert np.allclose(
        polymer.r[0],
        [9.7, 5.0, 5.0]
    )

    assert np.array_equal(
        polymer.images[0],
        [-1, 0, 0]
    )


def test_wrap_positions_multiple_boxes(polymer):
    polymer.r[0] = [25.3, -20.1, 30.4]

    polymer.wrapPositions()

    assert np.allclose(
        polymer.r[0],
        [5.3, 9.9, 0.4]
    )

    assert np.array_equal(
        polymer.images[0],
        [2, -3, 3]
    )


def test_wrapping_twice_does_not_change_image(polymer):
    polymer.r[0] = [10.2, 5.0, 5.0]

    polymer.wrapPositions()

    first_images = np.asarray(
        polymer.images
    ).copy()

    polymer.wrapPositions()

    assert np.array_equal(
        polymer.images,
        first_images
    )


def test_unwrapped_positions_recover_original_coordinates(polymer):
    original = np.array([
        [10.2, -0.3, 20.1],
        [5.0, 5.0, 5.0]
    ])

    polymer.r[:2] = original

    polymer.wrapPositions()

    unwrapped = np.asarray(
        polymer.unwrappedPositions()
    )

    assert np.allclose(
        unwrapped[:2],
        original
    )