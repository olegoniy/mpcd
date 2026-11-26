import numpy as np
import pytest

import sys
from pathlib import Path

sys.path.append(str(Path("../../build-release").resolve()))

import mpcd_cpp as mpcd


@pytest.fixture
def box():
    return np.array([10.0, 10.0, 10.0])


@pytest.fixture
def system(box):
    system = mpcd.System(
        5000,
        box,
        1.0,
        0.1,
        1.0,
        0.75,
        np.pi / 2,
        42
    )

    system.initPositionsUniform()
    system.initVelocitiesNormal()
    system.removeDrift()

    return system


@pytest.fixture
def polymer(box):
    polymer = mpcd.Polymer(
        20,
        box,
        0.0001,
        0.25,
        1.0,
        100.0,
        0.75,
        42
    )

    return polymer