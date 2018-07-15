# -*- coding: utf-8 -*-

# Copyright (c) PyGLD Project Contributors
# https://github.com/jnsebgosselin/pygld
#
# This is part of PyGLD (Python Ground-Loop Designer).
# Licensed under the terms of the MIT License.


# ---- Standard imports

import os

# ---- Third party imports

import pytest
from numpy import pi

# ---- Local imports

from pygld import Pipe, PipeMaterial
from pygld.api.materials import PREDEFINED_MATERIALS


@pytest.fixture
def pipe():
    return Pipe(3.3, 4.2)


def test_init_pipe(pipe):
    """
    Test that Pipe is initialized correctly and that the default properties
    are set and calculated as expected.
    """
    assert pipe

    assert pipe.di == 3.3 and pipe.do == 4.2
    assert round(pipe.wt, 6) == 0.450000
    assert round(pipe.Ai, 6) == 8.552986
    assert round(pipe.Ao, 6) == 13.854424
    assert round(pipe.Rcond, 6) == 0.095955

    assert pipe.material._category == 'Pipe'
    assert pipe.material._material == 'HDPE'


def test_print_pipe(pipe):
    """
    Test that printing an instance of a Pipe is working as expected.
    """
    expected_content = (
        "Inner diameter: 3.30 cm\n"
        "Outer diameter: 4.20 cm\n"
        "Wall thickness: 0.45 cm\n"
        "Inner area: 8.55 cm²\n"
        "Outer area: 13.85 cm²\n"
        "Pipe material: HDPE\n"
        "Thermal conductivity: 0.400 W/m·k\n"
        "Volumetric heat capacity: 1500 J/m³·K\n"
        "Thermal diffusivity: 2.67e-04 m²/s\n"
        "Conductive thermal resistance: 0.09596 m·K/W")

    assert pipe.__str__() == expected_content


def test_set_material():
    """
    Test that the interface to set and get the pipe's material is working as
    expected.
    """
    # Assert that the material is set correctly from the optional argument.

    pipe = Pipe(3.3, 4.2, PipeMaterial.init_as(1))
    assert pipe.material._category == 'Pipe'
    assert pipe.material._material == 'Geoperformx'

    # Assert that set_material is working as expected.

    pipe.set_material(PipeMaterial.init_as(0))
    assert pipe.material._category == 'Pipe'
    assert pipe.material._material == 'HDPE'

    # Assert that the material thermal properties are set correctly.

    pipe.material.kth = 0.63
    pipe.material.Cp = 1230
    assert pipe.material.kth == 0.63
    assert pipe.material.Cp == 1230


if __name__ == "__main__":
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
