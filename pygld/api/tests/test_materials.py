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

# ---- Local imports

from pygld import GroutMaterial, PipeMaterial, GroundMaterial
from pygld.api.materials import PREDEFINED_MATERIALS


def test_init_as_for_groutmaterial():
    """Test that the init_as is working as expected for all materials."""
    matnames = ['Grout', 'Ground', 'Pipe']
    classes = [GroutMaterial, GroundMaterial, PipeMaterial]

    for (matname, class_) in zip(matnames, classes):
        predmats = PREDEFINED_MATERIALS[matname]
        keys = list(predmats.keys())

        mat = class_.init_as(1)
        assert (mat.kth, mat.Cp) == predmats[keys[1]]

        # Assert no error is raised when the material index is out of range.

        mat = class_.init_as(-1)
        assert (mat.kth, mat.Cp) == predmats[keys[0]]
        mat = class_.init_as(len(keys)+1)
        assert (mat.kth, mat.Cp) == predmats[keys[-1]]

        assert mat._category == matname
        assert mat._material == keys[-1]


def test_get_predefined_materials():
    """
    Test that the get_predefined_materials is working as expected for
    all materials.
    """
    matnames = ['Grout', 'Ground', 'Pipe']
    classes = [GroutMaterial, GroundMaterial, PipeMaterial]
    for (matname, class_) in zip(matnames, classes):
        assert (class_.get_predefined_materials() ==
                PREDEFINED_MATERIALS[matname])


def test_getter_setter():
    """
    Test that the interface for setting and gettings the values of the
    properties is working as expected for all materials.
    """
    matnames = ['Grout', 'Ground', 'Pipe']
    classes = [GroutMaterial, GroundMaterial, PipeMaterial]
    for (matname, class_) in zip(matnames, classes):
        mat = class_()
        assert mat._category == matname
        assert mat._material is None
        assert mat.kth is None
        assert mat.Cp is None
        assert mat.al is None

        mat.kth = -5
        mat.Cp = -3567
        assert mat.kth == 5 and type(mat.kth) is float
        assert mat.Cp == 3567 and type(mat.Cp) is float
        assert mat.al == 5/3567


def test_print_raise_no_error():
    """
    Assert that the print_predefined_materials method raise no error. We do not
    actually test the content of the print because that would be too much work.
    """
    matnames = ['Grout', 'Ground', 'Pipe']
    classes = [GroutMaterial, GroundMaterial, PipeMaterial]
    for (matname, class_) in zip(matnames, classes):
        class_.print_predefined_materials()


def test_print_material():
    """
    Test that printing an instance of a class material is working as
    expected.
    """
    matnames = ['Grout', 'Ground', 'Pipe']
    classes = [GroutMaterial, GroundMaterial, PipeMaterial]
    for (matname, class_) in zip(matnames, classes):
        expected_content = (
            "%s material: User Defined\n"
            "Thermal conductivity: 0.670 W/m·k\n"
            "Volumetric heat capacity: 1500 J/m³·K\n"
            "Thermal diffusivity: 4.47e-04 m²/s") % matname
        mat = class_()
        mat.kth = 0.670
        mat.Cp = 1500

    assert mat.__str__() == expected_content


if __name__ == "__main__":
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
