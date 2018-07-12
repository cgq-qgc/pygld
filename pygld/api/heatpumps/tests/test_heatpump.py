# -*- coding: utf-8 -*-

# Copyright © 2018 Jean-Sébastien Gosselin
# https://github.com/jnsebgosselin/qwatson
#
# This file is part of PyGLD.
# Licensed under the terms of the GNU General Public License.

# ---- Standard imports

import os

# ---- Third party imports

import pytest
import numpy as np

# ---- Local imports

from pygld import HeatPump
from pygld.api.materials.heatcarrierfluids import FLUIDS


def test_heatpump_init():
    """
    Test that the heatpump class is initializing as expected with the
    proper default values.
    """
    heatpump = HeatPump()
    assert heatpump

    # Assert the default values for the independent variables.

    assert np.array_equal(heatpump.TinHP, np.array([28, 0]))
    assert np.array_equal(heatpump.qbat, np.array([-16.5, 14.5]))
    assert np.array_equal(heatpump.Vf, np.array([0.94635295, 0.94635295]))
    assert heatpump.fluid == 'water'
    assert heatpump.fr == 0

    # Assert the default values for the dependent variables.

    assert np.array_equal(np.round(heatpump.COP, 6),
                          np.array([3.849727, 3.372328]))
    assert np.array_equal(np.round(heatpump.CAP, 6),
                          np.array([20.926771, 17.196401]))
    assert np.array_equal(np.round(heatpump.ToutHP, 6),
                          np.array([33.273623, -2.555170]))
    assert np.array_equal(np.round(heatpump.Tm, 6),
                          np.array([30.636812, -1.277585]))


def test_heatpump_need_update_flags():
    """
    Test that the flag that indicate that an update of the dependent variables
    is required is resetted correctly when the value of an independent
    variable change.
    """
    heatpump = HeatPump()

    for key, value in heatpump._need_update.items():
        assert value is True

    # We print heatpump so that all dependent properties ares accessed.

    print(heatpump)

    for key, value in heatpump._need_update.items():
        assert value is False

    heatpump.TinHP = 34

    for key, value in heatpump._need_update.items():
        assert value is True


def test_lenght_notequal_error():
    """
    Test that an error is raised when TinHP, qtbat, and Vf length does not
    match perfectly.
    """
    heatpump = HeatPump()
    heatpump.TinHP = 28

    with pytest.raises(ValueError):
        print(heatpump.ToutHP)

    heatpump.TinHP = [28, 0]
    heatpump.qbat = [-16.5, 14.5, 12]
    with pytest.raises(ValueError):
        print(heatpump.ToutHP)


def test_single_value():
    """
    Test that everything is working as expected when working with
    single values
    """
    heatpump = HeatPump()
    heatpump.TinHP = 28
    heatpump.Vf = 0.94635295
    heatpump.qbat = -16.5

    assert np.round(heatpump.ToutHP, 6) == 33.273623
    assert np.round(heatpump.COP, 6) == 3.849727
    assert np.round(heatpump.CAP, 6) == 20.926771
    assert np.round(heatpump.Tm, 6) == 30.636812


def test_public_interface_insulation():
    """
    Test that it is not possible to modify the private variable outside of
    the public interface.
    """
    heatpump = HeatPump()

    # Independent variables.

    TinHP = np.array([1, 2, 3, 4])
    heatpump.TinHP = TinHP
    assert np.array_equal(heatpump.TinHP, TinHP)
    assert id(TinHP) != id(heatpump.TinHP)
    assert id(TinHP) != id(heatpump._TinHP)
    assert id(heatpump.TinHP) != id(heatpump._TinHP)

    Vf = np.array([5, 4, 3, 2])
    heatpump.Vf = Vf
    assert np.array_equal(heatpump.Vf, Vf)
    assert id(Vf) != id(heatpump.Vf)
    assert id(Vf) != id(heatpump._Vf)
    assert id(heatpump.Vf) != id(heatpump._Vf)

    qbat = np.array([10, -14, 30, -32])
    heatpump.qbat = qbat
    assert np.array_equal(heatpump.qbat, qbat)
    assert id(qbat) != id(heatpump.qbat)
    assert id(qbat) != id(heatpump._qbat)
    assert id(heatpump.qbat) != id(heatpump._qbat)

    fluid = 'ethyl_glycol'
    heatpump.fluid = fluid
    assert heatpump.fluid == 'ethyl_glycol'
    fluid = 'water'
    assert heatpump.fluid == 'ethyl_glycol'

    fr = 0.5
    heatpump.fr = fr
    assert heatpump.fr == 0.5
    fluid = 0.25
    assert heatpump.fr == 0.5

    # Dependent variables.

    assert id(heatpump.ToutHP) != id(heatpump._ToutHP)
    assert id(heatpump.Tm) != id(heatpump._Tm)
    assert id(heatpump.COP) != id(heatpump._COP)
    assert id(heatpump.CAP) != id(heatpump._CAP)


def test_get_set_model():
    """
    Test that the method to set the model of the heatpump is working as
    expected.
    """
    heatpump = HeatPump()
    model_list = heatpump.get_avail_heatpump_models()

    # Change the model with an index.

    print(heatpump)
    for key, value in heatpump._need_update.items():
        assert value is False

    heatpump.set_model(3)
    assert heatpump.model == model_list[3]
    for key, value in heatpump._need_update.items():
        assert value is True

    # Change the model with a str.

    print(heatpump)
    for key, value in heatpump._need_update.items():
        assert value is False

    heatpump.set_model(model_list[1])
    assert heatpump.model == model_list[1]
    for key, value in heatpump._need_update.items():
        assert value is True

    # Test errors.

    expected_model = heatpump.model

    with pytest.raises(IndexError):
        heatpump.set_model(999999)
    assert heatpump.model == expected_model

    with pytest.raises(ValueError):
        heatpump.set_model('dummy_model_name')
    assert heatpump.model == expected_model

    with pytest.raises(TypeError):
        heatpump.set_model(12.34)
    assert heatpump.model == expected_model


def test_avail_fluid_types():
    """
    Test that the method to get a list of available heat carrier fluid types
    is working as expected.
    """
    heatpump = HeatPump()
    assert heatpump.get_avail_fluid_types() == FLUIDS
    assert id(heatpump.get_avail_fluid_types()) != id(FLUIDS)


if __name__ == "__main__":
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
