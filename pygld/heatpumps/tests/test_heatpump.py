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

from pygld.heatpumps.heatpump import HeatPump


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


if __name__ == "__main__":
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
