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

    propnames = ['TinHP', 'qbat', 'Vf']
    expected_vals = [(28, 0), (16.5, 14.5), (0.94635295, 0.94635295)]
    for (propname, expected_val) in zip(propnames, expected_vals):
        prop = getattr(heatpump, propname)
        assert id(prop) == id(getattr(heatpump, propname))

        assert (prop.c == prop.cooling == prop['c'] == prop['cooling'] ==
                prop[0] == expected_val[0])
        assert (prop.h == prop.heating == prop['h'] == prop['heating'] ==
                prop[1] == expected_val[1])

    assert heatpump.fluid == 'water'
    assert heatpump.fr == 0

    # Assert the default values for the dependent variables.

    propnames = ['COP', 'CAP', 'ToutHP', 'Tm']
    expected_vals = [(3.849727, 3.372328), (20.926771, 17.196401),
                     (33.273623, -2.555170), (30.636812, -1.277585)]
    for (propname, expected_val) in zip(propnames, expected_vals):
        prop = getattr(heatpump, propname)
        assert id(prop) == id(getattr(heatpump, propname))

        assert (round(prop.c, 6) == round(prop.cooling, 6) ==
                round(prop['c'], 6) == round(prop['cooling'], 6) ==
                round(prop[0], 6) == expected_val[0])

        assert (round(prop.h, 6) == round(prop.heating, 6) ==
                round(prop['h'], 6) == round(prop['heating'], 6) ==
                round(prop[1], 6) == expected_val[1])


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

    heatpump.TinHP.c = 34

    for key, value in heatpump._need_update.items():
        assert value is True


def test_calcul_Tm():
    """Test that the average fluid temperature is calculated correctly."""
    heatpump = HeatPump()
    assert 0.5 * (heatpump.TinHP.c + heatpump.ToutHP.c) == heatpump.Tm.c
    assert 0.5 * (heatpump.TinHP.h + heatpump.ToutHP.h) == heatpump.Tm.h


if __name__ == "__main__":
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
