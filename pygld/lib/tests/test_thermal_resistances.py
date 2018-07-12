# -*- coding: utf-8 -*-

# Copyright (c) PyGLD Project Contributors
# https://github.com/jnsebgosselin/pygld
#
# This is part of PyGLD (Python Ground-Loop Designer).
# Licensed under the terms of the MIT License.


# ---- Standard imports

import sys
import os

# ---- Third party imports

import numpy as np
from numpy import pi
import pytest


# ---- Local imports

from pygld.lib.thermal_resistances import calcul_Rb_multipoles


# ---- Test calcul_Rb_multipoles

def test_multipole_with_claesson_2011():
    """
    Validation calcul of Rb with the multipoles method by comparing the
    results with Table 1 in Claesson and Hellstrom (2011).
    """
    kb = 1.5   # thermal conductivity of the grout
    k = 2.5    # thermal conductivity of the pipes
    Rp = 1.2/(2*pi*kb)  # thermal resistance of the pipes
    rp = 0.02  # outside radius of the pipes
    zcp = [0.03+0j, -0.03+0.02j]  # position of the pipe center

    rb = 0.07  # radius of the borehole

    # Published values in Table 1 of Claesson and Hellstrom, 2011 :
    Kb = np.array([7.439, 7.410, 7.405, 7.404])
    K1b = np.array([3.698, 3.683, 3.681, 3.680])
    K2b = np.array([3.742, 3.727, 3.724, 3.724])
    K12 = np.array([0.240, 0.243, 0.242, 0.242])

    expected_Rb = 1/Kb
    expected_Ra = (1/K12 * (1/K1b + 1/K2b))/(1/K1b + 1/K2b + 1/K12)

    for Jp in range(4):
        Rb, Ra = calcul_Rb_multipoles(kb, k, rb, rp, Rp, Jp, zcp)
        assert np.round(Rb, 4) == np.round(expected_Rb[Jp], 4)
        assert np.round(Ra, 2) == np.round(expected_Ra[Jp], 2)


if __name__ == "__main__":
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
