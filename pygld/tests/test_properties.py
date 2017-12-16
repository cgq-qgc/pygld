# -*- coding: utf-8 -*-

# copyright (c) 2016 Louis Lamarche
# copyright (c) 2016-2017 Jean-Sébastien
# https://github.com/jnsebgosselin/pygld
#
# This is part of PyGLD (Python Ground-Loop Designer).
# Licensed under the terms of the MIT License.


# ---- Imports: standard libraries

import sys
import os
from itertools import product


# ---- Imports: third parties

import numpy as np
from numpy import nan
import pytest


# ---- Imports: local

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from pygld import HeatCarrierFluid


# Test RawDataDownloader
# -------------------------------

def test_antifreeze():
    hcfluid = HeatCarrierFluid('water', 28)
    hcfluid.fr = 0.3
    with pytest.raises(ValueError):
        hcfluid.fr = -1
    with pytest.raises(ValueError):
        hcfluid.fr = 2
    assert hcfluid.fr == 0

    hcfluid.fluid = 'prop_glycol'
    assert hcfluid.fr == 0.3


def test_water():
    """ Test water properties."""
    hcfluid = HeatCarrierFluid('water', 28)
    expected_results = [[-10, 998.13, 4272, nan, 0.0026477, 4264011.36,
                         nan, 2.6526604750884153e-06, nan],
                        [10, 999.70, 4195, 0.5820, 0.0013059, 4193741.5,
                         9.412801546391753, 1.30629188756627e-06,
                         1.3877822464737036e-07]]

    for expected_result in expected_results:
        hcfluid.Tref = expected_result[0]

        # Assert freezing and boiling point temperature.
        assert hcfluid.Tfp == 0
        assert hcfluid.Tbp == 100

        # Assert primay and derived properties.
        result = [hcfluid.Tref, hcfluid.rho, hcfluid.cp, hcfluid.k,
                  hcfluid.mu, hcfluid.Cp, hcfluid.Pr, hcfluid.nu,
                  hcfluid.al]
        for val, exp_val in zip(result, expected_result):
            if np.isnan(exp_val):
                assert np.isnan(val)
            else:
                assert abs(val - exp_val)/exp_val*100 < 0.01


if __name__ == "__main__":
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
    # pytest.main()
