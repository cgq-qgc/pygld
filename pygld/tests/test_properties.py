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
    expected_results = [[-10, 0, 100, 998.13, 4272, nan, 0.0026477,
                         4264011.36, nan, 2.6526604750884153e-06, nan],
                        [10, 0, 100, 999.70, 4195, 0.5820, 0.0013059,
                         4193741.5, 9.412801546391753, 1.30629188756627e-06,
                         1.3877822464737036e-07]]

    for expected_result in expected_results:
        hcfluid.Tref = expected_result[0]
        result = [hcfluid.Tref, hcfluid.Tfp, hcfluid.Tbp, hcfluid.rho,
                  hcfluid.cp, hcfluid.k, hcfluid.mu, hcfluid.Cp, hcfluid.Pr,
                  hcfluid.nu, hcfluid.al]
        for val, exp_val in zip(result, expected_result):
            if np.isnan(exp_val):
                assert np.isnan(val)
            elif exp_val == 0:
                assert abs(val - exp_val) < 0.001
            else:
                assert abs(val - exp_val)/exp_val*100 < 0.01


def test_freezing_point_propglycol():
    """Test propylene glycol freezing and boiling point."""
    hcfluid = HeatCarrierFluid('prop_glycol')
    arr_fr = [0.304, 0.204, 0.096, 0]
    arr_Tfp = [-13.4, -7.6, -3.3, 0]
    arr_Tbp = [102.2, 100.6, 100, 100]
    for fr, Tfp, Tbp in zip(arr_fr, arr_Tfp, arr_Tbp):
        hcfluid.fr = fr
        if Tfp == 0:
            assert hcfluid.Tfp == 0
        else:
            assert abs(hcfluid.Tfp - Tfp)/Tfp*100 < 0.01
        assert abs(hcfluid.Tbp - Tbp)/Tbp*100 < 0.01


def test_freezing_point_ethylglycol():
    """Test ethylene glycol freezing and boiling point."""
    hcfluid = HeatCarrierFluid('ethyl_glycol')
    arr_fr = [0.306, 0.201, 0.089, 0]
    arr_Tfp = [-16.2, -8.9, -3.2, 0]
    arr_Tbp = [104.4, 102.2, 101.1, 100]
    for fr, Tfp, Tbp in zip(arr_fr, arr_Tfp, arr_Tbp):
        hcfluid.fr = fr
        if Tfp == 0:
            assert hcfluid.Tfp == 0
        else:
            assert abs(hcfluid.Tfp - Tfp)/Tfp*100 < 0.01
        assert abs(hcfluid.Tbp - Tbp)/Tbp*100 < 0.01



if __name__ == "__main__":
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
    # pytest.main()
