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
from numpy import nan
import pytest


# ---- Local imports

from pygld.fluidproperties import HeatCarrierFluid


def test_antifreeze_and_fluid_error():
    """Test that error are raised correctly for Tref and fr."""
    hcfluid = HeatCarrierFluid('water', 28)
    hcfluid.fr = 0.3
    with pytest.raises(ValueError):
        hcfluid.fr = -1
    with pytest.raises(ValueError):
        hcfluid.fr = 2
    assert hcfluid.fr == 0

    hcfluid.set_fluid('prop_glycol')
    assert hcfluid.fr == 0.3

    with pytest.raises(ValueError):
        hcfluid.set_fluid('test')


def test_properties_of_water():
    """ Test properties of Pure Water."""
    hcfluid = HeatCarrierFluid('water', 28)
    expected_results = [[-10, 0, 100, 998.13, 4272, nan, 0.0026477,
                         4264011.36, nan, 2.6526604750884153e-06, nan],
                        [10, 0, 100, 999.70, 4195, 0.5820, 0.0013059,
                         4193741.5, 9.412801546391753, 1.30629188756627e-06,
                         1.3877822464737036e-07]]

    for expected_result in expected_results:
        hcfluid.Tref = expected_result[0]
        result = [hcfluid.Tref, hcfluid.Tfp, hcfluid.Tbp, hcfluid.rho,
                  hcfluid.cp, hcfluid.kth, hcfluid.mu, hcfluid.Cp, hcfluid.Pr,
                  hcfluid.nu, hcfluid.al]
        for val, exp_val in zip(result, expected_result):
            if np.isnan(exp_val):
                assert np.isnan(val)
            elif exp_val == 0:
                assert abs(val - exp_val) < 0.001
            else:
                assert abs(val - exp_val)/exp_val*100 < 0.01


def test_freezing_point_of_propglycol():
    """Test the freezing and boiling point of Propylene Glycol."""
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


def test_freezing_point_of_ethylglycol():
    """Test the freezing and boiling point of Ethylene Glycol."""
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


def test_properties_of_propglycol():
    """Test the properties of Propylene Glycol."""
    hcfluid = HeatCarrierFluid('prop_glycol')
    hcfluid.Tref = [-5, 10]
    hcfluid.fr = 0.3

    assert hcfluid.rho.tolist() == [1037.89, 1032.55]
    assert hcfluid.cp.tolist() == [3779, 3820]
    assert hcfluid.kth.tolist() == [0.403, 0.421]
    assert hcfluid.mu.tolist() == [0.00907, 0.00452]


def test_properties_of_ethylglycol():
    """Test the properties of Ethylene Glycol."""
    hcfluid = HeatCarrierFluid('ethyl_glycol')
    hcfluid.Tref = -5
    hcfluid.fr = 0.3

    assert abs(hcfluid.rho - 1053.11) < 10**-6
    assert abs(hcfluid.cp - 3574) < 10**-6
    assert abs(hcfluid.kth - 0.417) < 10**-6
    assert abs(hcfluid.mu - 0.00503) < 10**-6


if __name__ == "__main__":
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
    # pytest.main()
