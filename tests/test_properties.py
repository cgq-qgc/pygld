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


# Qt Test Fixtures
# --------------------------------

@pytest.fixture
def hcf_bot(qtbot):
    hcf = HeatCarrierFluid()
    return hcf, qtbot


# Test RawDataDownloader
# -------------------------------

def test_water(hcf_bot):
    hcf, qtbot = hcf_bot
    assert hcf

    hcfluid = HeatCarrierFluid('water', 28)
    expected_results = [
            [-10, 998.13, 4272, nan, 0.0026477]]

    for er in expected_results:
        hcfluid.Tref = er[0]
        # Assert primary properties.
        assert np.round(hcfluid.rho, 2) == er[1]
        assert np.round(hcfluid.cp, 3) == er[2]
        if np.isnan(er[3]):
            assert np.isnan(hcfluid.k)
        else:
            assert np.round(hcfluid.k, 4) == er[3]
        assert hcfluid.mu == er[4]

        # Assert derived properties.
        assert np.round(hcfluid.Cp, 2) == 4264011.36


if __name__ == "__main__":
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
    # pytest.main()
