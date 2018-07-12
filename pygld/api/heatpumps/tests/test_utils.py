# -*- coding: utf-8 -*-

# Copyright © 2018 Jean-Sébastien Gosselin
# https://github.com/jnsebgosselin/qwatson
#
# This file is part of PyGLD.
# Licensed under the terms of the GNU General Public License.

# ---- Standard imports

import os
import os.path as osp

# ---- Third party imports

import pytest

# ---- Local imports

from pygld.api.heatpumps.utils import build_heatpump_database
from pygld.utils.fileio import delete_file_safely
from pygld.api.heatpumps import __datadir__


def test_build_database():
    """
    Test that the database is building as expected.
    """
    expected_filename = osp.join(__datadir__, 'hp_database.npy')
    delete_file_safely(expected_filename)
    assert not osp.exists(expected_filename)

    # Assert that an error is raised if the provided directory name does not
    # exist.

    with pytest.raises(FileNotFoundError):
        build_heatpump_database('dummy')

    # Build the database and assert that a file was created as expected.

    filename = build_heatpump_database(__datadir__)
    assert filename == expected_filename
    assert osp.exists(expected_filename)


if __name__ == "__main__":
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
