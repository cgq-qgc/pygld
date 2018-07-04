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

from pygld.heatpumps.utils import build_database
from pygld.utils.fileio import delete_file_safely
from pygld import __rootdir__

DBDIRNAME = osp.join(__rootdir__, 'heatpumps', 'data')
DBFNAME = osp.join(DBDIRNAME, 'hp_database.npy')


def test_build_database(qtbot):
    """
    Test that the database is building as expected.
    """
    delete_file_safely(DBFNAME)
    assert not osp.exists(DBFNAME)

    # Assert that an error is raised if the provided directory name does not
    # exist.

    with pytest.raises(FileNotFoundError):
        build_database('dummy')

    # Build the database and assert that a file was created as expected.

    filename = build_database(DBDIRNAME)
    assert filename == DBFNAME
    assert osp.exists(DBFNAME)


if __name__ == "__main__":
    pytest.main(['-x', os.path.basename(__file__), '-v', '-rw'])
