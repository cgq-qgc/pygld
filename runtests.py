# -*- coding: utf-8 -*-

# copyright (c) 2016 Louis Lamarche
# copyright (c) 2016-2017 Jean-Sébastien
# https://github.com/jnsebgosselin/pygld
#
# This is part of PyGLD (Python Ground-Loop Designer).
# Licensed under the terms of the MIT License.


"""
File for running tests programmatically.
"""

import pytest


def main():
    """
    Run pytest tests.
    """
    errno = pytest.main(['-x', 'pygld',  '-v', '-rw', '--durations=10',
                         '--cov=pygld'])

    if errno != 0:
        raise SystemExit(errno)


if __name__ == '__main__':
    main()
