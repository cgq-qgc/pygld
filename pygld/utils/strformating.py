# -*- coding: utf-8 -*-

# Copyright © 2018 Jean-Sébastien Gosselin
# https://github.com/jnsebgosselin/pygld
#
# This file is part of PyGLD.
# Licensed under the terms of the MIT License.


def array_to_str(array, str_fmt='{0:.1f}', sep=', '):
    """
    Take a numpy array and format each element into a nicely formatted
    string.
    """
    prefix = '[' if len(array) > 1 else ''
    suffix = ']' if len(array) > 1 else ''
    return (prefix +
            sep.join([str_fmt.format(v, i) for i, v in enumerate(array)]) +
            suffix)
