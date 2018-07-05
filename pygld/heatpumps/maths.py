# -*- coding: utf-8 -*-

# Copyright © 2018 Jean-Sébastien Gosselin
# https://github.com/jnsebgosselin/qwatson
#
# This file is part of PyGLD.
# Licensed under the terms of the GNU General Public License.

# ---- Standard imports

# ---- Third party imports

import numpy as np

# ---- Local imports


def multi_polyfit2rd(y, x1, x2):
    """
    Calculate the coefficient of a second order polynomial expression in two
    variables of the form :

    y = a1 + a2*x1 + a3*x1^2 + a4*x2, a5*x2^2 + a6*x1*x2

    where y is the dependent variable, x1 and x2 are the independent variables,
    and ai are the coefficients of the model.
    """

    # Remove nan values in the dataset.

    indx = np.where(~np.isnan(y))[0]
    y = y[indx]
    x1 = x1[indx]
    x2 = x2[indx]

    # Organize the data in the form : Ax = y

    x = np.column_stack([np.ones(len(y)), x1, x1**2, x2, x2**2, x1*x2])
    A = np.linalg.lstsq(x, y)[0]

    return A


def eval_polyfid2rd(A, x1, x2):
    """
    Use the list of coefficients A to compute values of the dependent variable
    from arrays of values of the independent variables x1 and x2.
    """
    try:
        X = np.hstack([np.ones(len(x1)), x1, x1**2, x2, x2**2, x1*x2])
    except TypeError:
        X = np.array([1, x1, x1**2, x2, x2**2, x1*x2])
    return np.dot(A, X)
