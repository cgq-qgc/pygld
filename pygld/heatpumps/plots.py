# -*- coding: utf-8 -*-

# Copyright © 2018 Jean-Sébastien Gosselin
# https://github.com/jnsebgosselin/qwatson
#
# This file is part of PyGLD.
# Licensed under the terms of the GNU General Public License.

# ---- Standard imports

import os
import csv
from collections import OrderedDict

# ---- Third party imports

import numpy as np
from scipy import interpolate
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay

# ---- Local imports

from pygld.heatpumps.heatpump import HeatPump


def eval_model():

    plt.close('all')

    w = HeatPump()

    for name in ['TCHV072', 'TCHV096', 'TCHV120', 'TCHV160', 'TCHV192',
                 'TCHV240', 'TCHV300', 'Generic']:
        w.hpname = name

        y = w.hpdata['CAPh']/w.hpdata['Wh']
        x1 = w.hpdata['EWT']
        x2 = w.hpdata['GPM']
        A = linalg_hp(y, x1, x2)

        yp = (A[0] +
              A[1]*x1 + A[2]*x1**2 +
              A[3]*x2 + A[4]*x2**2 +
              A[5]*x1*x2
              )

        fig, ax = plt.subplots()
        fig.canvas.set_window_title(name)
        ax.plot([np.floor(np.min(yp)), np.ceil(np.max(yp))],
                [np.floor(np.min(yp)), np.ceil(np.max(yp))], '--k')
        ax.plot(y, yp, 'o')

    plt.show()


def plot_cop():
    plt.close('all')

    w = HeatPump()

    w.hpname = 'Generic'
    w.hpname = 'TCHV300'
    # y = w.hpdata['CAPh']/w.hpdata['Wh']
    y = w.hpdata['CAPc']/w.hpdata['Wc']
    x1 = w.hpdata['EWT']
    x2 = w.hpdata['GPM']
    A = linalg_hp(y, x1, x2)

    print(w.get_flowRange())

    fig, ax = plt.subplots()

    ewt = np.arange(-1, 35)
#    flowrange = np.arange(w.get_flowRange()[0], w.get_flowRange()[1], 0.25)
    flowrange = np.arange(2.37, 4.73, 0.25)

    for val in flowrange:
        gpm = np.ones(len(ewt)) * val

        cop = (A[0] +
               A[1]*ewt + A[2]*ewt**2 +
               A[3]*gpm + A[4]*gpm**2 +
               A[5]*ewt*gpm
               )
        ax.plot(ewt, cop, '-')

    plt.show()


if __name__ == '__main__':
    import matplotlib.pyplot as plt

#    eval_model()
#    plot_cop()
#
    w = HeatPump()
#    w.setFixedSize(w.size())
    w.setCurrentUnitSystem('SI')

    w.set_Vftot({'cooling': 12.5/15.8503230745/1000,
                 'heating': 12.5/15.8503230745/1000})
