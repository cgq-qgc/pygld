# -*- coding: utf-8 -*-

# Copyright © 2018 Jean-Sébastien Gosselin
# https://github.com/jnsebgosselin/qwatson
#
# This file is part of PyGLD.
# Licensed under the terms of the GNU General Public License.

# ---- Standard imports

# ---- Third party imports

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# ---- Local imports

from pygld.heatpumps.utils import load_heatpump_database
from pygld.heatpumps.heatpump import HeatPump


def plot_fitmodel_eval_from(hpdata):
    """
    Compared the measured COP and CAP values with those evaluated with the
    equation-fit model.
    """
    hpname = hpdata['model']
    title = ("Evaluation of equation-fit models\n"
             "for heatpump %s") % hpname
    fig, axes = plt.subplots(2, 2)
    fig.set_size_inches(6, 6)
    fig.canvas.set_window_title(title)
    fig.suptitle(title)

    varnames = ['CAPc', 'CAPh', 'COPc', 'COPh']
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#ff7f00']
    iterables = zip(varnames, axes.flatten(), colors)
    for i, (varname, axe, color) in enumerate(iterables):
        y = hpdata[varname]

        # Remove the nan values from the dataset.

        indx = np.where(~np.isnan(y))[0]
        y = y[indx]
        x1 = hpdata['EWT'][indx]
        x2 = hpdata['GPM'][indx]

        # Predict the heatpumps COP or CAP

        A = hpdata['models'][varname]
        yp = (A[0] +
              A[1]*x1 + A[2]*x1**2 +
              A[3]*x2 + A[4]*x2**2 +
              A[5]*x1*x2
              )

        # Plot the comparison between the predicted and measured data.

        axe.plot(y, yp, 'o', ms=3, color=color, clip_on=False, zorder=10)
        axe.set_xlabel('measured %s' % varname)
        axe.set_ylabel('predicted %s' % varname)

        axe.set_aspect('equal')
        axe.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))
        axe.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))

        ticks_pos = axe.get_xticks()
        axe.set_yticks(ticks_pos)
        axe.set_xticks(ticks_pos)

        axe.plot(axe.get_xlim(), axe.get_ylim(), '--k', clip_on=False, lw=1,
                 zorder=1)

    fig.subplots_adjust(
        left=0.09, right=1-0.01, bottom=0.08, wspace=0.3, hspace=0.3,
        top=0.89)

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
    plt.close('all')

    database = load_heatpump_database()
    hpnames = list(database.keys())
    for hpname in hpnames:
        plot_fitmodel_eval_from(database[hpname])
