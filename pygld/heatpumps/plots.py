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
from pygld.heatpumps.maths import eval_polyfid2rd


def plot_fitmodel_eval_from(hpdata):
    """
    Compare the measured COP and CAP values with those evaluated with the
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
        indx = np.where(~np.isnan(y))[0]

        # Predict the heatpumps COP or CAP

        A = hpdata['eqfit_models'][varname]
        yp = eval_polyfid2rd(A, hpdata['EWT'][indx], hpdata['GPM'][indx])

        # Plot the comparison between the predicted and measured data.

        axe.plot(y[indx], yp, 'o', ms=3, color=color, clip_on=False, zorder=10)
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


if __name__ == '__main__':
    from pygld.heatpumps.heatpump import HeatPump
    plt.close('all')
    database = load_heatpump_database()
    hpnames = list(database.keys())
    for hpname in hpnames:
        plot_fitmodel_eval_from(database[hpname])
