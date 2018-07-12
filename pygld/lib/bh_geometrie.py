# -*- coding: utf-8 -*-

# Copyright (c) PyGLD Project Contributors
# https://github.com/jnsebgosselin/pygld
#
# This is part of PyGLD (Python Ground-Loop Designer).
# Licensed under the terms of the MIT License.

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


TLB = {'NPS': np.array([0.5, 0.75, 1, 1.25, 1.5, 2,
                        2.5, 3, 3.5, 4, 5, 6, 8]),
       'OD': np.array([0.840, 1.050, 1.315, 1.660, 1.900, 2.375, 2.875,
                       3.500, 4.000, 4.500, 5.563, 6.625, 8.625, 10.750,
                       12.750, 14.000, 16.000, 18.000]),
       'ODtol': np.array([0.004, 0.004, 0.005, 0.005, 0.006, 0.006, 0.007,
                          0.008, 0.009, 0.009, 0.011, 0.011, 0.013, 0.015,
                          0.017, 0.063, 0.072, 0.081])}


def calcul_sdr(NPS, SDR):
    """
    Take a NPS (in inches) and SDR values and return the outside OD and
    inside ID diameters of the pipe in meters.

    The outside diameters (OD) are specified as in Table 2 of:

        ANSI B36.10-1979: American National Standard for Welded and
            Seamless Wrought Steel Pipe.

    The outside diameter tolerances (ODtol) are specified as in Table 2 of:

        ASTM D3035-15: Standard Specification for Polyethylene (PE)
            Plastic Pipe (DR-PR) Based on Controlled Outside Diameter.

    The outside diameters and tolerances for the IPS values not
    listed in Table 2 of ASTM D3035-1 (2 1/2, 3 1/2, and 5)
    are specified such that their tolerances correspond to the same
    percentage of the outside diameter as those for the closest
    listed diameter.

    ID, OD = calcul_sdr(NPS, SDR)

    Input:
    ------
    NPS = Nominal Pipe Size in inches
    SDR = Standard Dimension Ratio

    Output:
    -------
    ID = Inside diameter in m
    OD = Outside diameter in m
    """

    # Looking up in the table to find the corresponding sizes :

    if NPS in TLB['NPS']:
        indx = np.where(NPS == TLB['NPS'])[0][0]

        OD = TLB['OD'][indx]

        tmin = max(0.060, OD/SDR)
        ttol = max(0.02, tmin*0.12)
        tmean = tmin + ttol/2

        ID = max(0, OD - (2*tmean))

        return ID*0.0254, OD*0.0254

    else:
        nps = np.round(TLB['NPS'], 2).astype(str).tolist()
        nps = 'in, '.join(nps).strip('[]')
        raise Exception('Wrong NPS value. Valid values for NPS are: '
                        '%s in.' % nps)


def plot_borehole(rb, zcp, rpin, rpext):

    fig = plt.figure(FigureClass=FigBoreholeGeo)
    fig.update_geo(rb, zcp, rpin, rpext)

    return fig


class FigBoreholeGeo(mpl.figure.Figure):

    def __init__(self, *args, **kwargs):
        super(FigBoreholeGeo, self).__init__(*args, **kwargs)

        self.set_facecolor('white')
        self.add_axes([0.05, 0.05, 0.9, 0.9], frameon=True, aspect='equal')
        self.set_size_inches(5, 5)

    def update_geo(self, rb, zcp, rpin, rpext):  # ============================

        ax = self.axes[0]
        patches = ax.patches

        xcp = np.real(zcp)
        ycp = np.imag(zcp)

        if len(ax.patches) == (2*len(zcp)+1):

            # The number of artists haven't changed. Simply update the
            # coordinate and radius.

            # Update borehole:
            patches[0].set_radius(rb)

            # Update outer pipe radius:
            for i, patch in enumerate(patches[1::2]):
                patch.set_radius(rpext[i])
                patch.center = xcp[i], ycp[i]

            # Update inner pipe radius:
            for i, patch in enumerate(patches[2::2]):
                patch.set_radius(rpin[i])
                patch.center = xcp[i], ycp[i]

        else:
            ax.clear()
            ax.axis('off')

            # Create borehole:
            ax.add_patch(mpl.patches.Circle(
                (0, 0), rb, fc='0.65', ec='0.25', ls='--', clip_on=False))

            # Create pipes:
            for i in range(len(zcp)):
                ax.add_patch(mpl.patches.Circle(
                    (xcp[i], ycp[i]), rpext[i], fc='0.25', ec='None',
                    clip_on=False))

                ax.add_patch(mpl.patches.Circle(
                    (xcp[i], ycp[i]), rpin[i], color='#9ECAE1', clip_on=False))

        ax.axis([-rb, rb, -rb, rb])

        self.canvas.draw()


# ---- test

def test_module():
    print('\n' + '-' * 79)
    print('Validation "calcul_sdr" with Table2 in Remund, 1999')
    print('-' * 79)
    for NPS in [0.75, 1, 1.25, 1.5, 2.]:
        for SDR in [9, 11, 15.5]:
            Din, Dout = calcul_sdr(NPS, SDR)
            print(('SDR = %02d; NPS = %0.2f in ; '
                   'ID = %0.4f m ; OD = %0.4f m'
                   ) % (SDR, NPS, Din, Dout))

    # -------------------------------------------------------------------------

    plt.close('all')

    # ---- Double U-Pipe 1.25 in ----

    rb = 152.4/2
    zcp = [0+50j, 0-50j, 50+0j, -50-0j]
    rpin = [33.99/2] * 4
    rpext = [42.16/2] * 4

    plot_borehole(rb, zcp, rpin, rpext)

    # ---- Coaxial ----

    rb = 152.4
    zcp = [0, 0]
    rpin = [100.03, 48.69]
    rpext = [114.29, 60.32]

    plot_borehole(rb, zcp, rpin, rpext)

    plt.show()

if __name__ == '__main__':

    test_module()
