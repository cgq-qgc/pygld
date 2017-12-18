# -*- coding: utf-8 -*-

# Copyright (c) PyGLD Project Contributors
# https://github.com/jnsebgosselin/pygld
#
# This is part of PyGLD (Python Ground-Loop Designer).
# Licensed under the terms of the MIT License.


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

try:
    from geothermie.g_function import G_function_ILS, G_function_FLS
except:
    import sys
    from os.path import dirname, realpath
    root = dirname(dirname(realpath(__file__)))
    sys.path.append(root)
    from geothermie.g_function import G_function_ILS, G_function_FLS


# =============================================================================


def calcul_Tp_ils(x, y, qp, al, ks, tf):
    # Calcul the penality temperature with the superposition principle using
    # the infinite line source model.

    # x, y = horizontal position of the wells in m.
    # al = soil thermal diffusivity (m²/day)
    # ks = soil thermal conductivity (W/mK)
    # tf = time in days

    N = len(x)

    Tpp = np.zeros(N)
    rb_mem = []
    GILS_mem = []
    for i in range(N):
        for j in range(N):
            if i != j:
                rb = ((x[i]-x[j])**2 + (y[i]-y[j])**2)**0.5
                Fo = al * tf / rb**2
                indx = np.where(rb == rb_mem)[0]
                if len(indx) == 0:
                    GILS = G_function_ILS(Fo)
                    GILS_mem.append(GILS)
                    rb_mem.append(rb)
                else:
                    GILS = GILS_mem[indx[0]]

                Tpp[i] = Tpp[i] + qp/ks*GILS

    Tp = np.mean(Tpp)

    return Tp


# =============================================================================


def calcul_Tp_fls(x, y, qp, al, ks, tf, H):
    # Calcul the penality temperature with the superposition principle using
    # the finite line source model.

    # x, y = horizontal position of the wells in m.
    # al = soil thermal diffusivity (m²/day)
    # ks = soil thermal conductivity (W/mK)
    # tf = final time in days
    # H = depth of the wells in m

    N = len(x)

    Tpp = np.zeros(N)
    rb_mem = []
    GFLS_mem = []
    for i in range(N):
        for j in range(N):
            if i != j:
                rb = ((x[i]-x[j])**2 + (y[i]-y[j])**2)**0.5
                Fo = al * tf / rb**2
                indx = np.where(rb == rb_mem)[0]
                if len(indx) == 0:
                    GFLS = G_function_FLS(Fo, rb/H)
                    GFLS_mem.append(GFLS)
                    rb_mem.append(rb)
                else:
                    GFLS = GFLS_mem[indx[0]]

                Tpp[i] = Tpp[i] + qp/ks * GFLS

    Tp = np.mean(Tpp)

    return Tp


# =============================================================================


def design_borefield(Nbh, dbb, mode='sequence'):
    # Nbh = Total number of boreholes
    # dbb = Horizontal distance between boreholes

    # ---- Core boreholes ---- #

    nc = int(np.floor(Nbh**0.5))
    x, y = [], []
    for row in range(nc):
        for col in range(nc):
            y.append(row*dbb)
            x.append(col*dbb)

    # ---- Edge boreholes ---- #

    if mode == 'sequence':
        col = 0
        row = nc
        while True:
            if len(x) == Nbh:
                break
            y.append(row*dbb)
            x.append(col*dbb)
            if col >= nc:
                row += -1
            else:
                col += 1

    if mode == 'alternate':
        count = 0
        while True:
            if len(x) == Nbh:
                break
            y.append(nc*dbb)
            x.append(count*dbb)
            if len(x) == Nbh:
                break
            y.append(count*dbb)
            x.append(nc*dbb)

            count += 1

    # ---- Return Results ---- #

    x = np.array(x).astype(float)
    y = np.array(y).astype(float)

    if len(x) != Nbh:
        raise ValueError

    return x, y


###############################################################################


def exemple5_1():

    # Exemple 5.1 of the course of Louis Lamarche at ETS, Montreal, Qc, Can.

    al = 0.07    # soil thermal diffusivity in m2/day
    ks = 2.25    # soil thermal conductivity in W/m K
    tf = 10*365  # total time in days

    qp = 7.2  # W/m
    nx = 6    # column count
    ny = 5    # row count
    d = 6     # horizontal distance between well in m
    H = 100   # depth of the wells in m

    # Building a grid:

    N = nx*ny
    x = np.zeros(N)
    y = np.zeros(N)
    k = 0
    for i in range(nx):
        for j in range(ny):
            x[k] = d*i
            y[k] = d*j
            k += 1

    x, y = design_borefield(30, 6)

    Tp3 = calcul_Tp_ils(x, y, qp, al, ks, tf)
    print('Tp3 =', Tp3, ' vs ', 8.61082077493)

    Tp4 = calcul_Tp_fls(x, y, qp, al, ks, tf, H)
    print('Tp4 =', Tp4, ' vs ', 7.18406563754)

if __name__ == '__main__':
    exemple5_1()
