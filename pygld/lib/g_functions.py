# -*- coding: utf-8 -*-

# Copyright (c) PyGLD Project Contributors
# https://github.com/jnsebgosselin/pygld
#
# This is part of PyGLD (Python Ground-Loop Designer).
# Licensed under the terms of the MIT License.

from __future__ import division, unicode_literals

from numpy import pi, log, log10, sqrt, exp, Inf, array
from scipy import integrate
from scipy import special
from scipy.special import j0, y0, j1, y1, erfc


# ---- Short-time analysis of vertical boreholes

def G_function_lamarche(tt, rt, kt, gam, Rt, v):
    # Lamarche, L. 2015. Short-time analysis of vertical boreholes, new
    # analytic solutions and choice of equivalent radius. International
    # Journal of Heat and Mass Transfer, 91, 800-807.
    # Equation #53

    I = integrate.quad(g_int_lamarche, 0, Inf,
                       args=(tt, Rt, rt, kt, gam, v))[0]

    G = 8*kt/pi**5/rt**2 * I

    return G


def g_int_lamarche(x, tt, Rt, rt, kt, gam, v):
    # tt = Fourier number
    # re = outer radius of the equivalent tube for concentric cylinders
    # rb = borehole radius
    # ks = soil thermal conductivity
    # kg = grout thermal conductivity
    # als = soil thermal diffusivity
    # alg = grout thermal diffusivity

    zeta1 = (1 - v*Rt*x**2)*y1(x) - v*x*y0(x)
    zeta2 = (1 - v*Rt*x**2)*j1(x) - v*x*j0(x)

    psi = (zeta2*(j0(x*rt*gam)*y1(x*rt) - j1(x*rt*gam)*y0(x*rt)*kt*gam) -
           zeta1*(j0(x*rt*gam)*j1(x*rt) - j1(x*rt*gam)*j0(x*rt)*kt*gam)
           )

    phi = (zeta1*(y0(x*rt*gam)*j1(x*rt) - y1(x*rt*gam)*j0(x*rt)*kt*gam) -
           zeta2*(y0(x*rt*gam)*y1(x*rt) - y1(x*rt*gam)*y0(x*rt)*kt*gam)
           )

    I = (1 - exp(-x**2*tt)) / (x**5 * (phi**2 + psi**2))

    return I


# ---- Infinite cylindrical heat source solution


def G_function_ICS(Fo, rt):
    """
    Infinite cylindrical heat source solution

    rt = r/rb
    Fo = Fourier number

    Reference:
    ----------
    Lamarche and Beauchamp (2007). A new contribution to the finite
        line-source model for geothermal boreholes. Energy and
        Building, 39:188-198.

    """

    I = integrate.quad(g_int, 0, Inf, args=(Fo, rt))[0]
    G = I/pi**2

    return G


def g_int(x, Fo, rt):
    num = exp(-Fo*x**2) - 1
    den = x**2 * (j1(x)**2 + y1(x)**2)
    I = num/den * (j0(rt*x)*y1(x)-y0(rt*x)*j1(x))

    return I


def g_function_cooper(Fo):
    """
    Cylindrical heat source approximation solution of Cooper (1976)

    Fo = Fourier number (al*t/rb**2)

    Reference:
    ----------
    Louis Lamarche course #3, slide 10, Eq. 3.3
    """

    Ca = 1.128379
    C0 = -0.5
    C1 = 0.2756227
    C2 = -0.1499385
    C3 = 0.0617932
    C4 = -0.01508767
    C5 = 0.001566857

    if Fo <= 6.124633:
        G = Fo**(0.5)/(2*pi) * (Ca + C0*Fo**(0.5) + C1*Fo + C2*Fo**(1.5) +
                                C3*Fo**(2) + C4*Fo**(2.5) + C5*Fo**(3))
    else:
        euler = 0.57721566490153286060651209008240243104215933593992
        z = log(4*Fo/exp(euler))
        G = (2*z*(8*Fo*(1+2*Fo)-1.0-3*z)+16*Fo+pi**2+3)/(128.0*pi*Fo**2)

    return G


def g_function_bernier(Fo, rba):
    """
    Bernier (2001) approximation of Ingersol et al. (1954) G-function of
    the cylindrical heat source problem.

    Fo = Fourier number
    rba = ratio of the radius where the temperature is calculated over the
          external radius of the borehole.

    Reference:
    ----------
    Bernier, M.A., 2001. Ground-Coupled Heat Pump System Simulation.
        ASHRAE Transactions, 107(1):1-12.
    """

    A = array([[-0.89129, 0.36081, -0.05508, 3.59617e-3, 0, 0],
               [-1.454099, 0.8993336, -0.311928, 0.061119, -0.00478046, 0],
               [-3.007691, 2.256059, -0.7928093, 0.134293, -0.00858244, 0],
               [-9.141771, 11.7025, -7.09574, 2.269837, -0.3669166, 0.023587]])
    if rba == 1:
        a = A[0, :]
    elif rba == 2:
        a = A[1, :]
    elif rba == 5:
        a = A[2, :]
    elif rba == 10:
        a = A[3, :]
    else:
        raise Exception('Second input argument must be either 1, 2, 5 or 10')

    x = log10(Fo)
    arg = a[0] + a[1]*x + a[2]*x**2 + a[3]*x**3 + a[4]*x**4 + a[5]*x**5
    G = 10**arg

    return G


def G_function_ILS(Fo):
    """
    Infinite Line Source solution

    Fo = Fourier number

    Reference:
    ----------
    Lamarche and Beauchamp (2007). A new contribution to the finite
        line-source model for geothermal boreholes. Energy and
        Building, 39:188-198.
    """

    if Fo == 0:
        G_ILS = 0
    else:
        G_ILS = 1/(4*pi) * special.expn(1, 1/(4*Fo))  # Exponential integral

    return G_ILS


def G_function_FLS(Fo, rba, dba=0):
    """
    Finite Line Source solution

    rba  = rb/H
    dba = D/H : Buried depth to borehole ratio
    Fo = Fourier number (al*t/rb^2)

    Reference:
    ----------
    Lamarche and Beauchamp (2007). A new contribution to the finite
        line-source model for geothermal boreholes. Energy and
        Building, 39:188-198.

    Cimmino M., Bernier M. and Adams F. (2013) A contribution towards the
        determination of g-functions using the finite line source. Applied
        Thermal Engineering, 51:401-412.
    """

    if Fo == 0:
        G_FLS = 0
    else:
        tt = 9*rba*rba*Fo
        G_FLS = g_function_fls(tt, rba, dba)/(2*pi)

    return G_FLS


def g_function_fls(tt, rba, dba=0):
    if tt == 0:
        y = 0
    else:
        beta = 3/(2*sqrt(tt))

        # ---- A-Integral ---- #

        am = sqrt(rba**2 + 1)

        DA1 = am*erfc(beta*am) - rba*erfc(beta*rba)
        DA2 = (exp(-beta**2*am**2) - exp(-beta**2*rba**2))/(beta*sqrt(pi))
        DA = DA1 - DA2

        A = integrate.quad(afunc, rba, am, args=(rba, beta))[0] - DA

        # ---- B-Integral ---- #

        a1 = sqrt(rba**2 + (1+2*dba)**2)
        a2 = sqrt(rba**2 + (2+2*dba)**2)
        a3 = sqrt(rba**2 + 4*dba**2)

        DB1 = a1*erfc(beta*a1) - 0.5*(a3*erfc(beta*a3) + a2*erfc(beta*a2))
        DB2 = (exp(-beta**2*a1**2) -
               0.5*(exp(-beta**2*a3**2) + exp(-beta**2*a2**2))
               )/(beta*sqrt(pi))
        DB = DB1 - DB2

        y2b = integrate.quad(afunc, a1, a2, args=(rba, beta))[0]
        y3b = integrate.quad(afunc, a3, a1, args=(rba, beta))[0]
        B = (1+dba)*y2b - dba*y3b + DB

        y = A-B
    return y


def afunc(z, rba, gam):
    return erfc(gam*z)/sqrt(z**2-rba**2)


def test_module():

    rb = 0.15/2   # borehole radius in m
    ks = 2.5                # soil thermal conductivity in W/m.K
    Cps = 2.256e6           # vol. heat capacity of the soil in J/m³.K
    als = ks/Cps * 3600*24  # soil thermal diffusivity in m2/day
    To = 5.0      # soil temperature in Celcius
    qb = -40.0    # radial flux at borehole wall in W/m
    H = 5        # borehole length
    D = 0         # buried depth

    for t in [8760, 240/24, 30/60/24]:  # time in days

        print('-'*25)
        print('for t = %0.1f hour(s) :' % (t*24))

        Fo = als*t/rb**2  # Fourier number

        G1 = G_function_ICS(Fo, 1)          # LSC (Ingersol et al., 1954)
        G4 = G_function_ILS(Fo)             # ILS: Infinite Line Source
        G5 = G_function_FLS(Fo, rb/H, D/H)  # FLS: Finite Line Source

        print(u'\nT0 = %0.3f ºC' % To)
        print(u'G_LSC = %0.8f' % G1)
        print(u'G_ILS = %0.8f' % G4)
        print(u'G_FLS = %0.8f' % G5)
        print()

        Tb1 = To - qb*G1/ks
        Tb4 = To - qb*G4/ks
        Tb5 = To - qb*G5/ks

        print(u'G_LSC: Tb = %0.5f ºC' % Tb1)
        print(u'G_ILS: Tb = %0.5f ºC' % Tb4)
        print(u'G_FLS: Tb = %0.5f ºC' % Tb5)

if __name__ == '__main__':
    test_module()
