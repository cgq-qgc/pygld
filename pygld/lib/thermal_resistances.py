# -*- coding: utf-8 -*-

# Copyright (c) PyGLD Project Contributors
# https://github.com/jnsebgosselin/pygld
#
# This is part of PyGLD (Python Ground-Loop Designer).
# Licensed under the terms of the MIT License.

import numpy as np
from numpy import (pi, log, conj, zeros, diag)


# ---- Head Loss in Pipes

def calcul_dp_pipe(V, di, rho, mu):
    """
    Calcul the hydraulic head loss per unit of length in a circular pipe.

    V: volumetric flow rate in the pipe in m3/s
    di: inside diameter of the pipe in m
    rho: fluid density in kg/m3
    mu: fluid viscosity in Pa/s
    """
    # Calcul the Reynolds number :

    A = np.pi * di**2/4
    w = V / A
    Re = rho * w / mu * di

    # Calcul Friction factor :

    if Re < 2300:
        # The flow is laminar.
        f_pipe = 64/Re
    else:
        # The flow is turbulent.
        f_pipe = 0.3164 * Re**(-1/4)

    # Calcul the pressure drop in Pa/m :

    dp = f_pipe * (1/di) * (rho*w**2/2)

    return dp


def calcul_dp_ann(V, d1, d2, rho, mu):
    """
    Calcul the hydraulic head loss per unit of length in the annulus for
    the concentric borehole design

    d1: outer diameter of the annulus
    d2: inner diameter of the annulus
    V: Volumetric flow inside the annulus in m3/s
    rho: fluid density in kg/m3
    mu: fluid viscosity in Pa/s
    """

    if d2 >= d1:
        raise ValueError('d1 must be greater than d2.')

    # ---- Calcul the Reynolds number :

    dh = d1 - d2

    A = np.pi * (d1**2 - d2**2)/4
    w = V / A
    Re = rho * w / mu * dh

    # ---- Calcul Friction factor :

    a = d2/d1
    Ree = Re * ((1+a**2)*np.log(a) + (1-a**2)) / ((1-a**2)*np.log(a))
    if Ree < 2300:
        # The flow is laminar.
        f = 64/Ree
    else:
        # The flow is turbulent.
        f = (1.8*np.log10(Ree)-1.5)**-2

    # ---- Calcul the pressure drop in Pa/m :

    dp = f * (1/dh) * (rho*w**2/2)

    return dp


# ---- Nusselt Number


def calcul_Nu_pipe(Re, Pr):
    """
    Compute mean Nusselt number for a flow in a circular pipe.

    Num = calcul_Nu_tube(Re, Pr)

    Re = Reynolds number formed with the tube diameter.
    Pr = Prandtl number
    Num = mean Nusselt number

    References:
    Baehr H.D. and K. Stephan, 2006. Heat and Mass Transfer. Second ed.
    Springer-Verlag, Berlin, Heidelberd. p.370-372.

    - Valid for thermally and hydrodynamic fully developped flow
    - Valid for small temperature differences between the flow and the wall
    """

    if Re < 2300:  # Laminar flow
        Num = 3.657  # Asymptotic value assuming constant wall temperature

    elif Re >= 2300 and Re < 10**4:  # Transition zone
        Num = 0.037 * (Re**0.75 - 180) * Pr**0.42

    elif Re >= 10**4:  # Turbulent flow
        f = (0.78*log(Re) - 1.5)**-2

        num = (f/8) * Re * Pr
        den = 1 + 12.7*np.sqrt(f/8) * (Pr**(2/3.)-1)
        Num = num / den

    else:
        print('Error: Reynolds number is nan')
        Num = np.nan

    # Formulation de Hellström (1991)

    # elif Re >= 2300:  # Turbulent + Laminar flow
    #     f = (1.58*np.log(Re) - 3.28)**-2

    #   num = (f/2.) * (Re-1000) * Pr
    #   den = 1 + 12.7*np.sqrt(f/2.)*(Pr**(2/3.) - 1)
    #   Num = num / den

    return Num


def calcul_Nu_annulus(Re, Pr, a):
    """
    Compute mean Nusselt number for a flow in an annular space between two
    concentric circular pipes.

    Num = calcul_Nu_annulus(Re, Pr, a)

    Re = Reynolds number formed with the hydraulic diameter dh = do-di, where
         do and di are, respectively, the outer and inner diameter of the
         annulus.
    Pr = Prandtl number
    a = ratio of the inner and outer diameter of the annulus (a = di/do),
        where 0 >= a >= 1.
    Num = mean Nusselt number

    References:
    Baehr H.D. and K. Stephan, 2006. Heat and Mass Transfer. Second ed.
    Springer-Verlag, Berlin, Heidelberd. p.370-372.

    - Valid for thermally and hydrodynamic fully developped flow
    - Valid for small temperature differences between the flow and the wall
    """

    Nu_tube = calcul_Nu_pipe(Re, Pr)

    if Re < 2300:  # Laminar flow
        # Heat is transferred at both the inner and outer tubes.
        Num = 3.657 + (4 - 0.102/(0.02 + a)) * a**0.04  # Asymptotic value

    elif Re >= 2300:  # Transition + turbulent flow
        # Heat is transferred at the inner tube. The outer tube is insulated.
        Nuii = Nu_tube * (0.86 * a**-0.16)

        # Heat is transferred at the outer tube. The inner tube is insulated.
        Nuoo = Nu_tube * (1 - 0.14*a**0.6)

        # Heat is transferred at both the inner and outer tubes.
        Num = (Nuii + Nuoo) / (1 + a)

    return Num


def calcul_Nu_annulus2(Re, Pr, a):
    """
    Using VDI Waermeatlas inster of Baehr

    Compute mean Nusselt number for a flow in an annular space between two
    concentric circular pipes.

    Num = calcul_Nu_annulus(Re, Pr, a)

    Re = Reynolds number formed with the hydraulic diameter dh = do-di, where
         do and di are, respectively, the outer and inner diameter of the
         annulus.
    Pr = Prandtl number
    a = ratio of the inner and outer diameter of the annulus (a = di/do),
        where 0 >= a >= 1.
    Num = mean Nusselt number


    - Valid for thermally and hydrodynamic fully developped flow
    - Valid for small temperature differences between the flow and the wall
    """

    Nu_tube = calcul_Nu_tube(Re, Pr)

    if Re < 2300:  # Laminar flow
        # Heat is transferred at both the inner and outer tubes.
        Num = 3.66 + (4 - 0.102/(0.02 + a)) * a**0.04  # Asymptotic value

    elif Re >= 2300:  # Transition + turbulent flow
        # Heat is transferred at the inner tube. The outer tube is insulated.
        Nuii = Nu_tube * (0.86 * a**-0.16)

        # Heat is transferred at the outer tube. The inner tube is insulated.
        Nuoo = Nu_tube * (1 - 0.14*a**0.6)

        # Heat is transferred at both the inner and outer tubes.
        Num = (Nuii + Nuoo) / (1 + a)

    return Num


# ---- Pipes Thermal Resistance

def calcul_Rconv_pipe(ri, Vf, muf, kf, rhof, Pr):
    """
    Compute the convective thermal resistance in a tube.

    Rconv, hf = calcul_Rconv_pipe(ri, Vf, Tkf, fluid, fr=0)

    ri = inner radius of the pipe in m
    Vf = volumetric flow of the fluid in the pipe in m3/s

    muf = dynamic viscosity of the fluid in Pa.s
    kf = thermal conductivtiy of the fluid in W/m.K
    rhof = density of the fluid in kg/m3
    Pr = Prandtl number

    Rconv = convective thermal resistance in m.K/W
    hf = convective heat transfer coefficient
    """

    # Calcul fluid average velocity in m/s:

    A = pi*ri**2   # total flow cross-section in m^3
    vf = Vf/A      # mean fluid velocity in m/s

    # Calcul Reynold and Nusselt numbers :

    Re = (2*ri) * vf * rhof / muf  # Reynolds number
    Num = calcul_Nu_pipe(Re, Pr)   # mean Nusselt number

    # Convective heat transfer coefficient :
    hf = Num * kf / (2*ri)

    # Convective heat resistance :
    Rconv = 1. / (2*pi*ri*hf)

    print(('Re = %0.0f ; Nu = %0.1f ; hf = %0.5f W/m².K '
           '; Rconv = %0.6f m.K/W') % (Re, Num, hf, Rconv))

    return Rconv, hf


def calcul_Rcond_pipe(kp, ri, ro):
    """
    Compute the conductive thermal resistance of a circular pipe.

    Rcond = calcul_Rcond_pipe(kp, ri, ro)

    kp = thermal conductivity of the pipe in W/m.K
    ri = inner radius of the pipe in m
    ro = outer radius of the pipe in m
    Rcond = conductive thermal resistance in m.K/W
    """
    return np.log(ro/ri) / (2*pi*kp)


def calcul_Rp(kp, rpi, rpo, Vf, muf, kf, rhof, Pr):
    """
    Compute the conductive and convective thermal resistance of a
    circular pipe.

    Rp = calcul_Rpipe(kp, ri, ro, Vf, Tref, fluid, fr)

    kp = thermal conductivity of the pipe in W/m.K
    rpi = inner radius of the pipe in m
    rpo = outer radius of the pipe in m
    Vf = Volumetric flow of the fluid in m3/s
    Tref = mean fluid temperature in Celcius
    fluid = can be either 'water', 'prop_glycol', or 'ethy_glycol'
    fr = antifreeze volumetric fraction, 0<fr<1
    Rcond = conductive thermal resistance in m.K/W
    Rconv = convective thermal resistance in m.K/W
    Rpipe = conductive + convective thermal resistance in m.K/W
    """

    Rconv, hf = calcul_Rconv_pipe(rpi, Vf, muf, kf, rhof, Pr)
    Rcond = calcul_Rcond_pipe(kp, rpi, rpo)

    Rp = Rconv + Rcond

    print('Rpipe = %0.6f m.K/W' % Rp)

    return Rp


# ---- Coaxial Thermal Resistance

def calcul_hf_annulus(ri, ro, Vf, muf, kf, rhof, Pr):

    # Compute Reynold number and Nusselt number :

    dh = 2*(ro-ri)  # hydraulic diameter of the annular duct
    a = ri / ro     # inner to outer diameter ratio

    A = pi*(ro**2 - ri**2)   # flow cross-section in m³
    vf = Vf/A                # mean fluid velocity in m/s

    Re = vf * dh * rhof / muf   # Reynolds number
    Num = calcul_Nu_annulus(Re, Pr, a)  # mean Nusselt number

    hf = kf * Num / dh

    print()
    if Re > 2300:
        print('The flow in the annulus is turbulent :')
    else:
        print('The flow in the annulus is laminar :')
    print('Re = %0.2f ; Nu = %0.2f ; hf = %0.2f W/m².K' % (Re, Num, hf))

    return hf


def calcul_Rconv_annulus(ri, ro, Vf, muf, kf, rhof, Pr):
    """
    Ri, Ro, hf = calcul_Rconv_annulus(ri, ro, Vf, Tref, fluid, fr=0)

    Vf = volumetric flow of the fluid in m³/s
    ri = inner radius of the annulus in m
    ro = outer radius of the annulus in m
    Ri = convective resistance of the inner wall of the annulus in m.K/W
    Ro = convective resistance of the outer wall of the annulus in m.K/W
    hf = mean heat transfer coefficient in W/m².K

    References:
    Baehr and Stephan, 2006. Heat and Mass Transfer. Second ed. Springer-
        Verlag, Berlin, Heidelberd. pp370-372.
    """

    # Mean heat transfer coefficient in W/m².K :

    hf = calcul_hf_annulus(ri, ro, Vf, muf, kf, rhof, Pr)

    # Convective thermal resistance at inner and outer wall of annulus :

    Ri = 1. / (hf*2*pi*ri)
    Ro = 1. / (hf*2*pi*ro)

    print(('Rconv_in = %0.6f m.K/W ; Rconv_out = %0.6f m.K/W') % (Ri, Ro))

    return Ri, Ro, hf


def calcul_Rb_coaxial(rb, kp, rpi, rpo, kgrout, Vf, muf, kf, rhof, Pr):
    """
    Compute fluid to ground thermal resistance of coaxial ground heat
    exchangers.

    rb = borehole radius in m
    kp = [kp1, kp2]
    kp1 = inner pipe thermal conductivity in W/m.K
    kp2 = outer pipe thermal conductivity in W/m.K
    kg = thermal conductivity of the grout in W/m.K

    Vf = fluid flow rate in m³/s

    rpi = inside radius of pipes in m
    rp1i = inside radius of inner pipe in m
    rp1o = outside radius of inner pipe in m

    rp2i = inside radius of outer pipe in m
    rp2o = outside radius of outer pipe in m

    Tref = reference temperature in Kelvin
    fluid = type of fluid (water, prop_glycol, ethy_glycol)
    """

    indx_outpipe = np.argmax(rpo)
    indx_intpipe = 1 - indx_outpipe

    kp1 = kp[indx_intpipe]
    kp2 = kp[indx_outpipe]

    rp1i = rpi[indx_intpipe]
    rp1o = rpo[indx_intpipe]

    rp2i = rpi[indx_outpipe]
    rp2o = rpo[indx_outpipe]

    if rp2i < rp1o:
        raise ValueError('Inner pipe ouside diameter is greater than '
                         'outer pipe inside diameter.')

    # ---- Internal thermal resistance, Ra

    # The thermal resistance between the inner and the outer flow channel
    # (internal thermal resistance) consists of:

    # convective heat transfer resistance between the bulk fluid in
    # the inner flow channel and the inner surface of the inner pipe:
    Rconv_p1_int, _ = calcul_Rconv_pipe(rp1i, Vf, muf, kf, rhof, Pr)

    # conductive thermal resistance of the inner pipe:
    Rcond_p1 = calcul_Rcond_pipe(kp1, rp1i, rp1o)

    # Convective heat transfer resistance between the outer surface of
    # the inner pipe and the bulk fluid in the outer flow channel and
    # convective heat transfer resistance between the bulk fluid in the
    # annular flow channel and the outer surface of the channel.
    Rconv_p1_ext, Rconv_p2, _ = calcul_Rconv_annulus(rp1o, rp2i, Vf, muf,
                                                     kf, rhof, Pr)

    # The internal thermal resistance is:
    Ra = Rconv_p1_int + Rcond_p1 + Rconv_p1_ext

    # ---- External thermal resistance, Rb

    # The thermal resistance between the outer flow channel and the borehole
    # wall is composed of three parts:

    # conductive thermal resistance of the outer pipe:
    Rcond_p2 = calcul_Rcond_pipe(kp2, rp2i, rp2o)

    # conductive thermal resistance of the grout:
    Rgrout = calcul_Rcond_pipe(kgrout, rp2o, rb)

    # The external borehole thermal resistance is:
    Rb = Rconv_p2 + Rcond_p2 + Rgrout

    return Rb, Ra, Rcond_p2, Rconv_p2, Rgrout


# ---- U-pipe Borehole Thermal Resistance

def calcul_Rb_linesource(ki, k, rb, r1, xc):
    # ki : grout thermal conductivity
    # k  : soil thermal conductivity
    # rb : borehole radius
    # r1 : pipe outside radius
    # xc : demi-distance des tuyaux

    sig = (ki-k)/(ki+k)
    b11 = xc/rb
    # b22 = b11
    # b12 = 2*xc/rb

    la1 = rb/r1
    la2 = rb/xc
    # la3 = la2/(2*la1)

    Rb = (log(la1) +
          log(la2/2) +
          sig*log(la2**4/(la2**4-1))
          ) / (4*pi*ki)

    Ra = (log(2*xc/r1) + sig*log((1+b11**2) / (1-b11**2))) / (pi*ki)

    return Rb, Ra


def calcul_Rb_Paul(kg, rb, rp, cas):
    # calcul resistance Redmund
    if cas.lower() == 'a':
        beta0 = 20.10
        beta1 = -0.94447
    elif cas.lower() == 'b':
        beta0 = 17.44
        beta1 = -0.6052
    elif cas.lower() == 'c':
        beta0 = 21.90587
        beta1 = -0.3796
    else:
        print('error, cas = a b or c')

    Sb = beta0 * (rb/rp)**beta1
    Rg = 1/(kg*Sb)

    return Rg


def calcul_Rb_Sharqawi(kg, rb, rp, xc):
    # calcul resistance Redmund

    Rb = 1/(2*pi*kg)*(-1.49*(xc/rb)+0.656*log(rb/rp)+0.436)

    return Rb


def calcul_Rb_multipoles(kb, ks, rb, rp, Rp, Jp, z):

    # kb : grout thermal conductivity
    # ks : soil thermal conductivity
    # rb : borehole radius
    # rp : pipe external radius
    # z  : pipe coordinates in complex notation where z = x + iy
    # Jp : Number of multipoles at each pipes. If Jp == 0, only line sources at
    #      the pipes are used.
    # Rp : pipe thermal resistance, which corresponds to the sum of the thermal
    #      resistance of the pipe wall and the fluid boundary layer.

    beta = 2 * pi * kb * Rp  # dimensionless thermal resistance
    sig = (kb-ks)/(kb+ks)      # thermal conductivity parameter
    r = np.abs(z)            # polar coordinate

    N = len(z)  # Number of pipes within the borehole

    # -------------------------------------------------------------------------
    # Constant component of angular-dependent fluid temperatures
    # (Eq.33 in Claesson and Hellstrom, 2012)

    Ro = np.zeros((N, N))  # Borehole thermal resistance for J=0
    for n in range(N):  # m == n
        Ro[n, n] = (np.log(rb/rp) +
                    beta +
                    sig * np.log(rb**2/(rb**2-r[n]**2)))

    for n in range(N-1):  # m != n
        for m in range(n+1, N):
            rmn = np.abs(z[n] - z[m])
            Ro[n, m] = (np.log(rb/rmn) +
                        sig * np.log(rb**2 / np.abs(rb**2-np.conj(z[n])*z[m])))
            Ro[m, n] = Ro[n, m]

    Ro = Ro / (2*pi*kb)

    # -------------------------------------------------------------------------
    # Components j>=1 of angular-dependent fluid temperatures

    if Jp == 0:
        R = Ro
    else:
        Pold = np.zeros((N, Jp)).astype(complex)
        P = np.zeros((N, Jp)).astype(complex)
        F = np.zeros((N, Jp)).astype(complex)
        R = np.zeros((N, N)).astype(complex)
        for m in range(N):
            q = np.zeros(N)
            q[m] = 1
            compt = 0  # safeguard in case convergence is not reached
            while 1:
                for n in range(N):
                    for ik in range(Jp):
                        F[n, ik] = calcul_F(q, P, rp, rb, kb, sig,
                                            n, ik, Jp, N, z)
                        P[n, ik] = ((-1 + (ik+1)*beta) /
                                    (1 + (ik+1)*beta) *
                                    np.conj(F[n, ik]))

                # Eq.38 in Claesson and Hellstrom (2012)
                err = np.max(np.abs(P - Pold))
                if err < 1e-6:
                    break
                else:
                    Pold = np.copy(P)
                    compt = compt + 1
                    if compt > 100:
                        break

            # second and third sums of Eq.32 in Claesson and Hellstrom (2012)
            for im in range(m, N):
                s1 = 0
                s2 = 0
                for i in range(N):
                    if i != im:
                        for ij in range(Jp):
                            s1 = s1 + P[i, ij]*rp**(ij+1)/(z[im]-z[i])**(ij+1)

                    for ij in range(Jp):
                        num = rp**(ij+1) * conj(z[im])**(ij+1)
                        den = (rb**2 - conj(z[im])*z[i])**(ij+1)
                        s2 = s2 + P[i, ij] * num/den

                R[m, im] = Ro[m, im] + np.real(s1 + sig*s2)
                R[im, m] = R[m, im]

    # -------------------------------------------------------------------------
    # Calculating borehole resistance, Rb

    Rm = np.linalg.inv(R)  # Eq.48 in Claesson and Hellstrom (2012)
    K = zeros((N, N))

    for i in range(N):  # Eq.49 in Claesson and Hellstrom (2012)
        K[i, i] = np.sum(np.real(Rm[i, :]))
    for i in range(N-1):
        for j in range(i+1, N):
            K[i, j] = -np.real(Rm[i, j])
            K[j, i] = K[i, j]

    Kb = np.sum(diag(K))  # Eq.54 in Claesson and Hellstrom (2012)
    Rb = 1./Kb

    # -------------------------------------------------------------------------
    # Calculating internal resistance, Ra

    if N == 2:
        Ka = (K[0, 0] + 2*K[0, 1])/2.0
        Ra = 1./Ka
    elif N == 4:
        Ka1 = (K[0, 0] + 2*K[0, 1] + 2*K[0, 2])
        # Ka2 = (K[0, 0] + 2*K[0, 1] + 2*K[0, 3])
        Ra = 1./Ka1

    # In cases pipe 1 and 2 are not symmetric, Ra must be calculated as:
#    DetR = R[0, 0] * R[1, 1] - R[0, 1]**2
#    K12 = np.real(R[0, 1]/DetR)
#    K1b = np.real((R[1, 1] - R[0, 1])/DetR)
#    K2b = np.real((R[0, 0] - R[0, 1])/DetR)
#    R12 = 1/K12
#    R1b = 1/K1b
#    R2b = 1/K2b
#    R2b = R1b
#    Raj = R12 * (R1b+R2b) / (R1b+R2b+R12)
    #print(K1b, K12, K2b)
    #print(Raj-Ra)

#    Ka = (2*K[0, 0] + K[1, 1] + K[0, 1]) / (K[0, 1]*)

    return Rb, Ra


def calcul_F(q, P, rp, rb, kb, sig, m, km, Jp, N, z):
    """Calcul Eq.34 in Claesson and Hellstrom (2011)."""
    k = km+1
    sa = 0
    for i in range(N):
        if i != m:
            den = 2*pi*kb*k*(z[i]-z[m])**k
            sa = sa + q[i]*rp**k/den
    sc = 0
    for i in range(N):
        if i != m:
            for ij in range(Jp):
                num = P[i, ij]*combinaisons(ij+k, ij)*rp**(ij+1)*(-rp)**k
                den = (z[m] - z[i])**(ij+k+1)
                sc = sc+num/den
    sb = 0
    for i in range(N):
        num = q[i]*rp**k*conj(z[i])**k
        den = 2*pi*kb*k*(rb**2-z[m]*conj(z[i]))**k
        sb = sb + num/den
    sd = 0
    for i in range(N):
        for ij in range(Jp):
            nj = min(ij+1, k)
            sj = 0
            for jp in range(nj+1):
                num1 = (combinaisons(ij+1, jp) *
                        combinaisons(ij+k-jp, ij))
                num2 = rp**(ij+k+1)*z[m]**(ij+1-jp)*conj(z[i])**(k-jp)
                den = (rb**2-z[m]*conj(z[i]))**(k+ij+1-jp)
                sj = sj+num1*num2/den
            sd = sd + conj(P[i, ij])*sj
    F = sa + sig*(sb+sd) + sc

    return F


def combinaisons(n, k):
    if k == 0:
        return 1
    elif k == 1:
        return n
    else:
        p = n-k+1
        for i in range(2, k+1):
            p = p*(n-k+i)/(i*1.0)
        return p


# ---- 3D Borehole thermal resistance

def compute_Rb3D(Rb2d, Ra, Vf, L, Cpf):
    # Ra = Borehole internal resistance in m.K/W
    # L = Borehole length in m
    # Cpf = fluid vol heat capacity in J/(m3⋅K)
    # Vf = volumetric outflow in m3/s

    Rb3d = Rb2d + 1/(3*Ra) * (L/(Vf*Cpf))**2

    return Rb3d


if __name__ == '__main__':
    from pygld.properties import HeatCarrierFluid

    kp = 0.225
    # kp = 0.4
    ri = 1.34/2 * 0.0833333
    ro = 1.66/2 * 0.0833333

    # conductive thermal resistance :

    Rcond =  np.log(ro/ri)/(2*pi*kp)
    print('Rcond =', Rcond)

    # convective thermal resistance :

    hcfluid = HeatCarrierFluid(fluid='water', Tref=20, fr=0)
    Pr = hcfluid.Pr
    kf = hcfluid.k * 0.5781759824
    Re = 2300

    f = (1.58*np.log(Re) - 3.28)**-2
    num = (f/2) * (Re-1000) * Pr
    den = 1 + 12.7*np.sqrt(f/2)*(Pr**(2/3) - 1)
    Num = num/den

    hf = Num*kf/(2*ri)
    print(hf, kf)

    Rconv = 1/(2*pi*ri*hf)
    print('Rconv =', Rconv)

    print('Rtot', (Rcond+Rconv)/2)
