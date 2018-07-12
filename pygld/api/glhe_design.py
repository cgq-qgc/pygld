# -*- coding: utf-8 -*-

# Copyright (c) PyGLD Project Contributors
# https://github.com/jnsebgosselin/pygld
#
# This is part of PyGLD (Python Ground-Loop Designer).
# Licensed under the terms of the MIT License.

from __future__ import division, unicode_literals

from copy import copy
import numpy as np
from numpy import pi
from time import time, strftime

# Local Imports :

if __name__ == '__main__':
    import sys
    from os.path import dirname, realpath
    root = dirname(dirname(realpath(__file__)))
    sys.path.append(root)

import geothermie as gtm


# =============================================================================


class OpenGLDesignerBase(object):
    def __init__(self):

        # Ground thermal loads in kW (+ for cooling, - for heating) :

        self.qa = -800
        self.qm = {'cooling': 7500, 'heating': -7500}
        self.qh = {'cooling': 15500, 'heating': -15500}

        self.tph = {'cooling': 6, 'heating': 3}  # peak load duration in h

        self.Nbh = 3  # Number of boreholes in the borefield
        self.dbb = 6  # linear distance between two adjacent boreholes
        self.xcb, self.ycb = gtm.design_borefield(3, 6, mode='sequence')

        self.Tm = {'cooling': 28, 'heating': 0}

        # ---- Borehole geometry ----

        self.btype = '1U-Pipe'

        dpi, dpo = gtm.calcul_sdr(1.25, 11)
        self.rb = 0.1542/2  # borehole radius in m
        self.rpi = [dpi/2] * 2  # inner radius of the pipes in m
        self.rpo = [dpo/2] * 2  # outer radius of the pipes in m

        self.zcp = [0.05, -0.05]
        # cartesian coordinate of the centers of the pipes in complex numbers,
        # where zcp = xcp + ycp*j

        # ---- Thermal properties ----

        self.kpipe = [0.4, 0.4]  # thermal conductivity of the pipes in W/m K

        self.Tg = 15  # undisturbed ground temperature in deg C

        self.Cgrt = 2.256e6  # vol. heat capacity of the grout in J/m³.K
        self.kgrt = 1.40  # grout thermal conductivity in W/m K

        self.Cpg = 2.256e6   # vol. heat capacity of soil in J/m³.K
        self.kgnd = 3.0  # thermal conductivity of soil in W/m.K

        # ---- Fluid properties ----

        # total volumetric flow in borehole in m3/s
        self.Vftot = {'cooling': 0.000650, 'heating': 0.000650}
        self.fluid = 'prop_glycol'
        self.fr = 0.25

        # ---- Advanced Settings ----

        self.Jp = 4  # Multipoles order
        self.Gfunc = {'Gf': 'ICS', 'G1': 'ICS', 'G2': 'STS'}

        # ICS: Infinite cylindrical heat source solution
        # ILS: Infinite Line Source solution
        # FLS: Finite Line Source solution
        # STS: Short-term analysis solution (Lamarche)

        # ---- Results ----

        self.__Ltot = None
        self.Rb2d = None
        self.Rb3d = None
        self.Ra = None
        self.dp = None

        # ---- Costs ----

        self.costs = {}
        self.costs['pipes'] = 0         # $/kg
        self.costs['grout'] = 0         # $/L
        self.costs['fluid'] = 0         # $/L
        self.costs['spacer'] = 0        # $/clips
        self.costs['drilling'] = 0      # $/m
        self.costs['installation'] = 0  # $/bore

        self.__initAttr__()

    def __initAttr__(self):  # attributes that are not linked with the UI
        self.ta = 10  # Annual average load duration in years
        self.rho_hdpe = 950  # density of HDPE material in kg/m3

    # =========================================================================

    @property                                                   # Borehole type
    def btype(self):
        return self.__btype

    @btype.setter
    def btype(self, x):
        if x in ['1U-Pipe', '2U-Pipe', 'Coaxial']:
            self.__btype = x
        else:
            raise NameError("Entry for btype must either be '1U-Pipe', "
                            "'2U-Pipe', or 'Coaxial'")

    # ============================================================ Volumes ====

    @property                         # Volume of fluid per unit length in m3/m
    def Afluid(self):
        if self.btype in ['1U-Pipe', '2U-Pipe']:
            Afluid = 0
            for rp in self.rpi:
                Afluid += pi*rp**2
        elif self.btype == 'Coaxial':
            rp2i = self.rpi[0]
            rp1o = self.rpo[1]
            rp1i = self.rpi[1]
            if rp2i <= rp1o:
                raise ValueError('Outside pipe inner radius is smaller'
                                 'than inside pipe outside radius.')
            Afluid = pi*(rp2i**2 - rp1o**2 + rp1i**2)

        return Afluid

    @property                         # Volume of grout per unit length in m3/m
    def Agrout(self):
        if self.btype in ['1U-Pipe', '2U-Pipe']:
            Abore = pi*self.rb**2
            Apipes = 0
            for rp in self.rpo:
                Apipes += pi*rp**2
            Agrout = Abore - Apipes
        elif self.btype == 'Coaxial':
            Agrout = pi*(self.rb**2 - self.rpo[0]**2)
            if self.rpo[0] <= self.rpi[1]:
                raise ValueError('Outside pipe inner radius is smaller'
                                 'than inside pipe outside radius.')

        return Agrout

    @property                 # Volume of pipe material per unit length in m3/m
    def Apipe(self):
        x = 0
        for rpi, rpo in zip(self.rpi, self.rpo):
            x = x + pi*(rpo**2 - rpi**2)

        return x

    @property
    def mpipe(self):            # Mass of pipe material per unit length in kg/m
        return self.Apipe * self.rho_hdpe

    # -------------------------------------------------------------------------

    @property                                      # Total GLHE borehole length
    def Ltot(self):
        return self.__Ltot

    @property                                            # Total length of pipe
    def Lpipe(self):
        if self.btype == '1U-Pipe':
            lp = self.Ltot*2
        elif self.btype == '2U-Pipe':
            lp = self.Ltot*4
        elif self.btype == 'Coaxial':
            lp = self.Ltot*2

        return lp

    # =========================================================================

    @property                              # soil thermal diffusivity in m²/day
    def dsoil(self):
        return self.kgnd/self.Cpg * 3600 * 24

    @property                             # grout thermal diffusivity in m²/day
    def dgrout(self):
        return self.kgrt/self.Cgrt * 3600 * 24

    @property                                # Volumetric flow in each borehole
    def Vfbh(self):
        x = {}
        for mode in ['heating', 'cooling']:
            x[mode] = self.Vftot[mode]/self.Nbh
        return x

    @property                                    # Volumetric flow in each pipe
    def Vfp(self):
        x = {}
        for mode in ['heating', 'cooling']:
            if self.btype == '2U-Pipe':
                x[mode] = self.Vfbh[mode]/2
            else:
                x[mode] = self.Vfbh[mode]
        return x

    # ============================================== head loss in the loop ====

    def calcul_dp_loop(self, mode):
        hcfluid = gtm.HCFluid(self.fluid, self.Tm[mode], self.fr)
        V = self.Vfp[mode]
        rho = hcfluid.rho
        mu = hcfluid.mu
        if self.btype in ['1U-Pipe', '2U-Pipe']:
            di = 2*self.rpi[0]  # we assume all pipes are the same size
            dp = 2*gtm.calcul_dp_pipe(V, di, rho, mu)
        elif self.btype == 'Coaxial':
            # pressure drop in the inner pipe:
            di = 2 * np.min(self.rpi)
            dp = gtm.calcul_dp_pipe(V, di, rho, mu)
            # pressure drop in the annulus:
            d1 = 2 * np.max(self.rpi)
            d2 = 2 * np.min(self.rpo)
            dp = dp + gtm.calcul_dp_ann(V, d1, d2, rho, mu)

        return dp

    # ================================================== Reynolds Number ======

    def calcul_Re(self, mode):

        # ------------------------------------------- Fluid properties --------

        hcfluid = gtm.HCFluid(self.fluid, self.Tm[mode], self.fr)
        rhof = hcfluid.rho
        muf = hcfluid.mu

        # --------------------------------------- U-Pipe or Inner Pipe --------

        ri = self.rpi[1]

        A = pi*ri**2      # total flow cross-section in m^3
        vf = self.Vfbh[mode]/A  # mean fluid velocity in m/s
        if self.btype == '2U-Pipe':
            vf = vf/2

        Re1 = (2*ri) * vf * rhof/muf  # Reynolds number

        # ---------------------------------------------------- Annulus --------

        ri = self.rpo[0]
        ro = self.rpi[1]
        if ri < ro:
            raise ValueError('ri < ro')

        dh = 2*(ro-ri)            # hydraulic diameter of the annular duct
        A = pi*(ro**2 - ri**2)    # flow cross-section in m³
        vf = self.Vfp[mode]/A     # mean fluid velocity in the conduit m/s
        Re2 = vf * dh * rhof/muf  # Reynolds number

        return Re1, Re2

    # ======================================== Borehole Thermal Resistance ====

    def calcul_Rb(self, mode):

        # Fluid properties:

        hcfluid = gtm.HCFluid(self.fluid, self.Tm[mode], self.fr)
        rhof = hcfluid.rho
        kf = hcfluid.k
        Pr = hcfluid.Pr
        muf = hcfluid.mu
        Vfp = self.Vfp[mode]

        R = {}
        if self.btype in ['1U-Pipe', '2U-Pipe']:
            if self.btype == '1U-Pipe':
                Np = 2
            elif self.btype == '2U-Pipe':
                Np = 4

            # Properties of the pipes (all pipes are assumed identical) :

            kp = self.kpipe[0]
            rpi = self.rpi[0]
            rpo = self.rpo[0]

            # Conductive thermal resistance of a single pipe :

            Rcond = gtm.calcul_Rcond_pipe(kp, rpi, rpo)

            # Convective thermal resistance of the fluid :

            Rconv, hf = gtm.calcul_Rconv_pipe(rpi, Vfp, muf, kf, rhof, Pr)

            # Thermal resistance of a single pipe:

            Rp = Rconv + Rcond
            print('Rp = %0.6f m.K/W' % Rp)

            # Thermal resistance of the borehole:

            Rb, Ra = gtm.calcul_Rb_multipoles(self.kgrt, self.kgnd, self.rb,
                                              rpo, Rp, self.Jp, self.zcp)

            R['Rp_cond'] = Rcond/Np
            R['Rp_conv'] = Rconv/Np
            R['Rp'] = Rp/Np

            R['Rb'] = Rb
            R['Ra'] = Ra
            R['Rg'] = Rb - Rp/Np
        elif self.btype == 'Coaxial':
            RR = gtm.calcul_Rb_coaxial(self.rb, self.kpipe, self.rpi, self.rpo,
                                       self.kgrt, Vfp, muf, kf, rhof, Pr)

            R['Rb'] = RR[0]
            R['Ra'] = RR[1]
            R['Rp_cond'] = RR[2]
            R['Rp_conv'] = RR[3]
            R['Rp'] = RR[2] + RR[3]
            R['Rg'] = RR[4]

        print()
        print('Rb = %0.5f m.K/W' % R['Rb'])
        print('Ra = %0.5f m.K/W' % R['Ra'])
        print('Rg = %0.5f m.K/W' % R['Rg'])
        print('Rp = %0.5f m.K/W + %0.5f m.K/W = %0.5f m.K/W' %
              (R['Rp_cond'], R['Rp_conv'], R['Rp']))
        print()

        return R

    # ======================================================== G-functions ====

    def calcul_Gf(self, mode, H=None):

        # Times in days :
        tph = self.tph[mode]/24
        t1m = 365.25/12
        t10y = self.ta * 365.25

        # Fourier number :
        Fof = self.dsoil * (t10y + t1m + tph)/self.rb**2

        # G-function :
        if self.Gfunc['Gf'] == 'ICS':
            Gf = gtm.G_function_ICS(Fof, 1)
        elif self.Gfunc['Gf'] == 'ILS':
            Gf = gtm.G_function_ILS(Fof)
        elif self.Gfunc['Gf'] == 'FLS':
            if H is None:
                Gf = gtm.G_function_ICS(Fof, 1)
            else:
                Gf = gtm.G_function_FLS(Fof, self.rb/H)

        return Gf

    # -------------------------------------------------------------------------

    def calcul_G1(self, mode, H=None):

        # Times in days:
        tph = self.tph[mode]/24
        t1m = 365.25/12

        # Fourier number:
        Fo1 = self.dsoil * (t1m + tph)/self.rb**2

        # G-function :
        if self.Gfunc['G1'] == 'ICS':
            G1 = gtm.G_function_ICS(Fo1, 1)
        elif self.Gfunc['G1'] == 'ILS':
            G1 = gtm.G_function_ILS(Fo1)
        elif self.Gfunc['G1'] == 'FLS':
            if H is None:
                G1 = gtm.G_function_ICS(Fo1, 1)
            else:
                G1 = gtm.G_function_FLS(Fo1, self.rb/H)

        return G1

    # -------------------------------------------------------------------------

    def calcul_G2(self, mode, Rg, Rp, Rp_cond):
        if self.Gfunc['G2'] == 'ICS':
            Fo2 = self.dsoil * self.tph[mode]/24/self.rb**2
            G2 = gtm.G_function_ICS(Fo2, 1)
        elif self.Gfunc['G2'] == 'ILS':
            Fo2 = self.dsoil * self.tph[mode]/24/self.rb**2
            G2 = gtm.G_function_ILS(Fo2)
        elif self.Gfunc['G2'] == 'STS':
            rb = self.rb        # borehole radius in m

            kg = self.kgrt      # thermal cond. of the grout in W/m.s
            Cpg = self.Cgrt     # volumetric heat capacity of the grout
            alg = self.dgrout   # thermal diff. of the grout in m2/day

            ks = self.kgnd      # thermal cond. of the soil in W/m.s
            als = self.dsoil    # thermal diff. of the soil in m2/day

            Af = self.Afluid
            Ag = self.Agrout
            kp = self.kpipe[0]  # pipe thermal conductivity
            Cpf = gtm.HCFluid(self.fluid, self.Tm[mode], self.fr).Cp

            # ---- calcul equivalent properties ----

            re = rb * np.exp(-2*pi*kg*Rg)             # eq.62
            Cpg_eq = (Ag * Cpg) / (pi*(rb**2-re**2))  # eq.63

            ri = re / np.exp(2*pi*kp*Rp_cond)         # eq.59
            Cpf_eq = (Af * Cpf) / (pi*ri**2)          # eq.60

            # re is the equivalent outside radius
            # ri is the equivalent inner radius
            # Cpg_eq is the equivalent grout vol. heat capacity
            # Cpf_eq is the equivalent fluid vol. heat capacity

            # ---- calcul G-function ----

            tt = alg * self.tph[mode]/24/re**2    # eq.13
            rt = rb/re                            # eq.13
            gam = (alg/als)**0.5                  # eq.13
            kt = ks/kg                            # eq.13

            Rt = 2*pi*kg*Rp                       # eq.21
            v = 1/2 * (ri/re)**2 * Cpf_eq/Cpg_eq  # eq.21

            G2 = gtm.G_function_lamarche(tt, rt, kt, gam, Rt, v)
        else:
            raise NameError("'mode' must be either 'ICS', 'ILS' or 'STS' "
                            " for G2.")

        return G2

    # ====================================================== Calcul Ltot ======

    def calcul_Ltot(self):
        Ltot = {'cooling': 0, 'heating': 0}
        Tp = {'cooling': 0, 'heating': 0}
        Rb3d = {'cooling': 0, 'heating': 0}
        msg = {'cooling': '', 'heating': ''}  # message to GUI
        for mode in ['cooling', 'heating']:

            t1 = time()

            print('-'*50)
            print(mode)
            print('-'*50)
            print()
            print('qa = %0.10f kW' % (self.qa/1000))
            print('qm = %0.10f kW' % (self.qm[mode]/1000))
            print('qh = %0.10f kW' % (self.qh[mode]/1000))
            print()

            # ------------------------------------------- Fluid properties ----

            hcfluid = gtm.HCFluid(self.fluid, self.Tm[mode], self.fr)
            if self.Tm[mode] < hcfluid.Tfp:
                self.__Ltot = None
                self.Rb2d, self.Ra = None, None
                self.Rb3d = None
                self.Tp = None
                return None, None, None, None

            Cpf = hcfluid.Cp
            Vf = self.Vfbh[mode]  # volumetric flow in m3/s

            # ---------------------- Effective borehole thermal resistance ----

            R = self.calcul_Rb(mode)
            Rb2d, Ra = R['Rb'], R['Ra']

            # ----------------------- Effective ground thermal resistances ----

            Gf = self.calcul_Gf(mode)
            G1 = self.calcul_G1(mode)
            G2 = self.calcul_G2(mode, R['Rg'], R['Rp'], R['Rp_cond'])

            print('R6h = %0.5f m.K/W' % (G2/self.kgnd))
            print('R1m = %0.5f m.K/W' % ((G1 - G2)/self.kgnd))
            print('R10y = %0.5f m.K/W' % ((Gf - G1)/self.kgnd))
            print()
            print('Gf =', Gf)
            print('G1 =', G1)
            print('G2 =', G2)

            # ----------------------------------- Total length calculation ----

            qa = self.qa
            qm = self.qm[mode]
            qh = self.qh[mode]
            Tm, Tg = self.Tm[mode], self.Tg
            ks, kg = self.kgnd, self.kgrt
            if self.Gfunc['G2'] in ['ICS', 'ILS']:
                L2d = (qa/ks*(Gf - G1) +
                       qm/ks*(G1 - G2) +
                       qh/ks*G2 + qh*Rb2d) / (Tm - Tg)
            elif self.Gfunc['G2'] == 'STS':
                L2d = (qa/ks*(Gf - G1) +
                       qm/ks*G1 + qm*Rb2d - qm/kg*G2 +
                       qh/kg*G2) / (Tm - Tg)

            print('Rb2d = %0.5f' % Rb2d)
            print('Ltot2d = %0.2f m\n' % L2d)

            if L2d <= 0:
                print('No solution in %s mode due to unbalanced loads' % mode)
                continue

            Lold = copy(L2d)
            icount = 1
            tf = 10*365.25  # long time scale in days
            while True:
                if icount > 101:
                    msg[mode] = ('Unable to converge '
                                 'after %d iterations' % icount)
                    print('\n' + msg[mode])
                    break

                # ------------------------------------------ Calcul ASHRAE ----

                Rb3d[mode] = gtm.compute_Rb3D(Rb2d, Ra, Vf, Lold/self.Nbh, Cpf)
                if Rb3d[mode] > 10**12:
                    print('Unable to converge due to Rb3D calculation in'
                          ' %s mode' % mode)
                    Ltot[mode] = np.inf
                    break

                Tp[mode] = gtm.calcul_Tp_fls(self.xcb, self.ycb, self.qa/Lold,
                                             self.dsoil, self.kgnd, tf,
                                             Lold/self.Nbh)
                if self.Gfunc['Gf'] == 'FLS':
                    Gf = self.calcul_Gf(mode, Lold/self.Nbh)
                if self.Gfunc['G1'] == 'FLS':
                    G1 = self.calcul_G1(mode, Lold/self.Nbh)
                if self.Gfunc['G2'] in ['ICS', 'ILS']:
                    Ltot[mode] = (qa/ks*(Gf - G1) +
                                  qm/ks*(G1 - G2) +
                                  qh/ks*G2 + qh*Rb3d[mode]
                                  ) / (Tm - Tg - Tp[mode])
                elif self.Gfunc['G2'] == 'STS':
                    Ltot[mode] = (qa/ks*(Gf - G1) +
                                  qm/ks*G1 + qm*Rb3d[mode] - qm/kg*G2 +
                                  qh/kg*G2
                                  ) / (Tm - Tg - Tp[mode])

                Ltot[mode] = np.max([np.round(Ltot[mode]), 0])

                print('iter %d : Rb=%0.3f ; Tp=%0.2f ; L=%d' %
                      (icount, Rb3d[mode], Tp[mode], Ltot[mode]))

                # ----------------------------------- PostProcess Results -----

                if Ltot[mode] <= 0:
                    msg[mode] = 'Unable to converge: '
                    if mode == 'cooling':
                        if Tp[mode] < 0:
                            msg[mode] += 'system is heating-dominated'
                        else:
                            msg[mode] += 'Tp is too large.\nTry increasing '
                            msg[mode] += 'the distance between the boreholes.'
                    elif mode == 'heating':
                        if Tp[mode] > 0:
                            msg[mode] + 'system is cooling-dominated'
                        else:
                            msg[mode] += 'Tp is too large.\nTry increasing '
                            msg[mode] += 'the distance between the boreholes.'
                    print('\n' + msg[mode])
                    break

                if abs(Ltot[mode]-Lold) <= 1:
                    t2 = time()
                    msg[mode] = ('Converged after %d iterations'
                                 ' in %0.3f sec' % (icount, t2-t1))
                    print('\n' + msg[mode])
                    print('Tp = ', Tp[mode])
                    print('Tm = ', self.Tm[mode])
                    print('Rb2d = ', Rb2d)
                    print('Rb3d = ', Rb3d[mode])
                    print('Ltot = %0.2f m' % Ltot[mode])
                    break
                else:
                    icount += 1
                    Lold = copy(Ltot[mode])

        print('-'*50)

        # Determine if system is cooling or heating dominated (dm):

        if Ltot['cooling'] > Ltot['heating']:
            dm = 'cooling'
        elif Ltot['cooling'] < Ltot['heating']:
            dm = 'heating'
        elif Ltot['cooling'] == Ltot['heating']:
            if self.qa > 0:
                dm = 'cooling'
            else:
                dm = 'heating'

        msg0 = 'GLHE design length is %s-controlled' % dm
        msg[dm] = msg0 + '\n' + msg[dm]
        print(msg0)
        print('Tp = ', Tp[dm])
        print('Ltot = %0.2f m' % Ltot[dm])
        print()

        # store results in class variables :

        self.__Ltot = Ltot[dm]
        self.Lbh = self.Ltot/self.Nbh

        R = self.calcul_Rb(dm)
        self.Rb2d = R['Rb']
        self.Ra = R['Ra']
        self.Rb3d = Rb3d[dm]

        # pressure drop in loop (kpa):
        self.dp = self.calcul_dp_loop(dm)*Ltot[dm]/self.Nbh/1000

        self.Tp = Tp[dm]
        self.calcul_cost()

        print('-'*50)

        return Ltot[dm], Tp[dm], dm, Rb3d[dm]

    # ============================================================== Costs ====

    def calcul_cost(self):
        cost_pipe = self.mpipe*self.costs['pipes'] * self.Ltot
        cost_fluid = self.Afluid*self.costs['fluid']*1000 * self.Ltot * self.fr
        cost_grout = self.Agrout*self.costs['grout']*1000 * self.Ltot
        if self.btype in ['1U-Pipe', '2U-Pipe']:
            cost_spacer = self.costs['spacer']*np.ceil(self.Lbh/10)*self.Nbh
        else:
            cost_spacer = 0

        cost_drilling = self.costs['drilling'] * self.Ltot
        cost_install = self.costs['installation'] * self.Nbh

        self.borefield_cost = (cost_pipe + cost_fluid + cost_grout +
                               cost_spacer + cost_drilling + cost_install)

#        print(cost_pipe, cost_fluid, cost_grout, cost_spacer, cost_drilling,
#              cost_install)

        return self.borefield_cost


if __name__ == '__main__':
    gld = OpenGLDesignerBase()

    # ---- Ground loads ----

    gld.qa = 1330
    gld.qm = {'cooling': 28598, 'heating': -21379}
    gld.qh = {'cooling': 134633, 'heating': -170781}
    gld.tph = {'cooling': 7, 'heating': 3}

    gld.Nbh = 11  # Number of boreholes in the borefield
    gld.dbb = 6  # linear distance between two adjacent boreholes
    gld.xcb, gld.ycb = gtm.design_borefield(12, 6, mode='sequence')

    gld.Tm = {'cooling': 29.2427058276, 'heating': -1.58684925111}

    # ---- Borehole geometry ----

    gld.btype = '1U-Pipe'
    gld.rb = 0.15240/2                # borehole radius in m
    gld.rpi = [0.03399/2, 0.03399/2]  # inner radius of the pipes in m
    gld.rpo = [0.04216/2, 0.03399/2]  # outer radius of the pipes in m
    gld.zcp = [0.05, -0.05]           # coordinates center of the pipes in m

    # ---- Thermal properties ----

    gld.kpipe = [0.7]*2  # thermal conductivity of the pipes in W/m K

    gld.Tg = 11  # undisturbed ground temperature in deg C

    gld.Cgrt = 3.9e6  # vol. heat capacity of the grout in J/m³.K
    gld.kgrt = 1.90   # grout thermal conductivity in W/m K

    gld.Cpg = 2.350e6   # vol. heat capacity of soil in J/m³.K
    gld.kgnd = 3.0      # thermal conductivity of soil in W/m.K

    # ---- Fluid properties ----

    # total volumetric flow in borehole in m3/s
    gld.Vftot = {'cooling': 13.5/1000, 'heating': 13.5/1000}
    gld.fluid = 'prop_glycol'
    gld.fr = 0.25

    # ---- Advanced Settings ----

    gld.Jp = 4  # Multipoles order
    gld.Gfunc = {'Gf': 'FLS',
                 'G1': 'FLS',
                 'G2': 'STS'}

    # ICS: Infinite cylindrical heat source solution
    # ILS: Infinite Line Source solution
    # FLS: Finite Line Source solution
    # STS: Short-term analysis solution (Lamarche)

    gld.calcul_Ltot()
