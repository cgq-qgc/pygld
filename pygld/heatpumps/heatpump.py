# -*- coding: utf-8 -*-

# Copyright © 2018 Jean-Sébastien Gosselin
# https://github.com/jnsebgosselin/qwatson
#
# This file is part of PyGLD.
# Licensed under the terms of the GNU General Public License.

# ---- Standard imports

# ---- Third party imports

import numpy as np
from scipy.spatial import Delaunay

# ---- Local imports

from pygld.fluidproperties import HeatCarrierFluid
from pygld.heatpumps.utils import load_heatpump_database
from pygld.heatpumps.plots import plot_fitmodel_eval_from

corrtbl_afreeze = {}

tbl_fluid0 = {}
tbl_fluid0['name'] = 'water'
tbl_fluid0['fr'] = [0]
tbl_fluid0['CAPh'] = [1]
tbl_fluid0['CAPc'] = [1]
tbl_fluid0['COPh'] = [1]
tbl_fluid0['COPc'] = [1]
tbl_fluid0['Wh'] = [1]
tbl_fluid0['Wc'] = [1]
tbl_fluid0['WPD'] = [1]

tbl_fluid1 = {}
tbl_fluid1['name'] = 'prop_glycol'
tbl_fluid1['fr'] = [0, 0.05, 0.15, 0.25]
tbl_fluid1['CAPh'] = np.array([1, 0.989, 0.968, 0.947])
tbl_fluid1['CAPc'] = np.array([1, 0.995, 0.986, 0.978])
tbl_fluid1['Wh'] = np.array([1, 0.997, 0.990, 0.983])
tbl_fluid1['Wc'] = np.array([1, 1.003, 1.009, 1.014])
tbl_fluid1['COPh'] = tbl_fluid1['CAPh']/tbl_fluid1['Wh']
tbl_fluid1['COPc'] = tbl_fluid1['CAPc']/tbl_fluid1['Wc']
tbl_fluid1['WPD'] = [1, 1.070, 1.210, 1.360]

tbl_fluid2 = {}
tbl_fluid2['name'] = 'ethyl_glycol'
tbl_fluid2['fr'] = [0, 0.05, 0.15, 0.25]
tbl_fluid2['CAPh'] = np.array([1, 0.993, 0.980, 0.966])
tbl_fluid2['CAPc'] = np.array([1, 0.998, 0.994, 0.988])
tbl_fluid2['Wh'] = np.array([1, 0.998, 0.994, 0.990])
tbl_fluid2['Wc'] = np.array([1, 1.002, 1.004, 1.008])
tbl_fluid2['COPh'] = tbl_fluid1['CAPh']/tbl_fluid1['Wh']
tbl_fluid2['COPc'] = tbl_fluid1['CAPc']/tbl_fluid1['Wc']
tbl_fluid2['WPD'] = [1, 1.040, 1.120, 1.200]

corrtbl_afreeze[tbl_fluid0['name']] = tbl_fluid0
corrtbl_afreeze[tbl_fluid1['name']] = tbl_fluid1
corrtbl_afreeze[tbl_fluid2['name']] = tbl_fluid2


class HeatPump(object):

    def __init__(self):
        """
        TinHP  : temperature of the fluid entering the HP in ºC.
        Nhp    : number of heat pump
        hpname : name of the heat pump
        """

        self.__initAttr__()

        self.TinHP = {'cooling': 28, 'heating': 0}
        self.Nhp = 1
        self.set_hpname(0)

    def __initAttr__(self):
        """
        Initialize the attributes that are not to be linked with the UI

        qbat  : building thermal load in kW (+ for cooling, - for heating)
        fluid : heat carrier fluid type
        fr    : antifreeze volumetric fraction
        Vftot : total volumetric flow in the system in L/s
        Tg    : undisturbed ground temperature in ºC
        """

        self._hpdb = load_heatpump_database()

        self.qbat = {'cooling': 16.5, 'heating': 14.5}
        self.fluid = 'water'
        self.fr = 0

        self.Vftot = {}
        self.Vftot['cooling'] = 0.05 * 16.5
        self.Vftot['heating'] = 0.05 * 14.5

        self.Tg = 12

    @property
    def hpdata(self):
        """Return the performance data table of the heat pump."""
        return self._hpdb[self.hpname]

    @property
    def hpnames(self):
        """Return a list of all available heatpumps in the database."""
        return list(self._hpdb.keys())

    @property
    def hpname(self):
        """Return the name of the current heat pump."""
        return self._hpname

    def set_hpname(self, value):
        """
        Set the name of the heatpump either from an index or a key. If the
        index or key is not found in the database, an error is raised.
        """
        if isinstance(value, int):
            self._hpname = self.hpnames[0]
        elif isinstance(value, str):
            if value in list(self._hpdb.keys()):
                self._name = value
            else:
                raise ValueError("Heatpump '%s' not found in the database"
                                 % value)
        else:
            raise TypeError("'value' must be either an int or a str")

    @property
    def ToutHP(self):
        """
        Return the temperature of the water leaving the heat pump (LWT) in ºC.
        """
        ToutHP = {}
        for mode in ['cooling', 'heating']:
            ToutHP[mode] = self.calcul_ToutHP(mode)

        return ToutHP

    @property
    def Tm(self):
        """Return the fluid mean temperature through the heat pump in ºC"""
        ToutHP = self.ToutHP
        TinHP = self.TinHP

        Tm = {}
        for mode in ['cooling', 'heating']:
            Tm[mode] = (TinHP[mode] + ToutHP[mode])/2

        return Tm

    @property
    def Vhp(self):
        """Return the volumetric flowrate per HP in L/s"""
        return {'heating': self.Vftot['heating'] / self.Nhp,
                'cooling': self.Vftot['cooling'] / self.Nhp}

    @property
    def COP(self):
        """Return the coefficient of performance of the heatpump."""
        coph = self.interp('COPh', self.TinHP['heating'], self.Vhp['heating'])
        copc = self.interp('COPc', self.TinHP['cooling'], self.Vhp['cooling'])
        return {'heating': coph, 'cooling': copc}

    @property
    def CAP(self):
        """Return the capacity of the heatpump in kW."""
        caph = self.interp('CAPh', self.TinHP['heating'], self.Vhp['heating'])
        capc = self.interp('CAPc', self.TinHP['cooling'], self.Vhp['cooling'])
        return {'heating': caph, 'cooling': capc}

    # ---- Calculs

    def calcul_ToutHP(self, mode):
        """
        Calcul the temperature of the fluid leaving the heat pump for the
        current mode of operation (cooling of heating).
        """

        # Calculate fluid properties :

        hcfluid = HeatCarrierFluid(self.fluid, self.TinHP[mode], self.fr)
        rhof = hcfluid.rho
        cpf = hcfluid.cp

        # Calculate ground load :

        qbat = self.qbat[mode]
        COP = self.get_COP(mode)

        if mode == 'cooling':
            qgnd = qbat * (COP + 1) / COP
        elif mode == 'heating':
            qgnd = -qbat * (COP - 1) / COP

        # Calculate outflow fluid temperature :

        ToutHP = self.TinHP[mode] + qgnd/(self.Vftot[mode]*rhof*cpf) * 10**6

        return ToutHP

    def interp(self, varname, ewt, gpm):
        """
        Evaluate the heatpump <varname> at the specifived ewt and gpm value
        with the equation-fit model.

        varname : COPc, COPh, CAPc, CAPh
        ewt : entering water temperature in the HP (ºC)
        gpm : volumetric flowrate in the HP (L/s)
        """

        A = self._hpdb[self.hpname]['models'][varname]

        y = (A[0] +
             A[1]*ewt + A[2]*ewt**2 +
             A[3]*gpm + A[4]*gpm**2 +
             A[5]*ewt*gpm)

        # Get the anti-freeze correction factor.

        afcorr = np.interp(self.fr,
                           corrtbl_afreeze[self.fluid]['fr'],
                           corrtbl_afreeze[self.fluid][varname])

        return y * afcorr

    def in_table(self, varname, p1, p2):
        """
        Check whether the (p1, p2) point is inside or outside the data table
        for the specified variable. COP or CAP values evaluated outside the
        data table (extrapolation) must be used with care.

        Based on this stackoverflow answer:
        http://stackoverflow.com/a/16898636/4481445
        """
        y = self._hpdb[self.hpname][varname]
        indx = np.where(~np.isnan(y))[0]

        x1 = self._hpdb[self.hpname]['EWT'][indx]
        x2 = self._hpdb[self.hpname]['GPM'][indx]

        hull = Delaunay(np.vstack((x1, x2)).T)
        return bool(hull.find_simplex((p1, p2)) >= 0)

    def get_flowRange(self):
        """
        Return the minimum and maximum operational flowrate of the heatpump.
        """
        vmax = np.max(self._hpdb[self.hpname]['GPM']) * self.Nhp
        vmin = np.min(self._hpdb[self.hpname]['GPM']) * self.Nhp
        return vmin, vmax

    def get_COP(self, mode):
        """
        Return the coefficient of performance (COP) of the heat pump for the
        specified mode of operation (cooling of heating).
        """
        if mode == 'cooling':
            return self.interp('COPc', self.TinHP[mode], self.Vhp[mode])
        elif mode == 'heating':
            return self.interp('COPh', self.TinHP[mode], self.Vhp[mode])

    def get_CAP(self, mode):
        """
        Return the capacity (CAP) of the heat pump in kW for the specified
        mode of operation (cooling of heating).
        """
        if mode == 'cooling':
            return self.interp('CAPc', self.TinHP[mode], self.Vhp[mode])
        elif mode == 'heating':
            return self.interp('CAPh', self.TinHP[mode], self.Vhp[mode])

    # ---- Utility methods

    def plot_heatpump_model_goodness(self):
        """
        Produce a graph that shows the goodness of fit of the equation-fit
        models used to evaluate the COP and CAP values of the heatpump as a
        function of the entering water temperature (EWT) and volumetric
        flowrate.
        """
        plot_fitmodel_eval_from(self.hpdata)


if __name__ == '__main__':
    heatpump = HeatPump()
    heatpump.plot_heatpump_model_goodness()
    hpdb = heatpump._hpdb
    print(heatpump.hpnames)
