# -*- coding: utf-8 -*-

# Copyright © 2018 Jean-Sébastien Gosselin
# https://github.com/jnsebgosselin/qwatson
#
# This file is part of PyGLD.
# Licensed under the terms of the GNU General Public License.

# ---- Standard imports

import copy

# ---- Third party imports

import numpy as np
from scipy.spatial import Delaunay

# ---- Local imports

from pygld.fluidproperties import HeatCarrierFluid
from pygld.heatpumps.utils import load_heatpump_database
from pygld.heatpumps.maths import eval_polyfid2rd
from pygld.heatpumps.plots import plot_fitmodel_eval_from

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

corrtbl_afreeze = {}
corrtbl_afreeze[tbl_fluid0['name']] = tbl_fluid0
corrtbl_afreeze[tbl_fluid1['name']] = tbl_fluid1
corrtbl_afreeze[tbl_fluid2['name']] = tbl_fluid2


class HeatPump(object):

    def __init__(self):
        self._hpdb = load_heatpump_database()
        self.set_model(0)

        # Independent properties :

        self._TinHP = np.array([28, 0])
        self._qbat = np.array([-16.5, 14.5])
        self._Vf = np.array([0.94635295, 0.94635295])
        self._fr = 0
        self._hcfluid = HeatCarrierFluid('water', self._TinHP, self._fr)

        # Setup the update dependent properties flag :

        self._dependent_props = ['Tm', 'ToutHP', 'COP', 'CAP']
        self._varstate_has_changed()

    def __str__(self):
        str_ = "model: %s" % self.model
        str_ += '\nfluid: %s' % self.fluid
        str_ += '\nfr   : %d' % self.fr
        str_ += '\nVf (L/s): ' + self.Vf.__str__()
        str_ += "\nqbat (kW): " + self.qbat.__str__()
        str_ += '\nTinHP (\u00B0C): ' + self.TinHP.__str__()
        str_ += '\nToutHP (\u00B0C): ' + self.ToutHP.__str__()
        str_ += '\nTm (\u00B0C): ' + self.Tm.__str__()
        str_ += '\nCOP : ' + self.COP.__str__()
        str_ += '\nCAP (kW): ' + self.CAP.__str__()
        return str_

    # ---- Model and data

    @property
    def hpdata(self):
        """Return the performance data table of the heatpump."""
        return self._hpdb[self.model]

    @property
    def model(self):
        """Return the name of the current heatpump."""
        return self._model

    def set_model(self, value):
        """
        Set the name of the heatpump either from an index or a key. If the
        index or key is not found in the database, an error is raised.
        """
        if isinstance(value, int):
            self._model = self.get_avail_heatpump_models()[value]
        elif isinstance(value, str):
            if value in list(self._hpdb.keys()):
                self._name = value
            else:
                raise ValueError("Heatpump '%s' not found in the database"
                                 % value)
        else:
            raise TypeError("'value' must be either an int or a str")

    # ---- Independent properties

    @property
    def fluid(self):
        """Type of heat carrier fluid."""
        return copy.copy(self._hcfluid.fluid)

    @fluid.setter
    def fluid(self, value):
        self._hcfluid.fluid = value
        self._varstate_has_changed()

    @property
    def fr(self):
        """Volumetric fraction of antifreeze of the heat carrier fluid."""
        return copy.copy(self._hcfluid.fr)

    @fr.setter
    def fr(self, value):
        self._hcfluid.fr = value
        self._varstate_has_changed()

    @property
    def TinHP(self):
        """
        Temperature of the water entering the heatpump (EWT) in ºC.
        """
        return np.copy(self._TinHP)

    @TinHP.setter
    def TinHP(self, values):
        values = values if isinstance(values, np.ndarray) else np.array(values)
        self._TinHP = values
        self._hcfluid.Tref = values
        self._varstate_has_changed()

    @property
    def qbat(self):
        """
        Part of the building thermal load applied to the heatpump in kW
        (- for cooling, + for heating).
        """
        return np.copy(self._qbat)

    @qbat.setter
    def qbat(self, values):
        values = values if isinstance(values, np.ndarray) else np.array(values)
        self._qbat = values
        self._varstate_has_changed()

    @property
    def Vf(self):
        """
        Volumetric flowrate of the fluid carrier through the heatpump in L/s
        """
        return np.copy(self._Vf)

    @Vf.setter
    def Vf(self, values):
        values = values if isinstance(values, np.ndarray) else np.array(values)
        self._Vf = values
        self._varstate_has_changed()

    def _varstate_has_changed(self):
        """
        Reset the '_need_update' flags for all the dependent variables
        when the value of an independent variable changes.
        """
        self._need_update = {p: True for p in self._dependent_props}

    # ---- Dependent properties

    @property
    def ToutHP(self):
        """
        Temperature of the water leaving the heatpump (LWT) in ºC.
        """
        if self._need_update['ToutHP']:
            self._ToutHP = self.calcul_ToutHP()
            self._need_update['ToutHP'] = False
        return np.copy(self._ToutHP)

    @property
    def Tm(self):
        """
        Mean temperature of the water circulating through the heatpump in ºC
        """
        if self._need_update['Tm']:
            self._Tm = self.calcul_Tm()
            self._need_update['Tm'] = False
        return np.copy(self._Tm)

    @property
    def COP(self):
        """Coefficient of performance of the heatpump."""
        if self._need_update['COP']:
            self._COP = self.calcul_COP()
            self._need_update['COP'] = False
        return np.copy(self._COP)

    @property
    def CAP(self):
        """Capacity of the heatpump in kW."""
        if self._need_update['CAP']:
            self._CAP = self.calcul_CAP()
            self._need_update['CAP'] = False
        return np.copy(self._CAP)

    # ---- Calculs

    def calcul_Tm(self):
        """
        Calcul the average temperature of the fluid circulating through the
        heatpump in ºC for every TinHP values.
        """
        return (self.TinHP + self.ToutHP)/2

    def calcul_ToutHP(self):
        """
        Calcul the temperature of the fluid leaving the heatpump in ºC for
        every TinHP values.
        """
        TinHP, Vf, qbat = self.TinHP, self.Vf, self.qbat
        COP = self.COP

        if len(TinHP) != len(Vf) != len(qbat):
            raise ValueError("The lenght of TinHP, Vf, and qbat must match"
                             " exactly.")

        # Get the fluid properties.

        rhof, cpf = self._hcfluid.rho, self._hcfluid.cp

        # Calcul the ground loads.

        qgnd = qbat * (COP - np.sign(qbat)) / COP

        # qbat is - for cooling and + for heating
        # qgnd = qbat * (COP + 1)/COP in cooling mode
        # qgnd = qbat * (COP - 1)/COP in heating mode

        # Calculate outflow fluid temperature.

        ToutHP = TinHP - qgnd/(Vf*rhof*cpf) * 10**6

        return ToutHP

    def calcul_COP(self):
        """
        Calcul the coefficient of performance of the heatpump for each pair
        of TinHP and Vf values.
        """
        return self._calcul_cop_or_cap('COP')

    def calcul_CAP(self):
        """
        Calcul the capacity of the heatpump in kW for each pair of TinHP and
        Vf values.
        """
        return self._calcul_cop_or_cap('CAP')

    def _calcul_cop_or_cap(self, varname):
        """Calcul the CAP or COP for each pair of TinHP and Vf values."""
        TinHP, Vf, qbat = self.TinHP, self.Vf, self.qbat

        if len(TinHP) != len(Vf) != len(qbat):
            raise ValueError("The lenght of TinHP, Vf, and qbat must match"
                             " exactly.")

        ndarray = np.empty(len(TinHP))

        # Calcul values when cooling.

        indx = np.where(qbat < 0)[0]
        ndarray[indx] = self.eval_fitmodel_for(
            varname + 'c', TinHP[indx], Vf[indx])

        # Calcul values when heating.

        indx = np.where(qbat >= 0)[0]
        ndarray[indx] = self.eval_fitmodel_for(
            varname + 'h', TinHP[indx], Vf[indx])

        return ndarray

    def eval_fitmodel_for(self, varname, ewt, vf):
        """
        Evaluate the heatpump COP or CAP at the specified ewt and vf value
        with the equation-fit model.

        varname : COPc, COPh, CAPc, CAPh
        ewt : entering water temperature in the HP (ºC)
        vf : volumetric flowrate in the HP (L/s)
        """
        # Get the equation-fit model coefficients and evaluate the values of
        # varname for the provided independent variables pair of points.

        A = self._hpdb[self.model]['eqfit_models'][varname]
        y = eval_polyfid2rd(A, ewt, vf)

        # Get the anti-freeze correction factor.

        afcorr = np.interp(self.fr,
                           corrtbl_afreeze[self.fluid]['fr'],
                           corrtbl_afreeze[self.fluid][varname])

        return y * afcorr

    # ---- Utility methods

    def in_table(self, varname, p1, p2):
        """
        Check whether the (p1, p2) point is inside or outside the data table
        for the specified variable. COP or CAP values evaluated outside the
        data table (extrapolation) must be used with care.

        Based on this stackoverflow answer:
        http://stackoverflow.com/a/16898636/4481445
        """
        y = self._hpdb[self.model][varname]
        indx = np.where(~np.isnan(y))[0]

        x1 = self._hpdb[self.model]['EWT'][indx]
        x2 = self._hpdb[self.model]['GPM'][indx]

        hull = Delaunay(np.vstack((x1, x2)).T)
        return bool(hull.find_simplex((p1, p2)) >= 0)

    def get_flowRange(self):
        """
        Return the minimum and maximum operational flowrate of the heatpump.
        """
        vmax = np.max(self._hpdb[self.model]['GPM'])
        vmin = np.min(self._hpdb[self.model]['GPM'])
        return vmin, vmax

    def plot_heatpump_model_goodness(self):
        """
        Produce a graph that shows the goodness of fit of the equation-fit
        models used to evaluate the COP and CAP values of the heatpump as a
        function of the entering water temperature (EWT) and volumetric
        flowrate.
        """
        plot_fitmodel_eval_from(self.hpdata)

    def get_avail_heatpump_models(self):
        """Return a list of all available heatpump models in the database."""
        return list(self._hpdb.keys())

    def print_avail_heatpump_models(self):
        """Print the list of the available heatpump models in the database."""
        models = self.get_avail_heatpump_models()
        N = len(models)
        max_indent = (N-1)//10
        for i, model in enumerate(models):
            indent = ' ' * (max_indent - i//100)
            print("%s%d - %s" % (indent, i, model))


if __name__ == '__main__':
    heatpump = HeatPump()
    print(heatpump)
    # print(heatpump)
    # heatpump.plot_heatpump_model_goodness()
    # heatpump.print_avail_heatpump_models()
