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

from pygld.fluidproperties import HeatCarrierFluid, FLUIDS
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
    """
    The :attr:`~pygld.HeatPump` holds all the properties and handles all the
    calculation and modelling related to a single heatpump coupled to a
    ground-loop heat exchanger system.

    An example is available in this
    `notebook <https://github.com/jnsebgosselin/pygld/blob/master/examples/
    example_heatpump.ipynb>`_.
    """

    def __init__(self):
        # Setup the update dependent properties flag :

        self._dependent_props = ['Tm', 'ToutHP', 'COP', 'CAP']
        self._varstate_has_changed()

        # Setup the independent properties :

        self._TinHP = np.array([28, 0])
        self._qbat = np.array([-16.5, 14.5])
        self._Vf = np.array([0.94635295, 0.94635295])
        self._fr = 0
        self._hcfluid = HeatCarrierFluid('water', self._TinHP, self._fr)

        # Load the database and set the model.

        self._hpdb = load_heatpump_database()
        self.set_model(0)

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
        """Return the performance data of the heatpump stored in a dict."""
        return self._hpdb[self.model]

    @property
    def model(self):
        """Return the fabrication model of the heatpump."""
        return self._model

    def set_model(self, value):
        """Set the fabrication model of the heatpump.

        Set the fabrication model of the heatpump, either from an `index` or a
        `str`. If the `index` or the `str` is not a valid key of the heatpump
        models database, an error is raised.
        The list of heatpump models available in the databse can be obtained
        with the :meth:`~pygld.HeatPump.get_avail_heatpump_models` method.
        """
        if isinstance(value, int):
            self._model = self.get_avail_heatpump_models()[value]
            self._varstate_has_changed()
        elif isinstance(value, str):
            if value in list(self._hpdb.keys()):
                self._model = value
                self._varstate_has_changed()
            else:
                raise ValueError("Heatpump '%s' not found in the database"
                                 % value)
        else:
            raise TypeError("'value' must be either an int or a str")

    # ---- Independent properties

    @property
    def fluid(self):
        """Heat carrier fluid type.

        Get or set the type of the heat carrier fluid used in the geothermal
        heat exchanger. By default, the heat carrier is set to pure water.
        The list of available heat carrier fluid types can be obtained with
        the :meth:`~pygld.HeatPump.get_avail_fluid_types` method.
        The heat carrier fluid is assumed to be pure water when
        :attr:`~pygld.HeatPump.fr` is set to 0.
        """
        return self._hcfluid.fluid

    @fluid.setter
    def fluid(self, x):
        self._hcfluid.fluid = x
        self._varstate_has_changed()

    @property
    def fr(self):
        """Antifreeze volumetric fraction of the heat carrier fluid in m³/m³
        (0 ≤ fr ≤ 1).

        Get or set the antifreeze volumetric fraction of the heat carrier
        fluid. The value of `fr` must be between 0 and 1 and is assumed to be
        0 when :attr:`~pygld.HeatPump.fluid` is set to 'water'.
        """
        return self._hcfluid.fr

    @fr.setter
    def fr(self, x):
        self._hcfluid.fr = x
        self._varstate_has_changed()

    @property
    def TinHP(self):
        """Temperature of the water entering the heatpump in °C.

        Get or set the temperature of the water entering the heatpump as a
        single value or a series of values. A :attr:`~pygld.HeatPump.ToutHP`,
        :attr:`~pygld.HeatPump.Tm`, :attr:`~pygld.HeatPump.COP`, and
        :attr:`~pygld.HeatPump.CAP` value is calculated for every value of
        :attr:`~pygld.HeatPump.TinHP`. A numpy array will always be returned
        when getting :attr:`~pygld.HeatPump.TinHP`, independently of the
        format used to set the attribute.
        """
        return np.copy(self._TinHP)

    @TinHP.setter
    def TinHP(self, x):
        x = np.array([x]) if not hasattr(x, '__iter__') else np.array(x)
        self._TinHP = x.astype(float)
        self._hcfluid.Tref = x
        self._varstate_has_changed()

    @property
    def qbat(self):
        """Building thermal load applied to the heatpump in kW
        (- for cooling, + for heating).

        Get or set the portion of the building thermal load that is applied to
        the heatpump as a single value or a series of values. The lenght of
        the data series used to set the attribute must match exactly that of
        the :attr:`~pygld.HeatPump.TinHP`, or an error will be raised when
        computing the dependent properties (:attr:`~pygld.HeatPump.ToutHP`,
        :attr:`~pygld.HeatPump.Tm`, :attr:`~pygld.HeatPump.COP`, and
        :attr:`~pygld.HeatPump.CAP`). Negative values are used to
        represent the cooling thermal loads of the building while positive
        values are used to represent the heating thermal loads.
        """
        return np.copy(self._qbat)

    @qbat.setter
    def qbat(self, x):
        x = np.array([x]) if not hasattr(x, '__iter__') else np.array(x)
        self._qbat = x.astype(float)
        self._varstate_has_changed()

    @property
    def Vf(self):
        """Volumetric flowrate of the heat carrier fluid in L/s.

        Get or set the volumetric flowrate of the heat carrier fluid
        circulating through the heatpump as a single value or a series of
        values. The lenght of the data series used to set the attribute must
        match exactly that of the :attr:`~pygld.HeatPump.TinHP`, or an error
        will be raised when computing the dependent properties
        (:attr:`~pygld.HeatPump.ToutHP`, :attr:`~pygld.HeatPump.Tm`,
        :attr:`~pygld.HeatPump.COP`, and :attr:`~pygld.HeatPump.CAP`).
        """
        return np.copy(self._Vf)

    @Vf.setter
    def Vf(self, x):
        x = np.array([x]) if not hasattr(x, '__iter__') else np.array(x)
        self._Vf = x.astype(float)
        self._varstate_has_changed()

    # ---- Dependent properties

    @property
    def ToutHP(self):
        """Temperature of the water leaving the heatpump in °C.

        Get the temperature of the water leaving the heatpump as a series of
        values stored in a numpy array of a length that match that of
        :attr:`~pygld.HeatPump.TinHP`.
        """
        return np.copy(self._calcul_ToutHP())

    @property
    def Tm(self):
        """Mean temperature of the water circulating through the
        heatpump in °C.

        Get the mean temperature of the water circulating through the heatpump
        as series of values stored in a numpy array of a length that match that
        of :attr:`~pygld.HeatPump.TinHP`.
        """
        return np.copy(self._calcul_Tm())

    @property
    def COP(self):
        """Coefficient of performance of the heatpump.

        Get the coefficient of performance of the heatpump as series of values
        stored in a numpy array of a length that match that of
        :attr:`~pygld.HeatPump.TinHP`. The coefficients are calculated either
        for the cooling or heating mode according to the sign of the values
        set for :attr:`~pygld.HeatPump.qbat`.
        """
        return np.copy(self._calcul_COP())

    @property
    def CAP(self):
        """Capacity of the heatpump in kW.

        Get the capacity of the heatpump as series of values stored in a numpy
        array of a length that match that of :attr:`~pygld.HeatPump.TinHP`.
        The capacities are calculated either for the cooling or heating mode
        according to the sign of the values set for
        :attr:`~pygld.HeatPump.qbat`.
        """
        return np.copy(self._calcul_CAP())

    # ---- Calculs

    def _varstate_has_changed(self):
        """
        Reset the '_need_update' flags for all the dependent variables
        when the value of an independent variable changes.
        """
        self._need_update = {p: True for p in self._dependent_props}

    def _calcul_Tm(self):
        """
        Calcul the average temperature of the fluid circulating through the
        heatpump in ºC for every TinHP values.
        """
        if self._need_update['Tm']:
            ToutHP = self._calcul_ToutHP()
            self._Tm = (self._TinHP + ToutHP)/2
            self._need_update['Tm'] = False
        return self._Tm

    def _calcul_ToutHP(self):
        """
        Calcul the temperature of the fluid leaving the heatpump in ºC for
        every TinHP values.
        """
        if self._need_update['ToutHP']:
            # Get the fluid properties.

            rhof, cpf = self._hcfluid.rho, self._hcfluid.cp

            # Calcul the ground loads.

            self._calcul_COP()
            qgnd = self._qbat * (self._COP - np.sign(self._qbat)) / self._COP

            # qbat is - for cooling and + for heating
            # qgnd = qbat * (COP + 1)/COP in cooling mode
            # qgnd = qbat * (COP - 1)/COP in heating mode

            # Calculate outflow fluid temperature.

            self._ToutHP = self._TinHP - qgnd/(self._Vf*rhof*cpf) * 10**6
            self._need_update['ToutHP'] = False
        return self._ToutHP

    def _calcul_COP(self):
        """
        Calcul the coefficient of performance of the heatpump for each pair
        of TinHP and Vf values.
        """
        if self._need_update['COP']:
            self._COP = self._calcul_cop_or_cap('COP')
            self._need_update['COP'] = False
        return self._COP

    def _calcul_CAP(self):
        """
        Calcul the capacity of the heatpump in kW for each pair of TinHP and
        Vf values.
        """
        if self._need_update['CAP']:
            self._CAP = self._calcul_cop_or_cap('CAP')
            self._need_update['CAP'] = False
        return self._CAP

    def _calcul_cop_or_cap(self, varname):
        """Calcul the CAP or COP for each pair of TinHP and Vf values."""
        if not (len(self._TinHP) == len(self._Vf) == len(self._qbat)):
            raise ValueError("The lenght of TinHP, Vf, and qbat must match"
                             " exactly.")

        ndarray = np.empty(len(self._TinHP))

        # Calcul values when cooling.

        indx = np.where(self._qbat < 0)[0]
        ndarray[indx] = self._eval_fitmodel_for(
            varname + 'c', self._TinHP[indx], self._Vf[indx])

        # Calcul values when heating.

        indx = np.where(self._qbat >= 0)[0]
        ndarray[indx] = self._eval_fitmodel_for(
            varname + 'h', self._TinHP[indx], self._Vf[indx])

        return ndarray

    def _eval_fitmodel_for(self, varname, ewt, vf):
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

    def in_table(self, varname, TinHP, Vf):
        """
        Return whether the (TinHP, Vf) pair of values is inside or outside the
        performance data table of the heatpump for the specified varname.
        COP or CAP values evaluated outside the data table (extrapolation)
        must be used with care.
        """
        # Based on this stackoverflow answer:
        # http://stackoverflow.com/a/16898636/4481445

        y = self._hpdb[self.model][varname]
        indx = np.where(~np.isnan(y))[0]

        x1 = self._hpdb[self.model]['EWT'][indx]
        x2 = self._hpdb[self.model]['GPM'][indx]

        hull = Delaunay(np.vstack((x1, x2)).T)
        return bool(hull.find_simplex((TinHP, Vf)) >= 0)

    def get_flowrange(self):
        """
        Return the minimum and maximum operational flowrate of the heatpump
        in L/s.
        """
        vmax = np.max(self._hpdb[self.model]['GPM'])
        vmin = np.min(self._hpdb[self.model]['GPM'])
        return vmin, vmax

    def plot_heatpump_model_goodness(self):
        """
        Produce a graph that shows the goodness of fit of the equation-fit
        models used to evaluate the COP and CAP values of the heatpump as a
        function of the entering water temperature (TinHP) and the volumetric
        flowrate (Vf).
        """
        plot_fitmodel_eval_from(self.hpdata)

    def get_avail_fluid_types(self):
        """Return a list of all available heat carrier fluid types."""
        return copy.copy(FLUIDS)

    def get_avail_heatpump_models(self):
        """Return a list of all available heatpump models in the database."""
        return list(self._hpdb.keys())

    def print_avail_heatpump_models(self):
        """Print the list of all available heatpump models in the database."""
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
