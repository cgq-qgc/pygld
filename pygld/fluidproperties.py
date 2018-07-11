# -*- coding: utf-8 -*-

# Copyright © 2017-2018 Jean-Sébastien Gosselin
# https://github.com/jnsebgosselin/pygld
#
# This file is part of PyGLD.
# Licensed under the terms of the MIT License.

# ---- Standard imports

import copy

# ---- Third party imports

import numpy as np
from scipy import interpolate
import os

# ---- Local imports

from pygld.utils.strformating import array_to_str

FLUIDS = ['prop_glycol', 'ethyl_glycol', 'water']


class HeatCarrierFluid(object):
    """
    The :attr:`~pygld.HeatCarrierFluid` holds all the thermophysical
    properties related to the fluid circulating in the ground-loop of the
    heat exchanger system. The :attr:`~pygld.HeatCarrierFluid` is initialized
    by default as a pure water :attr:`~pygld.HeatCarrierFluid.fluid` with an
    antifreeze volumetric fraction (:attr:`~pygld.HeatCarrierFluid.fr`) of 0
    and at a reference temperature (:attr:`~pygld.HeatCarrierFluid.Tref`)
    of 20°C.

    The fluid type of :attr:`~pygld.HeatCarrierFluid` can be changed with the
    :meth:`~pygld.HeatCarrierFluid.set_fluid` method and a list of all
    available heat carrier fluid types can be obtained from the
    :meth:`~pygld.HeatCarrierFluid.get_avail_fluid_types` method. The
    reference temperature and antifreeze volumetric fraction can be changed
    by setting directly the value of the :attr:`~pygld.HeatCarrierFluid.Tref`
    and :attr:`~pygld.HeatCarrierFluid.fr` attributes. Printing an instance of
    :attr:`~pygld.HeatCarrierFluid` will print a summary of the fluid's
    independent and dependent properties.

    The freezing and boiling points of the
    :attr:`~pygld.HeatCarrierFluid.fluid` are evaluated for a given value of
    :attr:`~pygld.HeatCarrierFluid.fr` from a piecewise one-dimensional cubic
    interpolation of the fluid's thermophysical properties table.
    The thermophysical properties of the
    :attr:`~pygld.HeatCarrierFluid.fluid` are determined for a given set of
    :attr:`~pygld.HeatCarrierFluid.Tref` and :attr:`~pygld.HeatCarrierFluid.fr`
    values from a two dimensional piecewise cubic interpolation of the
    fluid's thermophysical properties table.
    The derived thermophysical properties of the
    :attr:`~pygld.HeatCarrierFluid.fluid` are calculated from the interpolated
    values of the primary properties.

    An `Example`_ is available at the end of this section.

    .. note::
        Extrapolation is not allowed when evaluating the
        thermophysical properties of the fluid.
        A `nan` value is returned for sets of
        :attr:`~pygld.HeatCarrierFluid.Tref` and
        :attr:`~pygld.HeatCarrierFluid.fr` whose values fall outside of the
        thermophysical properties table of the fluid.
    """

    def __init__(self, fluid='water', Tref=20, fr=0):
        self.set_fluid(fluid)
        self.Tref = Tref
        self.fr = fr

    def __str__(self):
        str_ = "Type of fluid: %s" % self.fluid
        str_ += '\nAntifreeze volumetric fraction: %0.2f' % self.fr
        str_ += '\nFreezing point temperature (°C): %0.1f' % self.Tfp
        str_ += '\nBoiling point temperature (°C): %0.1f' % self.Tbp

        str_ += '\nTemperature of reference (°C): '
        str_ += array_to_str(self.Tref, "{:.1f}")
        str_ += '\nFluid density in (kg/m³): '
        str_ += array_to_str(self.rho, "{:.2f}")
        str_ += '\nCinematic viscosity (Pa·s): '
        str_ += array_to_str(self.mu, "{:.2e}")
        str_ += '\nThermal conductivity (W/m·k): '
        str_ += array_to_str(self.kth, "{:.3f}")
        str_ += '\nSpecific heat capacity (J/kg·K): '
        str_ += array_to_str(self.cp, "{:.1f}")

        str_ += '\nKynematic viscosity (m²/s): '
        str_ += array_to_str(self.nu, "{:.2e}")
        str_ += '\nPrantl number: '
        str_ += array_to_str(self.Pr, "{:.1f}")
        str_ += '\nThermal diffusivity (m²/s): '
        str_ += array_to_str(self.al, "{:.2e}")
        str_ += '\nVolumetric Heat Capacity (J/m³·K): '
        str_ += array_to_str(self.Cp, "{:.2e}")
        return str_

    @property                                              # Heat carrier fluid
    # ---- Fluid type and data

    @property
    def hcfdata(self):
        """
        Return a `dict` containing the tables of thermophysical properties
        of the heat carrier :attr:`~pygld.HeatCarrierFluid.fluid`.
        """
        return copy.copy(self.__TTP)

    def fluid(self):
        """Return the type of fluid of the heat carrier.

        The type of fluid of the heat carrier is set to pure water by default.
        The heat carrier fluid is assumed to be pure water when
        :attr:`~pygld.HeatCarrierFluid.fr` is set to 0.
        """
        return self._fluid

    def set_fluid(self, x):
        """Set the type of fluid of the heat carrier fluid and
        load the tables of thermophysical properties of the fluid.

        If the specified fluid type is not available, an error is raised.
        The list of available heat carrier fluid types can be obtained with
        the :meth:`~pygld.HeatCarrierFluid.get_avail_fluid_types` method.
        """
        if x == 'prop_glycol':
            self._fluid = x
            filename = 'proptables_propglycol.npy'
        elif x == 'ethyl_glycol':
            self._fluid = x
            filename = 'proptables_ethylglycol.npy'
        elif x == 'water':
            self._fluid = 'water'
            filename = 'proptables_purewater.npy'
        else:
            raise ValueError('Supported fluid value are', FLUIDS)

        dirname = os.path.dirname(os.path.realpath(__file__))
        pathname = os.path.join(dirname, 'tables', filename)

        # TTP: Table of Thermophysical Properties
        self.__TTP = np.load(pathname)


    @property                                               # Temperature in °C
    def Tref(self):
        return self.__Tref

    @Tref.setter
    def Tref(self, x):
        self.__Tref = x

    # =========================================================================

    @property                                  # antifreeze volumetric fraction
    def fr(self):                                              # (0 <= fr <= 1)
        if self.fluid == 'water':
            return 0
        else:
            return self.__fr

    @fr.setter
    def fr(self, x):
        if x > 1 or x < 0:
            raise ValueError('fr must be between 0 and 1 (0 <= fr < 1).')
        else:
            self.__fr = x

    # =========================================================================

    @property                                # Freezing point temperature in °C
    def Tfp(self):
        if self.fluid == 'water':
            return self.__TTP['freez_point']
        else:
            x = self.fr
            xp = self.__TTP['freez_point'][0]
            yp = self.__TTP['freez_point'][1]

            return np.interp(x, xp, yp)

    # =========================================================================

    @property                                 # Boiling point temperature in °C
    def Tbp(self):
        if self.fluid == 'water':
            return self.__TTP['boil_point']
        else:
            x = self.fr
            xp = self.__TTP['boil_point'][0]
            yp = self.__TTP['boil_point'][1]

            return np.interp(x, xp, yp)

    # =========================================================================

    @property                                          # Fluid density in kg/m3
    def rho(self):
        x = self.__TTP['density'][0]
        y = self.__TTP['density'][1]
        z = self.__TTP['density'][2]

        return self.interp(x, z, y)

    # =========================================================================

    @property                                     # Cinematic viscosity in Pa.s
    def mu(self):
        x = self.__TTP['viscosity'][0]
        y = self.__TTP['viscosity'][1]
        z = self.__TTP['viscosity'][2]/1000

        return self.interp(x, z, y)

    # =========================================================================

    @property                                 # Thermal conductivity in W/(m·k)
    def k(self):
        x = self.__TTP['ther_cond'][0]
        y = self.__TTP['ther_cond'][1]
        z = self.__TTP['ther_cond'][2]

        return self.interp(x, z, y)

    # =========================================================================

    @property                              # Specific Heat Capacity in J/(kg·K)
    def cp(self):
        x = self.__TTP['spec_heat'][0]
        y = self.__TTP['spec_heat'][1]
        z = self.__TTP['spec_heat'][2]*1000

        return self.interp(x, z, y)

    # =========================================================================

    @property                                     # Kynematic viscosity in m²/s
    def nu(self):
        return self.mu/self.rho

    # =========================================================================

    @property                                                   # Prantl number
    def Pr(self):
        return self.mu * self.cp / self.k

    # =========================================================================

    @property                                     # Thermal diffusivity in m²/s
    def al(self):
        return self.k / (self.cp * self.rho)

    # =========================================================================

    @property                            # Volumetric Heat Capacity in J/(m3.K)
    def Cp(self):
        return self.cp * self.rho

    # =========================================================================

    def interp(self, x, z, y=None):
        """
        Interpolate a value from table for aqueous solutions and pure water
        x is temperature
        y is the antifreeze volumetric fraction
        z is thermodynamic variable
        """
        if self.fluid == 'water':
            x1 = self.Tref
            z1 = interpolate.griddata(x, z, x1, method='cubic')
        else:
            x1 = self.Tref
            if np.size(x1) == 1:
                y1 = self.fr
            else:
                y1 = np.ones(np.size(x1)) * self.fr
            z1 = interpolate.griddata((x, y), z, (x1, y1), method='cubic')

        return z1


if __name__ == '__main__':
    hcf = HeatCarrierFluid('prop_glycol')
    print(hcf.rho)
