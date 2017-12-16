# -*- coding: utf-8 -*-

# copyright (c) 2016 Louis Lamarche
# copyright (c) 2016-2017 Jean-Sébastien
# https://github.com/jnsebgosselin/pygld
#
# This is part of PyGLD (Python Ground-Loop Designer).
# Licensed under the terms of the MIT License.

from __future__ import division, unicode_literals

import numpy as np
from scipy import interpolate
import os
import xlrd
from collections import OrderedDict


class HeatCarrierFluid(object):
    """
    Inputs:
    -------
    fluid = Name of the heat carrier fluid (water ; prop_glycol)
    Tref  = Reference temperature in °C
    fr    = Anti-freeze volumetric fraction (0 <= fr <= 1)
            If fr = 0, properties of pure water at Tref are returned.

    Primary Properties:
    -------------------
    rho = Fluid density in kg/m3
    mu  = Cinematic viscosity in Pa.s
    k   = Thermal conductivity in W/(m·k)
    cp  = Specific Heat Capacity in J/(kg·K)

    Derived Properties:
    -------------------
    nu = Kynematic viscosity in m²/s
    Pr = Prantl number
    al = Thermal diffusivity in m²/s
    Cp = Volumetric heat capacity in J/(m3.K)

    Others:
    -------
    Tfp = Freezing point temperature in °C
    Tbp = Boiling point temperature in °C

    Notes:
    ------
    Extrapolation is not allowed. A NaN value is returned if outside
    of fluid properties table.

    References:
    -----------
    Properties for pure water were taken from :

    VDI-Gesellschaft, 2010. VDI Heat Atlas. VDI-Verlag GmbH, Dusseldorf,
        Germany, Berlin; London.

    Properties for propylene and ethlyne glycol were taken from:

    ASHRAE, 2009. 2009 ASHRAE Handbook - Fundamentals. Si edition ed.,
        American Society of Heating, Refrigerating and Air-Conditioning
        Engineers, Atlanta, GA. chap.31, pp. 823-835.
    """

    def __init__(self, fluid='water', Tref=20, fr=0):

        self.fluid = fluid
        self.Tref = Tref
        self.fr = fr

    # =========================================================================

    @property                                              # Heat carrier fluid
    def fluid(self):
        return self.__fluid

    @fluid.setter
    def fluid(self, x):
        if x == 'prop_glycol':
            self.__fluid = x
            filename = 'proptables_propglycol.npy'
        elif x == 'ethyl_glycol':
            self.__fluid = x
            filename = 'proptables_ethylglycol.npy'
        else:
            self.__fluid = 'water'
            filename = 'proptables_purewater.npy'
            if x != 'water':
                print('Warning: fluid not defined, setting fluid '
                      'to pure water.')

        dirname = os.path.dirname(os.path.realpath(__file__))
        pathname = os.path.join(dirname, 'tables', filename)

        # TTP: Table of Thermophysical Properties
        self.__TTP = np.load(pathname)

    # =========================================================================

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
            print('0 > fr > 1')
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

        # interpolate value from table for aqueous solutions and pure water
        # x is temperature
        # y is the antifreeze volumetric fraction
        # z is thermodynamic variable

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


# =============================================================================


def purewaterprop2binary(rootname):
    """
    Load the property table of pure water from the Excel file, format
    the data and save the results in a binary file.
    """

    # ------------------------------- Extract data from the Excel Workbook ----

    wb = xlrd.open_workbook(rootname+'.xls')

    data = []
    sheet = wb.sheet_by_index(0)
    for row in range(3, sheet.nrows):
        data.append(sheet.row_values(row, start_colx=0, end_colx=None))

    wb.release_resources()

    # ------------------------------------- Format Data and Save to Binary ----

    data = np.array(data).astype(float)
    N, M = np.shape(data)

    x = data[:, 0]
    y = np.zeros(N)

    A = OrderedDict([('density', []), ('spec_heat', []),
                     ('ther_cond', []), ('viscosity', [])])

    # OrderedDict is needed since the order of the keys is important to
    # respect the structured of the Excel table.

    for i, key in enumerate(A.keys()):
        z = data[:, i+1]

        # Removes NaNs and store the results in a [x, y, z] matrix:

        indx = np.where(~np.isnan(z))[0]
        A[key] = np.array([x[indx], y[indx], z[indx]])

    A['freez_point'] = 0
    A['boil_point'] = 100

    # -------------------------------------------------- Produce a ndarray ----

    dtype = []
    for key in A.keys():
        dtype.append((key, 'float64', np.shape(A[key])))

    A2 = np.zeros((), dtype=dtype)
    for key in A.keys():
        A2[key] = A[key]

    np.save(rootname+'.npy', A2)
    wb.release_resources()

    np.save(rootname+'.npy', A2)


# =============================================================================


def brineprop2binary(rootname):
    """
    Load the property table of the brine from the Excel file, format
    the data and save the results in a binary file.
    """

    wb = xlrd.open_workbook(rootname+'.xls')

    # ------------------ density, spec. heat, thermal cond., and viscosity ----

    A = {'density': [], 'spec_heat': [], 'ther_cond': [], 'viscosity': []}
    for key in A.keys():

        # Retrieve data from sheet:

        sheet = wb.sheet_by_name(key)
        data = []
        for row in range(2, sheet.nrows):
            data.append(sheet.row_values(row, start_colx=0, end_colx=None))

        # Format data:

        x = np.array(data[1:])[:, 0]   # fluid temperature
        x = x.astype(float)

        y = np.array(data[0])[1:]      # anti-freeze volumetric fraction
        y = y.astype(float)

        z = np.array(data[1:])[:, 1:]  # thermophysical variable
        z = z.astype(float)

        # Flatten the results:

        N = len(x)
        M = len(y)

        x = np.repeat(x, M)
        y = np.tile(y, N)
        z = z.ravel()

        # Removes NaNs and store the results in a [x, y, z] matrix:

        indx = np.where(~np.isnan(z))[0]
        A[key] = np.array([x[indx], y[indx], z[indx]])

    # ----------------------------------------- freezing and boiling point ----

    sheet = wb.sheet_by_name('freezing and boiling point')
    data = []
    for row in range(3, sheet.nrows):
        data.append(sheet.row_values(row, start_colx=0, end_colx=None))
    data = np.array(data).astype(float)

    # Format data:

    fr = data[:, 1]/100
    Tfp = data[:, 2]
    Tbp = data[:, 3]

    # Removes NaNs and produces a [y, z] matrix:

    indx = np.where(~np.isnan(Tfp))[0]
    A['freez_point'] = np.array([fr[indx], Tfp[indx]])

    indx = np.where(~np.isnan(Tbp))[0]
    A['boil_point'] = np.array([fr[indx], Tbp[indx]])

    # -------------------------------------------------- Produce a ndarray ----

    dtype = []
    for key in A.keys():
        dtype.append((key, 'float64', np.shape(A[key])))

    A2 = np.zeros((), dtype=dtype)
    for key in A.keys():
        A2[key] = A[key]

    np.save(rootname+'.npy', A2)
    wb.release_resources()


if __name__ == '__main__':
    brineprop2binary('tables/proptables_propglycol')
    brineprop2binary('tables/proptables_ethylglycol')
    purewaterprop2binary('tables/proptables_purewater')
