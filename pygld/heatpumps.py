# -*- coding: utf-8 -*-
"""
copyright (C) 2016 INRS
contact: jean-sebastien.gosselin@outlook.com

This is part of OpenGLD

OpenGLD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it /will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import division, unicode_literals

from PySide import QtGui, QtCore
from PySide.QtGui import QLabel
import os
import numpy as np
from scipy import interpolate
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
import csv
from collections import OrderedDict
from geothermie.properties_mod import HCFluid

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


def build_database():
    dirname = os.path.dirname(os.path.realpath(__file__))
    dirname = os.path.join(dirname, 'GHP')

    db = OrderedDict()
    for f in os.listdir(dirname):
        if f.endswith('.csv'):
            try:
                dat = load_HP_table_fromfile(os.path.join(dirname, f))
                name = dat['model']
                db[name] = dat
            except Exception:
                print('unable to load data from %s' % f)

    np.save(os.path.join(dirname, 'hp_database.npy'), db)


def load_HP_table_fromfile(filename):
    with open(filename, 'r', encoding='utf8') as f:
        reader = list(csv.reader(f))

        dat = {'filename': filename}
        for i, row in enumerate(reader):
            if row[0] == 'Manufacturer':
                dat['manufacturer'] = row[1]
            elif row[0] == 'Model':
                dat['model'] = row[1]
                print('Loading %s data sheet...' % row[1])
            elif row[0] == 'Nominal CAP (tons)':
                dat['Nominal CAP'] = float(row[1])*3.51685
            elif row[0] == 'EWT (F)':
                A = np.array(reader[i+1:]).astype(float)

                dat['EWT'] = (A[:, 0]-32)/1.8       # temperature in ºC
                dat['GPM'] = A[:, 1]/15.8503230745  # volumetric flow in L/s

                dat['WPDc'] = A[:, 2]*6894.76   # pressure drop in coil in Pa
                dat['WPDh'] = A[:, 5]*6894.76

                dat['CAPc'] = A[:, 3]*0.293071  # capacity in kW
                dat['CAPh'] = A[:, 6]*0.293071

                dat['Wc'] = A[:, 4]  # hp work in kW
                dat['Wh'] = A[:, 7]

                break

    dat['models'] = {}

    x1 = dat['EWT']
    x2 = dat['GPM']

    dat['models']['CAPc'] = linalg_hp(dat['CAPc'], x1, x2)
    dat['models']['CAPh'] = linalg_hp(dat['CAPh'], x1, x2)

    dat['models']['COPc'] = linalg_hp(dat['CAPc']/dat['Wc'], x1, x2)
    dat['models']['COPh'] = linalg_hp(dat['CAPh']/dat['Wh'], x1, x2)

    return dat


def linalg_hp(y, x1, x2):

    indx = np.where(~np.isnan(y))[0]

    y = y[indx]
    x1 = x1[indx]
    x2 = x2[indx]

    B = y
    A = np.column_stack([np.ones(len(y)), x1, x1**2, x2, x2**2, x1*x2])

    A = np.linalg.lstsq(A, B)[0]

    return A


class HeatPump(object):
    """
    TinHP: Temperature of the fluid entering the HP in Celsius.
    Nhp: Number of heat pump
    hpname:  Name of the heat pump
    """
    def __init__(self):                # attributes that are linked with the UI

        self.__initAttr__()

        self.TinHP = {'cooling': 28, 'heating': 0}
        self.Nhp = 1
        self.hpname = self.hpDB.keys()[0]

    def __initAttr__(self):        # attributes that are not linked with the UI

        hpfile = os.path.dirname(os.path.realpath(__file__))
        hpfile = os.path.join(hpfile, 'GHP', 'hp_database.npy')
        self.hpDB = np.load(hpfile).item()

        # Building thermal load in kW (+ for cooling, - for heating):
        self.qbat = {'cooling': 16.5, 'heating': 14.5}

        self.fluid = 'water'  # Heat carrier fluid type
        self.fr = 0           # antifreeze volumetric fraction

        self.Vftot = {}  # Total volumetric flow in the system in L/s
        self.Vftot['cooling'] = 0.05 * 16.5
        self.Vftot['heating'] = 0.05 * 14.5

        self.Tg = 12  # undisturbed ground temperature

    @property
    def hpdata(self):
        return self.hpDB[self.hpname]

    @property
    def ToutHP(self):
        """Temperature of the water leaving the heat pump in ºC."""
        ToutHP = {}
        for mode in ['cooling', 'heating']:
            ToutHP[mode] = self.calcul_ToutHP(mode)

        return ToutHP

    @property
    def Tm(self):
        """Fluid mean temperature through the heat pump in ºC"""
        ToutHP = self.ToutHP
        TinHP = self.TinHP

        Tm = {}
        for mode in ['cooling', 'heating']:
            Tm[mode] = (TinHP[mode] + ToutHP[mode])/2

        return Tm

    @property
    def Vhp(self):
        """Volumetric flowrate per HP in L/s"""
        vhp = {'heating': self.Vftot['heating']/self.Nhp,
               'cooling': self.Vftot['cooling']/self.Nhp}
        return vhp

    def calcul_ToutHP(self, mode):

        # Calculate fluid properties :

        hcfluid = HCFluid(self.fluid, self.TinHP[mode], self.fr)
        rhof = hcfluid.rho
        cpf = hcfluid.cp

        # Calculate ground load :

        qbat = self.qbat[mode]
        COP = self.get_COP(mode)

        if mode == 'cooling':
            qgnd = qbat * (COP+1)/COP
        elif mode == 'heating':
            qgnd = -qbat * (COP-1)/COP

        # Calculate outflow fluid temperature :

        ToutHP = self.TinHP[mode] + qgnd/(self.Vftot[mode]*rhof*cpf) * 10**6

        return ToutHP

    def interp(self, varname, ewt, gpm):
        # ewt: entering water temperature in the HP (ºC)
        # gpm: columetric flowrate in the HP (L/s)

        A = self.hpDB[self.hpname]['models'][varname]

        var = (A[0] +
               A[1]*ewt + A[2]*ewt**2 +
               A[3]*gpm + A[4]*gpm**2 +
               A[5]*ewt*gpm
               )

        # Anti-freeze correction factor :

        afcorr = np.interp(self.fr,
                           corrtbl_afreeze[self.fluid]['fr'],
                           corrtbl_afreeze[self.fluid][varname]
                           )

        return var*afcorr

    def in_table(self, varname, x1, y1):
        x = self.hpDB[self.hpname]['EWT']
        y = self.hpDB[self.hpname]['GPM']
        z = self.hpDB[self.hpname][varname]

        # remove nan values :
        indx = np.where(~np.isnan(z))[0]
        x = x[indx]
        y = y[indx]

        # Check if point is inside the table :
        # http://stackoverflow.com/a/16898636/4481445

        hull = Delaunay(np.vstack((x, y)).T)

        return bool(hull.find_simplex((x1, y1)) >= 0)

    def get_flowRange(self):
        vmax = np.max(self.hpDB[self.hpname]['GPM']) * self.Nhp
        vmin = np.min(self.hpDB[self.hpname]['GPM']) * self.Nhp
        return vmin, vmax

    def get_COP(self, mode):
        if mode == 'cooling':
            return self.interp('COPc', self.TinHP[mode], self.Vhp[mode])
        elif mode == 'heating':
            return self.interp('COPh', self.TinHP[mode], self.Vhp[mode])

    def get_CAP(self, mode):
        if mode == 'cooling':
            return self.interp('CAPc', self.TinHP[mode], self.Vhp[mode])
        elif mode == 'heating':
            return self.interp('CAPh', self.TinHP[mode], self.Vhp[mode])
#
#    def get_Whp(self, mode):
#        if mode == 'cooling':
#            return self.interp('Wc', self.TinHP[mode], self.Vhp[mode])
#        elif mode == 'heating':
#            return self.interp('Wh', self.TinHP[mode], self.Vhp[mode])
#
#    def get_WPD(self, mode):
#        return self.interp('WPD', self.TinHP[mode], self.Vhp[mode])


def eval_model():

    plt.close('all')

    w = HeatPump()

    for name in ['TCHV072', 'TCHV096', 'TCHV120', 'TCHV160', 'TCHV192',
                 'TCHV240', 'TCHV300', 'Generic']:
        w.hpname = name

        y = w.hpdata['CAPh']/w.hpdata['Wh']
        x1 = w.hpdata['EWT']
        x2 = w.hpdata['GPM']
        A = linalg_hp(y, x1, x2)

        yp = (A[0] +
              A[1]*x1 + A[2]*x1**2 +
              A[3]*x2 + A[4]*x2**2 +
              A[5]*x1*x2
              )

        fig, ax = plt.subplots()
        fig.canvas.set_window_title(name)
        ax.plot([np.floor(np.min(yp)), np.ceil(np.max(yp))],
                [np.floor(np.min(yp)), np.ceil(np.max(yp))], '--k')
        ax.plot(y, yp, 'o')

    plt.show()


def plot_cop():
    plt.close('all')

    w = HeatPump()

    w.hpname = 'Generic'
    w.hpname = 'TCHV300'
    # y = w.hpdata['CAPh']/w.hpdata['Wh']
    y = w.hpdata['CAPc']/w.hpdata['Wc']
    x1 = w.hpdata['EWT']
    x2 = w.hpdata['GPM']
    A = linalg_hp(y, x1, x2)

    print(w.get_flowRange())

    fig, ax = plt.subplots()

    ewt = np.arange(-1, 35)
#    flowrange = np.arange(w.get_flowRange()[0], w.get_flowRange()[1], 0.25)
    flowrange = np.arange(2.37, 4.73, 0.25)

    for val in flowrange:
        gpm = np.ones(len(ewt)) * val

        cop = (A[0] +
               A[1]*ewt + A[2]*ewt**2 +
               A[3]*gpm + A[4]*gpm**2 +
               A[5]*ewt*gpm
               )
        ax.plot(ewt, cop, '-')

    plt.show()


if __name__ == '__main__':
    import matplotlib.pyplot as plt

#    eval_model()
#    plot_cop()
#
    w = HeatPump()
    w.show()
#    w.setFixedSize(w.size())
    w.setCurrentUnitSystem('SI')

    w.set_Vftot({'cooling': 12.5/15.8503230745/1000,
                 'heating': 12.5/15.8503230745/1000})
