# -*- coding: utf-8 -*-

# Copyright © 2018 Jean-Sébastien Gosselin
# https://github.com/jnsebgosselin/qwatson
#
# This file is part of PyGLD.
# Licensed under the terms of the GNU General Public License.

# ---- Standard imports

import os
import os.path as osp
import csv
from collections import OrderedDict

# ---- Third party imports

import numpy as np

# ---- Local imports

from .__init__ import __datadir__
from .maths import multi_polyfit2rd


def load_heatpump_database(dirname=None):
    """
    Load the database from the specified directory or create it if it does
    not exist.
    """
    dirname = __datadir__ if dirname is None else dirname

    hpfile = os.path.join(dirname, 'hp_database.npy')
    if not osp.exists(hpfile):
        build_heatpump_database(dirname)

    return np.load(hpfile).item()


def build_heatpump_database(dirname):
    """
    Build a database of heat pumps performance data from a list of
    formatted csv files located within the specified directory name.
    """
    db = OrderedDict()

    files = [osp.join(dirname, f) for
             f in os.listdir(dirname) if f.endswith('.hp')]
    for file in files:
        try:
            data = load_heatpump_table_fromfile(file)
            db[data['model']] = data
        except Exception:
            print('unable to load data from %s' % file)

    filename = osp.join(dirname, 'hp_database.npy')
    np.save(filename, db)

    return filename


def load_heatpump_table_fromfile(filename):
    """
    Load heat pump performance data from the specified csv file and format
    the data in a dictionary.

    EWT  : entering temperature in ºC
    GPM  : volumetric flow in L/s
    WPDc : pressure drop in the coil for the cooling mode in Pa
    WPDh : pressure drop in the coil for the heating mode in Pa
    Wc   : work of the heat pump for the cooling mode in kW
    Wh   : work of the heat pump for the heating mode in kW

    Note that the data in the input csv files are in imperial units, so we need
    to convert the values in SI format.
    """
    with open(filename, 'r', encoding='utf8') as f:
        reader = list(csv.reader(f))

        data = {'filename': filename}
        for i, row in enumerate(reader):
            if row[0] == 'Manufacturer':
                data['manufacturer'] = row[1]
            elif row[0] == 'Model':
                data['model'] = row[1]
                print('Loading %s data sheet...' % row[1])
            elif row[0] == 'Nominal CAP (tons)':
                data['Nominal CAP'] = float(row[1]) * 3.51685
            elif row[0] == 'EWT (F)':
                A = np.array(reader[i+1:]).astype(float)

                data['EWT'] = (A[:, 0]-32) / 1.8
                data['GPM'] = A[:, 1] / 15.8503230745

                data['WPDc'] = A[:, 2] * 6894.76
                data['WPDh'] = A[:, 5] * 6894.76

                data['CAPc'] = A[:, 3] * 0.293071
                data['CAPh'] = A[:, 6] * 0.293071

                data['Wc'] = A[:, 4]
                data['Wh'] = A[:, 7]

                data['COPc'] = data['CAPc'] / data['Wc']
                data['COPh'] = data['CAPh'] / data['Wh']
                break

    # Use the performance data to build an equation-fit model of the form :

    # y = a1 + a2*EWT + a3*EWT^2 + a4*GPM, a5*GPM^2 + a6*EWT*GPM

    # where ai are the coefficients, EWT is the entering temperature in the
    # heat pump, GPM is the flow rate and y is the variable that we want to
    # model, namely the coefficient of performance and capacity in cooling
    # and heating mode.

    # This is based on Equations 19 and 20 in :
    # Jin, H. and J.D.Spitler, 2002. A parameter estimation based model of
    # water-to-water heat pumps for use in energy calculation programs.
    # ASHRAE Transactions: Research, 108(1): 3-17

    data['eqfit_models'] = {}

    x1 = data['EWT']
    x2 = data['GPM']

    data['eqfit_models']['CAPc'] = multi_polyfit2rd(data['CAPc'], x1, x2)
    data['eqfit_models']['CAPh'] = multi_polyfit2rd(data['CAPh'], x1, x2)

    data['eqfit_models']['COPc'] = multi_polyfit2rd(data['COPc'], x1, x2)
    data['eqfit_models']['COPh'] = multi_polyfit2rd(data['COPh'], x1, x2)

    return data


if __name__ == '__main__':
    input_fname = osp.join(__datadir__, 'TCHV072.hp')
    data_TCHV072 = load_heatpump_table_fromfile(input_fname)
