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


def build_database(dirname):
    """
    Build a database of heat pumps performance data from a list of
    formatted csv files located within the specified directory name.
    """
    db = OrderedDict()
    for f in os.listdir(dirname):
        if f.endswith('.csv'):
            try:
                dat = load_HP_table_fromfile(os.path.join(dirname, f))
                name = dat['model']
                db[name] = dat
            except Exception:
                print('unable to load data from %s' % f)

    filename = osp.join(dirname, 'hp_database.npy')
    np.save(filename, db)

    return filename


def load_HP_table_fromfile(filename):
    """
    Load heat pump performance data from the specified csv file and format
    the data in a dictionary.

    EWT  : entering temperature in ºC
    GPM  : volumetric flow in L/s
    WPDc : pressure drop in the coil for the cooling mode in Pa
    WPDh : pressure drop in the coil for the heating mode in Pa
    Wc   : work of the heat pump for the cooling mode in kW
    Wh   : work of the heat pump for the heating mode in kW
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
                data['Nominal CAP'] = float(row[1])*3.51685
            elif row[0] == 'EWT (F)':
                A = np.array(reader[i+1:]).astype(float)

                data['EWT'] = (A[:, 0]-32)/1.8
                data['GPM'] = A[:, 1]/15.8503230745

                data['WPDc'] = A[:, 2]*6894.76
                data['WPDh'] = A[:, 5]*6894.76

                data['CAPc'] = A[:, 3]*0.293071
                data['CAPh'] = A[:, 6]*0.293071

                data['Wc'] = A[:, 4]
                data['Wh'] = A[:, 7]
                break

    data['models'] = {}

    x1 = data['EWT']
    x2 = data['GPM']

    data['models']['CAPc'] = linalg_hp(data['CAPc'], x1, x2)
    data['models']['CAPh'] = linalg_hp(data['CAPh'], x1, x2)

    data['models']['COPc'] = linalg_hp(data['CAPc']/data['Wc'], x1, x2)
    data['models']['COPh'] = linalg_hp(data['CAPh']/data['Wh'], x1, x2)

    return data


def linalg_hp(y, x1, x2):

    indx = np.where(~np.isnan(y))[0]

    y = y[indx]
    x1 = x1[indx]
    x2 = x2[indx]

    B = y
    A = np.column_stack([np.ones(len(y)), x1, x1**2, x2, x2**2, x1*x2])

    A = np.linalg.lstsq(A, B)[0]

    return A


if __name__ == '__main__':
    build_database()

#     import matplotlib.pyplot as plt

# #    eval_model()
# #    plot_cop()
# #
#     w = HeatPump()
# #    w.setFixedSize(w.size())
#     w.setCurrentUnitSystem('SI')

#     w.set_Vftot({'cooling': 12.5/15.8503230745/1000,
#                  'heating': 12.5/15.8503230745/1000})
