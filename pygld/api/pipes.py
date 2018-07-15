# -*- coding: utf-8 -*-

# Copyright © 2017-2018 Jean-Sébastien Gosselin
# https://github.com/jnsebgosselin/pygld
#
# This file is part of PyGLD.
# Licensed under the terms of the MIT License.

# ---- Standard imports

# ---- Third party imports

from numpy import pi

# ---- Local imports

from pygld.lib.thermal_resistances import calcul_Rcond_pipe
from pygld import PipeMaterial


class Pipe(object):
    """
    The :attr:`~pygld.Pipe` class holds the geometry and material properties
    of the pipes used in the construction of ground-loop heat exchangers.
    The inner and outer diameters of the pipe in cm must be passed as arguments
    when instantiating :attr:`~pygld.Pipe`. These values can afterward be
    accessed or changed with the :attr:`.di` and :attr:`.do` properties

    The thermophysical properties of the pipe are held in a
    :class:`~pygld.PipeMaterial` object that can be accessed with the
    :attr:`.material` property or set with :meth:`.set_material`.
    The optional argument 'material' can be used to pass a custom
    :class:`~pygld.PipeMaterial` when instantiating a :attr:`~pygld.Pipe`.
    If no :class:`~pygld.PipeMaterial` is provided when instantiating
    :attr:`~pygld.Pipe`, the pipe's :attr:`.material` will be set to that of
    a standard HDPE pipe by default, with a thermal conductivity of 0.4 W/m·k
    and a volumetric heat capacity of 1500 J/m³·K. These values can be changed
    afterward by setting directly the properties :attr:`PipeMaterial.kth` and
    :attr:`PipeMaterial.Cp` of the pipe's :attr:`.material`.

    An `Example`_ is available at the end of this section.
    """
    def __init__(self, di, do, material=None):
        self._di = di
        self._do = do

        material = PipeMaterial.init_as(0) if material is None else material
        self.set_material(material)

    def __str__(self):
        str_ = 'Inner diameter: %0.2f cm' % (self.di)
        str_ += '\nOuter diameter: %0.2f cm' % (self.do)
        str_ += '\nWall thickness: %0.2f cm' % (self.wt)
        str_ += '\nInner area: %0.2f cm²' % (self.Ai)
        str_ += '\nOuter area: %0.2f cm²\n' % (self.Ao)
        str_ += self._material.__str__()
        str_ += '\nConductive thermal resistance: %0.5f m·K/W' % self.Rcond
        return str_

    # ---- Geometry

    @property
    def di(self):
        """Get or set the inner diameter of the pipe in cm."""
        return self._di

    @di.setter
    def di(self, x):
        return self._di

    @property
    def do(self):
        """Get or set the outer diameter of the pipe in cm."""
        return self._do

    @do.setter
    def do(self, x):
        return self._do

    @property
    def do(self):
        """Get or set the outer diameter of the pipe in cm."""
        return self._do

    @property
    def wt(self):
        """
        Return the pipe wall thickness in cm calculated as:

        .. math::
            wt = (do - di) / 2

        where :math:`do` and :math:`di` are, respectively, the outer and inner
        diameter of the pipe in cm.
        """
        return (self.do - self.di)/2

    @property
    def Ai(self):
        """Return the inner area of the pipe in cm² calculated as:

        .. math::
            Ai = \\frac{\\pi \\cdot di^2}{4}

        where :math:`di` is the inner diameter of the pipe.
        """
        return pi*self.di**2/4

    @property
    def Ao(self):
        """Return the inner area of the pipe in cm² calculated as:

        .. math::
            Ao = \\frac{\\pi \\cdot do^2}{4}

        where :math:`do` is the outer diameter of the pipe.
        """
        return pi*self.do**2/4

    # ---- Material

    @property
    def material(self):
        """
        Return the :class:`~pygld.PipeMaterial` of the pipe.
        """
        return self._material

    def set_material(self, material):
        """
        Set the :class:`~pygld.PipeMaterial` of the pipe.
        """
        if isinstance(material, PipeMaterial):
            self._material = material
        else:
            raise TypeError("The 'material' argument must be an instance"
                            "of 'PipeMaterial'")

    @property
    def Rcond(self):
        """Return the conductive thermal resistance of the pipe in m·K/W
        calculated as:

        .. math::
            Rcond = \\frac{log(do/di)}{2 \\cdot \\pi \\cdot kth}

        where :math:`di` and :math:`do` are, repectively, the inner and outer
        diameter of the pipe in m and :math:`kth` is its thermal conductivity
        in W/m·k.
        """
        return calcul_Rcond_pipe(self._material.kth,
                                 self.di/2/100, self.do/2/100)


if __name__ == '__main__':
    pipe = Pipe(di=3.39852, do=4.2164)
    print(pipe)

    pipe.material.kth = 0.43
    pipe.material.Cp = 1540
    print(pipe.material)
