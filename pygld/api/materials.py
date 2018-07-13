# -*- coding: utf-8 -*-

# Copyright © 2017-2018 Jean-Sébastien Gosselin
# https://github.com/jnsebgosselin/pygld
#
# This file is part of PyGLD.
# Licensed under the terms of the MIT License.

# ---- Standard imports

# ---- Third party imports

from collections import OrderedDict

# ---- Local imports

PREDEFINED_MATERIALS = OrderedDict({
    'Ground': OrderedDict({
        'Amphibolite': (2.9, 2600),
        'Andesite': (2.2, 2400),
        'Anhydrite': (4.1, 2000),
        'Aplite': (3.1, 2400),
        'Arkose': (2.9, 2000),
        'Basalt': (1.7, 2400),
        'Breccia': (2.8, 2100),
        'Clay, dry': (0.4, 1600),
        'Clay, moist': (1.6, 2400),
        'Claystone': (2.2, 2300),
        'Coal': (0.3, 1800),
        'Concrete': (1.6, 1800),
        'Conglomerate': (2.8, 2100),
        'Diorite': (2.6, 2900),
        'Dolomite': (3.2, 2500),
        'Dunite': (4.2, 2900),
        'Eclogite': (2.9, 3100),
        'Gabbro': (1.9, 2600),
        'Gneiss': (2.9, 2100),
        'Granite': (3.4, 2400),
        'Granodiorite': (3.3, 2600),
        'Gravel, dry': (0.4, 1500),
        'Gravel, saturated': (1.8, 2400),
        'Gypsum': (1.6, 2000),
        'Gyttja (sapropel)': (0.6, 2000),
        'Ice': (2.2, 1900),
        'Ignimbrite': (3.0, 2100),
        'Lamprophyre': (2.6, 2100),
        'Limestone, massive': (2.8, 2300),
        'Limestone, marly': (2.2, 2300),
        'Limestone, oolitic': (2.4, 2300),
        'Marble': (2.6, 2000),
        'Marl': (2.1, 2300),
        'Marl, clayey': (2.0, 2200),
        'Marl, dolomitic': (2.2, 2300),
        'Metaquartzite': (5.8, 2100),
        'Micashist': (2.0, 2200),
        'Peat': (0.4, 2200),
        'Pegmatite': (3.0, 2100),
        'Peridotite': (4.0, 2700),
        'Quartzite': (6.0, 2100),
        'Rhyolite': (3.3, 2100),
        'Salt': (5.4, 1200),
        'Sand, dry': (0.4, 1400),
        'Sand, dry, compacted': (1.2, 1700),
        'Sand, moist': (1.0, 1800),
        'Sand, saturated': (2.4, 2500),
        'Sand, frozen': (2.0, 1500),
        'Sandstone': (2.3, 2000),
        'Serpentinite': (3.0, 2200),
        'Shale': (2.1, 2300),
        'Silt, dry': (0.4, 1600),
        'Silt, moist': (1.8, 2200),
        'Siltstone': (2.4, 2300),
        'Syenite': (2.6, 2400),
        'Till (moraine)': (2.0, 2100),
        'Tonalite': (2.7, 2400),
        'Trachyte': (2.8, 2100),
        'Tuff': (1.1, 1100),
        'Water': (0.6, 4200)
        }),
    'Grout': OrderedDict({
        'Bentonite 12%': (0.7, 3900),
        'Bentonite/Sand 12%/15%': (1.5, 3400),
        'Water': (0.6, 4200)
        }),
    'Pipe': OrderedDict({
        'HDPE': (0.4, 1500),
        'Geoperformx': (0.7, 1500)
        })
     })


class BaseMaterial(object):
    Grout = 'Grout'
    Pipe = 'Pipe'
    Ground = 'Ground'

    def __init__(self, kth=None, Cp=None):
        self._category = None
        self._material = None
        self.kth = kth
        self.Cp = Cp

    def __str__(self):
        kth = 'None' if self.kth is None else '{:.3f} W/m·k'.format(self.kth)
        Cp = 'None' if self.Cp is None else '{:.0f} J/m³·K'.format(self.Cp)
        al = 'None' if self.al is None else '{:.2e} m²/s'.format(self.al)

        category = ('Material' if self._category is None else
                    self._category + ' material')
        material = self._material or 'User Defined'
        str_ = '%s: %s\n' % (category, material)
        str_ += 'Thermal conductivity: %s' % kth
        str_ += '\nVolumetric heat capacity: %s' % Cp
        str_ += '\nThermal diffusivity: %s' % al
        return str_

    @property
    def Cp(self):
        """Volumetric heat capacity in J/m³·K.

        Get or set the volumetric heat capacity of the material as a single
        positive float value.
        """
        return self._Cp

    @Cp.setter
    def Cp(self, x):
        self._Cp = None if x is None else abs(float(x))
        self._material = None

    @property
    def kth(self):
        """Thermal conductivity in W/m·K.

        Get or set the thermal conductivity of the material as a single
        positive float value.
        """
        return self._kth

    @kth.setter
    def kth(self, x):
        self._kth = None if x is None else abs(float(x))
        self._material = None

    @property
    def al(self):
        """Thermal diffusivity in m²/s.

        Return the thermal diffusivity of the material as a single
        positive float value calculated as:

        .. math::
            al[i] = \\frac{kth[i]}{cp[i] \\cdot rho[i]}

        where :math:`kth` is the thermal conductivity in W/m·k and
        :math:`Cp` is the volumetric heat capacity of the material in J/m³·K.
        """
        if self.kth is None or self.Cp is None:
            return None
        else:
            return self.kth / self.Cp

    # ---- Utilities

    @classmethod
    def get_predefined_materials(cls, category):
        return PREDEFINED_MATERIALS.get(category, None)

    @classmethod
    def print_predefined_materials(cls, categories=None, end='\n'):
        if categories is None:
            categories = list(PREDEFINED_MATERIALS.keys())
        elif isinstance(categories, str):
            categories = [categories]

        str_ = ''
        for category in categories:
            str_ = str_ + '\n' if str_ else str_
            str_ += "Predefined materials for the '%s':" % category

            predefmats = PREDEFINED_MATERIALS.get(category, None)
            if predefmats is not None:
                maxlen = len(max(list(predefmats.keys()), key=len)) - 1
                fmt = '\n{:>5d} - {:<}{:s} kth={:.3f} W/m·k, Cp={:d} J/m³·K'
                for i, (key, value) in enumerate(predefmats.items()):
                    pad = '.' * (maxlen-len(key))
                    pad = ' ' + pad if pad else pad
                    str_ += fmt.format(i, key, pad, value[0], value[1])
            else:
                str_ += ' None'
        print(str_, end=end)

    @classmethod
    def init_as(cls, category, material):
        """
        Initializes the thermal properties based on the specified
        category and material. The material can be either a str or an index.
        The list of predefined categories and materials can be printed with the
        `Material.print_predefined_materials` method. If no argument is
        provided, the materials will be printed for all available categories.

        This is a convenience function; the member variables can also be
        initialized manually.
        """
        predefmats = PREDEFINED_MATERIALS.get(category, None)
        material = (list(predefmats.keys())[material] if
                    isinstance(material, int) else material)

        instance = cls(*predefmats[material])
        instance._material = material
        instance._category = category
        return instance


class GroutMaterial(BaseMaterial):

    @classmethod
    def get_predefined_materials(cls):
        """Return a dict containing the predefined :class:`GroutMaterial`."""
        return super().get_predefined_materials(BaseMaterial.Grout)

    @classmethod
    def print_predefined_materials(cls, end='\n'):
        """Print all the predefined :class:`GroutMaterial`."""
        super().print_predefined_materials(BaseMaterial.Grout, end)

    @classmethod
    def init_as(cls, material):
        """
        Initializes the thermal properties of the :class:`GroutMaterial`
        based on the specified material index.
        The methods :meth:`~GroutMaterial.print_predefined_materials` and
        :meth:`~GroutMaterial.get_predefined_materials` can be used to,
        respectively, print or get in a dict the predefined
        :class:`GroutMaterial` that are availables.

        This is a convenience function; the member variables can also be
        initialized manually.
        """
        return super().init_as(BaseMaterial.Grout, material)


class PipeMaterial(BaseMaterial):

    @classmethod
    def get_predefined_materials(cls):
        """Return a dict containing the predefined :class:`PipeMaterial`."""
        return super().get_predefined_materials(BaseMaterial.Pipe)

    @classmethod
    def print_predefined_materials(cls, end='\n'):
        """Print the predefined :class:`PipeMaterial`."""
        super().print_predefined_materials(BaseMaterial.Pipe, end)

    @classmethod
    def init_as(cls, material):
        """
        Initializes the thermal properties of the :class:`PipeMaterial`
        based on the specified material index.
        The methods :meth:`~PipeMaterial.print_predefined_materials` and
        :meth:`~PipeMaterial.get_predefined_materials` can be used to,
        respectively, print or get in a dict the predefined
        :class:`PipeMaterial` that are availables.

        This is a convenience function; the member variables can also be
        initialized manually.
        """
        return super().init_as(BaseMaterial.Pipe, material)


class GroundMaterial(BaseMaterial):

    @classmethod
    def get_predefined_materials(cls):
        """Return a dict containing the predefined :class:`GroundMaterial`."""
        return super().get_predefined_materials(BaseMaterial.Ground)

    @classmethod
    def print_predefined_materials(cls, end='\n'):
        """Print the predefined :class:`GroundMaterial`."""
        super().print_predefined_materials(BaseMaterial.Ground, end)

    @classmethod
    def init_as(cls, material):
        """
        Initializes the thermal properties of the :class:`GroundMaterial`
        based on the specified material index.
        The methods :meth:`~GroundMaterial.print_predefined_materials` and
        :meth:`~GroundMaterial.get_predefined_materials` can be used to,
        respectively, print or get in a dict the predefined
        :class:`GroundMaterial` that are availables.

        This is a convenience function; the member variables can also be
        initialized manually.
        """
        return super().init_as(BaseMaterial.Ground, material)


if __name__ == '__main__':
    # ---- Example Grout

    GroutMaterial.print_predefined_materials(end='\n\n')

    grout = GroutMaterial.init_as(1)
    print(grout, end='\n\n')

    grout.kth = 1.75
    print(grout, end='\n\n')

    grout = GroutMaterial()
    grout.kth = 2
    grout.Cp = 3500
    print(grout, end='\n\n')

    # ---- Example Ground

    GroundMaterial.print_predefined_materials(end='\n\n')

    ground = GroundMaterial.init_as(14)
    print(ground, end='\n\n')

    ground.kth = 3.35
    print(ground, end='\n\n')

    ground = GroundMaterial()
    ground.kth = 2.27
    ground.Cp = 3278
    print(ground, end='\n\n')

    # ---- Example Pipe

    PipeMaterial.print_predefined_materials(end='\n\n')
    pipe = PipeMaterial.init_as(1)
    print(pipe, end='\n\n')

    pipe.kth = 0.67
    print(pipe, end='\n\n')

    pipe = PipeMaterial()
    pipe.kth = 0.53
    pipe.Cp = 1767
    print(pipe)
