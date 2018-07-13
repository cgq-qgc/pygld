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
    """
    This is the base class from which all other ground-loop heat exchanger
    materials should derived.
    """
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
    def init_as(cls, category_name, material_index):
        """
        Initializes the thermal properties based on the specified
        category and material index.

        This is a convenience function; the member variables can also be
        initialized manually.
        """
        predefmats = PREDEFINED_MATERIALS.get(category_name, None)

        material_names = list(predefmats.keys())

        # Ensure there is no IndexError when accessing the material database.
        material_index = min(max(material_index, 0), len(material_names)-1)
        material_name = material_names[material_index]

        instance = cls(*predefmats[material_name])
        instance._material = material_name
        instance._category = category_name
        return instance


class GroutMaterial(BaseMaterial):
    """
    The :class:`GroutMaterial` class holds all the thermophysical properties
    relative to the grout that are required for the modeling of ground-loop
    heat exchanger systems.

    A database of predefined grout materials is available and can be obtained
    or printed with the :meth:`~GroutMaterial.get_predefined_materials` and
    :meth:`~GroutMaterial.print_predefined_materials` methods.
    The :class:`GroutMaterial` can be initialized directly from a predefined
    material with the :meth:`~GroutMaterial.init_as`, passing as argument the
    index of the material in the database.

    An `Example`_ is available at the end of this section.
    """

    def __init__(self, kth=None, Cp=None):
        super().__init__(kth, Cp)
        self._category = BaseMaterial.Grout

    @classmethod
    def get_predefined_materials(cls):
        """Return a dict containing the predefined :class:`GroutMaterial`."""
        return super().get_predefined_materials(BaseMaterial.Grout)

    @classmethod
    def print_predefined_materials(cls, end='\n'):
        """Print all the predefined :class:`GroutMaterial`."""
        super().print_predefined_materials(BaseMaterial.Grout, end)

    @classmethod
    def init_as(cls, material_index):
        """
        Initialize the thermal properties of the :class:`GroutMaterial`
        based on the specified material index.

        This is a convenience function; the member variables can also be
        initialized manually.
        """
        return super().init_as(BaseMaterial.Grout, material_index)


class PipeMaterial(BaseMaterial):
    """
    The :class:`PipeMaterial` class holds all the thermophysical properties
    relative to the pipes that are required for the modeling of ground-loop
    heat exchanger systems.

    A database of predefined pipe materials is available and can be obtained
    or printed with the :meth:`~PipeMaterial.get_predefined_materials` and
    :meth:`~PipeMaterial.print_predefined_materials` methods.
    The :class:`PipeMaterial` can be initialized directly from a predefined
    material with the :meth:`~PipeMaterial.init_as`, passing as argument the
    index of the material in the database.

    An `Example`_ is available at the end of this section.
    """

    def __init__(self, kth=None, Cp=None):
        super().__init__(kth, Cp)
        self._category = BaseMaterial.Pipe

    @classmethod
    def get_predefined_materials(cls):
        """Return a dict containing the predefined :class:`PipeMaterial`."""
        return super().get_predefined_materials(BaseMaterial.Pipe)

    @classmethod
    def print_predefined_materials(cls, end='\n'):
        """Print the predefined :class:`PipeMaterial`."""
        super().print_predefined_materials(BaseMaterial.Pipe, end)

    @classmethod
    def init_as(cls, material_index):
        """
        Initialize the thermal properties of the :class:`PipeMaterial`
        based on the specified material index.

        This is a convenience function; the member variables can also be
        initialized manually.
        """
        return super().init_as(BaseMaterial.Pipe, material_index)


class GroundMaterial(BaseMaterial):
    """
    The :class:`GroundMaterial` class holds all the thermophysical properties
    relative to the ground that are required for the modeling of ground-loop
    heat exchanger systems.

    A database of predefined ground materials is available and can be obtained
    or printed with the :meth:`~GroundMaterial.get_predefined_materials` and
    :meth:`~GroundMaterial.print_predefined_materials` methods.
    The :class:`GroundMaterial` can be initialized directly from a predefined
    material with the :meth:`~GroundMaterial.init_as`, passing as argument the
    index of the material in the database.

    An `Example`_ is available at the end of this section.
    """

    def __init__(self, kth=None, Cp=None):
        super().__init__(kth, Cp)
        self._category = BaseMaterial.Ground

    @classmethod
    def get_predefined_materials(cls):
        """Return a dict containing the predefined :class:`GroundMaterial`."""
        return super().get_predefined_materials(BaseMaterial.Ground)

    @classmethod
    def print_predefined_materials(cls, end='\n'):
        """Print the predefined :class:`GroundMaterial`."""
        super().print_predefined_materials(BaseMaterial.Ground, end)

    @classmethod
    def init_as(cls, material_index):
        """
        Initialize the thermal properties of the :class:`GroundMaterial`
        based on the specified material index.

        This is a convenience function; the member variables can also be
        initialized manually.
        """
        return super().init_as(BaseMaterial.Ground, material_index)


if __name__ == '__main__':
    grout = GroutMaterial.init_as(28)

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
