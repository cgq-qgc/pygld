Grout Material
---------------------------------

.. currentmodule:: pygld

.. autoclass:: pygld.GroutMaterial
   :no-members:
   :no-undoc-members:

Thermal Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoattribute:: pygld.GroutMaterial.Cp
.. autoattribute:: pygld.GroutMaterial.kth
.. autoattribute:: pygld.GroutMaterial.al

Utilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: pygld.GroutMaterial.get_predefined_materials
.. automethod:: pygld.GroutMaterial.print_predefined_materials
.. automethod:: pygld.GroutMaterial.init_as

Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Import and print the list of predefined :class:`~pygld.GroutMaterial`::

    >>> from pygld import GroundMaterial
    >>> GroundMaterial.print_predefined_materials()
    Predefined materials for the 'Grout':
        0 - Bentonite 12% ........ kth=0.700 W/m·k, Cp=3900 J/m³·K
        1 - Bentonite/Sand 12%/15% kth=1.500 W/m·k, Cp=3400 J/m³·K
        2 - Water ................ kth=0.600 W/m·k, Cp=4200 J/m³·K
        
Initialize :class:`~pygld.GroutMaterial` from the predefined material at index
#1 and print the results::

    >>> grout = GroundMaterial.init_as(1)
    >>> print(grout)
    Grout material: Bentonite/Sand 12%/15%
    Thermal conductivity: 1.500 W/m·k
    Volumetric heat capacity: 3400 J/m³·K
    Thermal diffusivity: 4.41e-04 m²/s
    
Change the thermal conductivity for another value and print the results::

    >>> grout.kth = 1.75
    >>> print(grout)
    Grout material: User Defined
    Thermal conductivity: 1.750 W/m·k
    Volumetric heat capacity: 3400 J/m³·K
    Thermal diffusivity: 5.15e-04 m²/s

Note that the grout material is now labeled as 'User defined'.

It is also possible to initialize an empty :class:`~pygld.GroutMaterial` and
set the value of the properties manually::

    >>> grout = GroutMaterial()
    >>> grout.kth = 2
    >>> grout.Cp = 3500
    >>> print(grout)
    Grout material: User Defined
    Thermal conductivity: 2.000 W/m·k
    Volumetric heat capacity: 3500 J/m³·K
    Thermal diffusivity: 5.71e-04 m²/s