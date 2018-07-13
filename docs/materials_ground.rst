Ground Material
---------------------------------

.. currentmodule:: pygld

.. autoclass:: pygld.GroundMaterial
   :no-members:
   :no-undoc-members:

Thermal Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoattribute:: pygld.GroundMaterial.Cp
.. autoattribute:: pygld.GroundMaterial.kth
.. autoattribute:: pygld.GroundMaterial.al

Utilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: pygld.GroundMaterial.get_predefined_materials
.. automethod:: pygld.GroundMaterial.print_predefined_materials
.. automethod:: pygld.GroundMaterial.init_as

Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Import and print the list of predefined :class:`~pygld.GroundMaterial`::

    >>> from pygld import GroundMaterial
    >>> GroundMaterial.print_predefined_materials()
    Predefined materials for the 'Ground':
        0 - Amphibolite ........ kth=2.900 W/m·k, Cp=2600 J/m³·K
        1 - Andesite ........... kth=2.200 W/m·k, Cp=2400 J/m³·K
        2 - Anhydrite .......... kth=4.100 W/m·k, Cp=2000 J/m³·K
        3 - Aplite ............. kth=3.100 W/m·k, Cp=2400 J/m³·K
        4 - Arkose ............. kth=2.900 W/m·k, Cp=2000 J/m³·K
        5 - Basalt ............. kth=1.700 W/m·k, Cp=2400 J/m³·K
        6 - Breccia ............ kth=2.800 W/m·k, Cp=2100 J/m³·K
        7 - Clay, dry .......... kth=0.400 W/m·k, Cp=1600 J/m³·K
        8 - Clay, moist ........ kth=1.600 W/m·k, Cp=2400 J/m³·K
        9 - Claystone .......... kth=2.200 W/m·k, Cp=2300 J/m³·K
       10 - Coal ............... kth=0.300 W/m·k, Cp=1800 J/m³·K
       11 - Concrete ........... kth=1.600 W/m·k, Cp=1800 J/m³·K
       12 - Conglomerate ....... kth=2.800 W/m·k, Cp=2100 J/m³·K
       13 - Diorite ............ kth=2.600 W/m·k, Cp=2900 J/m³·K
       14 - Dolomite ........... kth=3.200 W/m·k, Cp=2500 J/m³·K
       15 - Dunite ............. kth=4.200 W/m·k, Cp=2900 J/m³·K
       16 - Eclogite ........... kth=2.900 W/m·k, Cp=3100 J/m³·K
       17 - Gabbro ............. kth=1.900 W/m·k, Cp=2600 J/m³·K
       18 - Gneiss ............. kth=2.900 W/m·k, Cp=2100 J/m³·K
       19 - Granite ............ kth=3.400 W/m·k, Cp=2400 J/m³·K
       20 - Granodiorite ....... kth=3.300 W/m·k, Cp=2600 J/m³·K
       21 - Gravel, dry ........ kth=0.400 W/m·k, Cp=1500 J/m³·K
       22 - Gravel, saturated .. kth=1.800 W/m·k, Cp=2400 J/m³·K
       23 - Gypsum ............. kth=1.600 W/m·k, Cp=2000 J/m³·K
       24 - Gyttja (sapropel) .. kth=0.600 W/m·k, Cp=2000 J/m³·K
       25 - Ice ................ kth=2.200 W/m·k, Cp=1900 J/m³·K
       26 - Ignimbrite ......... kth=3.000 W/m·k, Cp=2100 J/m³·K
       27 - Lamprophyre ........ kth=2.600 W/m·k, Cp=2100 J/m³·K
       28 - Limestone, massive . kth=2.800 W/m·k, Cp=2300 J/m³·K
       29 - Limestone, marly ... kth=2.200 W/m·k, Cp=2300 J/m³·K
       30 - Limestone, oolitic . kth=2.400 W/m·k, Cp=2300 J/m³·K
       31 - Marble ............. kth=2.600 W/m·k, Cp=2000 J/m³·K
       32 - Marl ............... kth=2.100 W/m·k, Cp=2300 J/m³·K
       33 - Marl, clayey ....... kth=2.000 W/m·k, Cp=2200 J/m³·K
       34 - Marl, dolomitic .... kth=2.200 W/m·k, Cp=2300 J/m³·K
       35 - Metaquartzite ...... kth=5.800 W/m·k, Cp=2100 J/m³·K
       36 - Micashist .......... kth=2.000 W/m·k, Cp=2200 J/m³·K
       37 - Peat ............... kth=0.400 W/m·k, Cp=2200 J/m³·K
       38 - Pegmatite .......... kth=3.000 W/m·k, Cp=2100 J/m³·K
       39 - Peridotite ......... kth=4.000 W/m·k, Cp=2700 J/m³·K
       40 - Quartzite .......... kth=6.000 W/m·k, Cp=2100 J/m³·K
       41 - Rhyolite ........... kth=3.300 W/m·k, Cp=2100 J/m³·K
       42 - Salt ............... kth=5.400 W/m·k, Cp=1200 J/m³·K
       43 - Sand, dry .......... kth=0.400 W/m·k, Cp=1400 J/m³·K
       44 - Sand, dry, compacted kth=1.200 W/m·k, Cp=1700 J/m³·K
       45 - Sand, moist ........ kth=1.000 W/m·k, Cp=1800 J/m³·K
       46 - Sand, saturated .... kth=2.400 W/m·k, Cp=2500 J/m³·K
       47 - Sand, frozen ....... kth=2.000 W/m·k, Cp=1500 J/m³·K
       48 - Sandstone .......... kth=2.300 W/m·k, Cp=2000 J/m³·K
       49 - Serpentinite ....... kth=3.000 W/m·k, Cp=2200 J/m³·K
       50 - Shale .............. kth=2.100 W/m·k, Cp=2300 J/m³·K
       51 - Silt, dry .......... kth=0.400 W/m·k, Cp=1600 J/m³·K
       52 - Silt, moist ........ kth=1.800 W/m·k, Cp=2200 J/m³·K
       53 - Siltstone .......... kth=2.400 W/m·k, Cp=2300 J/m³·K
       54 - Syenite ............ kth=2.600 W/m·k, Cp=2400 J/m³·K
       55 - Till (moraine) ..... kth=2.000 W/m·k, Cp=2100 J/m³·K
       56 - Tonalite ........... kth=2.700 W/m·k, Cp=2400 J/m³·K
       57 - Trachyte ........... kth=2.800 W/m·k, Cp=2100 J/m³·K
       58 - Tuff ............... kth=1.100 W/m·k, Cp=1100 J/m³·K
       59 - Water .............. kth=0.600 W/m·k, Cp=4200 J/m³·K
        
Initialize :class:`~pygld.GroundMaterial` from the predefined material at index
#14 and print the results::

    >>> ground = GroundMaterial.ini_from(14)
    >>> print(ground)
    Ground material: Dolomite
    Thermal conductivity: 3.200 W/m·k
    Volumetric heat capacity: 2500 J/m³·K
    Thermal diffusivity: 1.28e-03 m²/s
    
Change the thermal conductivity for another value and print the results::

    >>> ground.kth = 3.35
    >>> print(ground)
    Ground material: User Defined
    Thermal conductivity: 3.350 W/m·k
    Volumetric heat capacity: 2500 J/m³·K
    Thermal diffusivity: 1.34e-03 m²/s

Note that the ground material is now labeled as 'User defined'.

It is also possible to initialize an empty :class:`~pygld.GroundMaterial` and
set the property values manually::

    >>> ground = GroundMaterial()
    >>> ground.kth = 2.27
    >>> ground.Cp = 3278
    >>> print(ground)
    Ground material: User Defined
    Thermal conductivity: 2.270 W/m·k
    Volumetric heat capacity: 3278 J/m³·K
    Thermal diffusivity: 6.92e-04 m²/s