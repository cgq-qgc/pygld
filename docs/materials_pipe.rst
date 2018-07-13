Pipe Material
---------------------------------

.. currentmodule:: pygld

.. autoclass:: pygld.PipeMaterial
   :no-members:
   :no-undoc-members:

Thermal Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoattribute:: pygld.PipeMaterial.Cp
.. autoattribute:: pygld.PipeMaterial.kth
.. autoattribute:: pygld.PipeMaterial.al

Utilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: pygld.PipeMaterial.get_predefined_materials
.. automethod:: pygld.PipeMaterial.print_predefined_materials
.. automethod:: pygld.PipeMaterial.init_as

Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Import and print the list of predefined :class:`~pygld.PipeMaterial`::

    >>> from pygld import PipeMaterial
    >>> PipeMaterial.print_predefined_materials()
    Predefined materials for the 'Pipe':
        0 - HDPE ...... kth=0.400 W/m·k, Cp=1500 J/m³·K
        1 - Geoperformx kth=0.700 W/m·k, Cp=1500 J/m³·K
        
Initialize :class:`~pygld.PipeMaterial` from the predefined material at index
#1 and print the results::

    >>> pipe = PipeMaterial.ini_from(1)
    >>> print(pipe)
    Pipe material: Geoperformx
    Thermal conductivity: 0.700 W/m·k
    Volumetric heat capacity: 1500 J/m³·K
    Thermal diffusivity: 4.67e-04 m²/s
    
Change the thermal conductivity for another value and print the results::

    >>> pipe.kth = 1.75
    >>> print(pipe)
    Pipe material: User Defined
    Thermal conductivity: 0.670 W/m·k
    Volumetric heat capacity: 1500 J/m³·K
    Thermal diffusivity: 4.47e-04 m²/s

Note that the grout material is now labeled as 'User defined'.

It is also possible to initialize an empty :class:`~pygld.PipeMaterial` and
set the value of the properties manually::

    >>> pipe = PipeMaterial()
    >>> pipe.kth = 2
    >>> pipe.Cp = 3500
    >>> print(pipe)
    Pipe material: User Defined
    Thermal conductivity: 0.530 W/m·k
    Volumetric heat capacity: 1767 J/m³·K
    Thermal diffusivity: 3.00e-04 m²/s
