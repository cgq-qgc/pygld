HeatCarrierFluid class
======================

.. currentmodule:: pygld

.. autoclass:: pygld.HeatCarrierFluid
   :no-members:
   :no-undoc-members:
   
Fluid type and data
---------------------------------

.. autoattribute:: pygld.HeatCarrierFluid.hcfdata
.. autoattribute:: pygld.HeatCarrierFluid.fluid
.. automethod:: pygld.HeatCarrierFluid.set_fluid

Fluid independent properties
---------------------------------

.. autoattribute:: pygld.HeatCarrierFluid.Tref
.. autoattribute:: pygld.HeatCarrierFluid.fr

Fluid freezing and boiling points
---------------------------------

.. autoattribute:: pygld.HeatCarrierFluid.Tfp
.. autoattribute:: pygld.HeatCarrierFluid.Tbp

Fluid dependent properties
----------------------------------

.. autoattribute:: pygld.HeatCarrierFluid.rho
.. autoattribute:: pygld.HeatCarrierFluid.mu
.. autoattribute:: pygld.HeatCarrierFluid.kth
.. autoattribute:: pygld.HeatCarrierFluid.cp

Fluid derived dependent properties
----------------------------------

.. autoattribute:: pygld.HeatCarrierFluid.nu
.. autoattribute:: pygld.HeatCarrierFluid.Pr
.. autoattribute:: pygld.HeatCarrierFluid.al
.. autoattribute:: pygld.HeatCarrierFluid.Cp

Utilities
---------------------------------

.. automethod:: pygld.HeatCarrierFluid.get_avail_fluid_types

Example
---------------------------------

Import and instantiate the :class:`~pygld.HeatCarrierFluid` class and print
a summary of the fluid's default properties::

    >>> from pygld import HeatCarrierFluid
    >>> hcfluid = HeatCarrierFluid()
    >>> print(hcfluid)
    Type of fluid: water
    Antifreeze volumetric fraction: 0.00
    Freezing point temperature (°C): 0.0
    Boiling point temperature (°C): 100.0
    Temperature of reference (°C): 20.0
    Fluid density in (kg/m³): 998.21
    Cinematic viscosity (Pa·s): 1.00e-03
    Thermal conductivity (W/m·k): 0.599
    Specific heat capacity (J/kg·K): 4185.0
    Kynematic viscosity (m²/s): 1.00e-06
    Prantl number: 7.0
    Thermal diffusivity (m²/s): 1.44e-07
    Volumetric Heat Capacity (J/m³·K): 4.18e+06
    
Print the list of available heat carrier fluids::

    >>> print(hcfluid.get_avail_fluid_types())
    ['prop_glycol', 'ethyl_glycol', 'water']
    
Change the :attr:`~pygld.HeatCarrierFluid.fluid` type,
:attr:`~pygld.HeatCarrierFluid.fr`, and :attr:`~pygld.HeatCarrierFluid.Tref`
values and print a summary of the fluid's updated properties::

    >>> hcfluid.set_fluid('prop_glycol')
    >>> hcfluid.fr = 0.3
    >>> hcfluid.Tref = [28, 14, 0]
    >>> print(hcfluid)
    Type of fluid: prop_glycol
    Antifreeze volumetric fraction: 0.30
    Freezing point temperature (°C): -13.1
    Boiling point temperature (°C): 102.2
    Temperature of reference (°C): [28.0, 14.0, 0.0]
    Fluid density in (kg/m³): [1024.65, 1030.92, 1036.24]
    Cinematic viscosity (Pa·s): [2.33e-03, 3.85e-03, 7.07e-03]
    Thermal conductivity (W/m·k): [0.439, 0.425, 0.409]
    Specific heat capacity (J/kg·K): [3869.9, 3831.1, 3793.0]
    Kynematic viscosity (m²/s): [2.28e-06, 3.73e-06, 6.82e-06]
    Prantl number: [20.6, 34.7, 65.6]
    Thermal diffusivity (m²/s): [1.11e-07, 1.08e-07, 1.04e-07]
    Volumetric Heat Capacity (J/m³·K): [3.97e+06, 3.95e+06, 3.93e+06]


References
---------------------------------

Properties for pure water were taken from p.154 of:

VDI-Gesellschaft, 2010. VDI Heat Atlas. VDI-Verlag GmbH, Dusseldorf,
    Germany, Berlin; London.

Properties for propylene and ethlyne glycol were taken from:

ASHRAE, 2009. 2009 ASHRAE Handbook - Fundamentals. Si edition ed.,
American Society of Heating, Refrigerating and Air-Conditioning
Engineers, Atlanta, GA. chap.31, pp. 823-835.
