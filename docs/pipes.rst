Pipe
==================

.. currentmodule:: pygld

.. autoclass:: pygld.Pipe
   :no-members:
   :no-undoc-members:

Geometry
---------------------------------

.. autoattribute:: pygld.Pipe.di
.. autoattribute:: pygld.Pipe.do
.. autoattribute:: pygld.Pipe.wt
.. autoattribute:: pygld.Pipe.Ai
.. autoattribute:: pygld.Pipe.Ao

Material and thermal resistance
--------------------------------------

.. autoattribute:: pygld.Pipe.material
.. automethod:: pygld.Pipe.set_material
.. autoattribute:: pygld.Pipe.Rcond

References
---------------------------------

Pipe conductive thermal resistance :

Hellström, G. 1991. Ground Heat Storage - Thermal Analyses of Duct
    Storage Systems. Ph.D. thesis. University of Lund. Department
    of mathematical physics, Lund, Sweden. Eq.8.2, p.75.

Example
---------------------------------

Import and instantiate the :class:`~pygld.HeatPump` class::

    >>> from pygld import Pipe
    >>> pipe = Pipe(di=3.39852, do=4.2164)
    >>> print(pipe)
    Inner diameter: 3.40 cm
    Outer diameter: 4.22 cm
    Wall thickness: 0.41 cm
    Inner area: 9.07 cm²
    Outer area: 13.96 cm²
    Pipe material: HDPE
    Thermal conductivity: 0.400 W/m·k
    Volumetric heat capacity: 1500 J/m³·K
    Thermal diffusivity: 2.67e-04 m²/s
    Conductive thermal resistance: 0.08580 m·K/W
    
Change and print the thermal properties of the pipe. This must be done through
the pipe's :attr:`.material` as::

    >>> pipe.material.kth = 0.43
    >>> pipe.material.Cp = 0.1540
    >>> print(pipe.material)
    Pipe material: User Defined
    Thermal conductivity: 0.430 W/m·k
    Volumetric heat capacity: 1540 J/m³·K
    Thermal diffusivity: 2.79e-04 m²/s
