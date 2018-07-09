HeatPump class
==================

.. currentmodule:: pygld

.. autoclass:: pygld.HeatPump
   :no-members:
   :no-undoc-members:

Model and Data
---------------------------------

.. autoattribute:: pygld.HeatPump.hpdata
.. autoattribute:: pygld.HeatPump.model
.. automethod:: pygld.HeatPump.set_model

Heat carrier fluid properties
---------------------------------

.. autoattribute:: pygld.HeatPump.fluid
.. autoattribute:: pygld.HeatPump.fr

Heatpump independent properties
---------------------------------

.. autoattribute:: pygld.HeatPump.qbat
.. autoattribute:: pygld.HeatPump.TinHP
.. autoattribute:: pygld.HeatPump.Vf 

.. note::
    A numpy array will always be returned when getting
    :attr:`~pygld.HeatPump.qbat`, :attr:`~pygld.HeatPump.TinHP`, or
    :attr:`~pygld.HeatPump.Vf` independently of the format used to
    set the attribute.

Heatpump dependent properties
---------------------------------

.. autoattribute:: pygld.HeatPump.Tm
.. autoattribute:: pygld.HeatPump.ToutHP
.. autoattribute:: pygld.HeatPump.CAP
.. autoattribute:: pygld.HeatPump.COP


Utilities
---------------------------------

.. automethod:: pygld.HeatPump.in_table
.. automethod:: pygld.HeatPump.get_flowRange
.. automethod:: pygld.HeatPump.get_avail_fluid_types
.. automethod:: pygld.HeatPump.get_avail_heatpump_models
.. automethod:: pygld.HeatPump.print_avail_heatpump_models
.. automethod:: pygld.HeatPump.plot_heatpump_model_goodness
