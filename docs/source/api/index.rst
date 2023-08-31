API reference
=============

.. toctree::
   :maxdepth: 1
   :caption: Contents:

Simulation
----------

.. autoclass:: ntdmc_trachoma.simulation.Simulation
   :members:

Population
----------

.. autoclass:: ntdmc_trachoma.state.Population
   :members:

Events
------

.. automodule:: ntdmc_trachoma.mda

Available event classes
-----------------------

.. autoclass:: ntdmc_trachoma.mda.MDA
.. autofunction:: ntdmc_trachoma.mda.start_sim
.. autofunction:: ntdmc_trachoma.mda.end_sim

Creating the events list from a scenario definition
---------------------------------------------------

The ``process_scenario_definition.py`` module exposes only one
function `process_scenario_definition`.  This function is used to
process a scenario definition file and construct a list of events a
:py:class:`Simulation <ntdmc_trachoma.simulation.Simulation>` instance
can work with. See :py:func:`Simulation.simulate
<ntdmc_trachoma.simulation.Simulation.simulate>`

.. autofunction:: ntdmc_trachoma.process_scenario_definition.process_scenario_definition
   
