Tutorial
========

This short tutorial will guide you through creating and running a
trachoma simulation using the NTDMC trachoma model.  It also
illustrates an example of post-processing simulation data to display
the number of infected individuals in a given age group over time.

Creating a simulation
---------------------

To start off with, make sure you have prepared valid parameters and
scenario files. See :doc:`parameters` and :doc:`scenarios`.
For this tutorial, wee will use the following parameters and scenario:

- :download:`parameters.json`
- :download:`scenario.json`

Next, create a new :py:class:`Simulation
<ntdmc_trachoma.simulation.Simulation>` instance:

.. code-block:: python

   from pathlib import Path
   from ntdmc_trachoma.simulation import Simulation

   param_path = Path("parameters.json")
   sim = Simulation(param_path)

Creating a :py:class:`Simulation
<ntdmc_trachoma.simulation.Simulation>` object automatically creates a
:py:class:`Population <ntdmc_trachoma.state.Population>` in which 30%
of the individuals are infected.  Infected individuals are randomly
(uniformly) selected.

Running a simulation
--------------------

Next, we can run scenario over a range of values of the spread
parameter :math:`\beta`.  Let's say these values are listed in a text
file ``beta_values.txt``, with one value per line.

.. code-block:: python

   betavals = [0.2, 0.35, 0.4]
   scenario_path = Path('scenario.json')
   sim.simulate(scenario_path, betavals, record=True)

Setting ``record=True``, enables recording of the population state at
each iteration of the model within the ``Simulation`` object's
:py:class:`output buffer <ntdmc_trachoma.output.Output>`.

Writing simulation output to disk
---------------------------------

Once the simulation is finished, we can write the content of the
output buffer on disk:

.. code-block:: python

   sim.output.write()

This creates (or overwrites) 4 binary files containing:

- The latent state (``latent_state.bin``).
- The infected state (``infected_state.bin``).
- The diseased state (``diseased_state.bin``).
- The individuals' age (``ages.bin``).

over the course of the simulation.  For instance, if you simulated the
model over 52 weeks (52 iterations) with a population size of 1000,
``ages.bin`` would contain 52000 integers. The first 1000 corresping
to the age distribution for the first iteration of the model, the next
1000 after corrsponding to the second iteration, and so on and so forth.

.. warning::

   Because the output files are binary files, the type of the data is
   important.  Infection state data is written as 8-bits unsigned
   integers (``numpy.uint8``) with one single bit per individual.  Age
   data is also written as 8-bits unsigned integers, but this time with
   one word (8-bits) per individual.  You check that the size of
   ``ages.bin`` is 8 times the size of any of the other 3 output files.

Post-processing simulation output: an example
---------------------------------------------

We can now use this model output for post-processing.  As an example,
let's plot the evolution of trachoma prevalence among 9 to 15 years
olds.

The plan is straightforward: for each value of the math:`\beta`
parameter, we read the corresponding records in the ``ages.bin`` and
``infection_state.bin`` ouput files.  Because individuals are sorted
by increasing age, we can easily determine which subset of the
population correspond to individuals ages between 9 and 15 years old.

We start by looping over the values of :math:`\beta`:

.. code-block:: python

   popsize = 1024
   nrecords = 1144
   nbytes = (nrecords * popsize) // 8
   for ibeta, beta in enumerate(betavals):
       ages = np.fromfile(
            "ages.bin",
            dtype=np.uint8,
            count=size * nrecords,
            offset=ibeta * size * nrecords,
       )
       inf = np.fromfile(
           "infection_state.bin",
           dtype=np.uint8,
           count=count,
           offset=ibeta * count,
       )


Because our :download:`scenario <scenario.json>` spans 1144 weeks
and records are made every weeks, we set ``nrecords=1144``.  The
infected state of individuals are packed together into groups of 8
bits, with one bit per individual.  To read all the age records
made for a given value of :math:`\beta`, we need to read
``nrecords * popsize`` integers. Because age data is represented by
1-byte unsigned integers, that's ``nrecords * popsize`` bytes to
read.  Similarly for the amount of memory to read all the infected
state records, expected with divide the bytes count by 8. This is
because the infected states of individual are packed into 1-byte
(8-bits) integers with one bit per individual.

Speaking of, let's unpack these bits into a boolean array of size
``popsize``.

.. code-block:: python

   rec_size = popsize // 8
   inf_records = [
       # Unpack infected state record into a boolean array
       np.unpackbits(
           inf_records_packed[i * rec_size:(i + 1) * rec_size]
       ).astype(np.bool_)
       for i in range(nrecords)
   ]

The above lines of Python generate a list of NumPy arrays, with
each array containing data for an unpacked record, that is a
list of ``popsize`` True/False values.

Next, we need to identify which of the individuals are aged
between 9 and 15 years old.  We need to do so for each record,
since age distribution varies over the course of a simulation,
and therefore along records.  We then simply count how many
individual in the sub-population are infected.

.. code-block:: python

   ninf = [
       np.count_nonzero(
           inf[(ages >= 9) & (ages <= 15 * 52)]
       )
       for inf, ages in zip(inf_records, ages_records)
   ]

Finally, we can plot the infection count over time:

.. code-block:: python

   plt.plot(ninf, label=f"beta = {beta})

.. image:: /static/tutorial_plot.png

You can download the full post-processing script :download:`here <plot_ninf.py>`.
