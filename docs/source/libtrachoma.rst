libtrachoma.c
=============

The dynamical core of the NTDMC trachoma model is implemented in the C
programming language. This C core consists of a function
:c:func:`step` and a few :ref:`helper functions <c-helper-functions>`
that Python code can call to initialise the value of various
variables.

Stepping function
-----------------

.. c:autofunction:: step
   :file: trachoma.c

.. _c-helper-functions:

Helper functions
----------------

.. c:autofunction:: set_base_periods
   :file: trachoma.c

.. c:autofunction:: set_infection_parameters
   :file: trachoma.c

.. c:autofunction:: set_background_mortality
   :file: trachoma.c

Structures
-----------

.. c:autostruct:: state
   :file: trachoma.h
   :members:

.. c:autostruct:: output
   :file: trachoma.h
   :members:

