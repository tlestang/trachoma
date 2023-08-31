import numpy as np
import ctypes
from ctypes import POINTER, c_ubyte, c_int


class Output:
    """Container for records of the population state at a given time.

    An instance of class ``Output`` contains various NumPy arrays
    acting as a in-memory storage area for recording various
    population properties over the course of a model simulation.  Each
    property (e.g. age or infected state) is stored into one
    corresponding long one-dimensionanal array in which records are
    appended to each others.  For instance, if K records are made,
    ``output.age`` will contain ``K * N`` integers, where ``N`` is the
    number of individuals in the population.

    In practice, the size storage arrays is not dynamic: it is set
    upon initialisation of an ``Ouput`` instance object.  The size of
    the array must be enough to hold all the records to be made during
    the object's lifetime.

    :param popsize: Number of individuals in the population
    :type popsize: int
    :param max_records: Maximum number of records the output buffer
      can hold.
    :type max_records: int

    Instances of class ``Output`` are storage buffers that are meant
    to be passed down to the :c:func:`C step function
    <step>`. ``Output`` instances are automatically converted to a
    :c:struct:`struct output <output>` object, of which a reference is
    passed to .  If an non-``NULL`` pointer is passed to
    :c:func:`step`, the records will be at every iteration of the
    model.

    .. seealso::

       :c:struct:`struct output <output>`

    Example
    ~~~~~~~

    .. code-block:: python

       pop = Population(...)
       betavals = [0.2, 0.3, 0.4]
       nweeks = 52

       # Instanciate Output class with the right amount of storage
       # space.
       out = Output(
                popsize=len(self.pop.ages),
                max_records=nsteps * len(betavals),
       )

       # Load C library and step model.  The ``step`` function writes
       # to the output buffer at every iteration of the model.
       clib = ctypes.CDLL("/path/to/trachoma/c/lib.so")
       for beta in betavals:
           clib.step(pop, out, nweeks, beta)

       # Finally, we can write the content of the output buffer on
       # disk.
       out.write()

    """
    def __init__(self, popsize, max_records):
        self.popsize = popsize
        infection_store_size = (popsize // 8) * max_records
        self._lat = np.zeros(infection_store_size, dtype=np.uint8)
        self._inf = np.zeros(infection_store_size, dtype=np.uint8)
        self._dis = np.zeros(infection_store_size, dtype=np.uint8)
        self.ages = np.zeros(popsize * max_records, dtype=np.int32)
        self.nrecords = c_int(0)
        self.nrecordsp = ctypes.pointer(self.nrecords)

    @property
    def _as_parameter_(self):
        offset = (self.popsize // 8) * self.nrecords.value
        offset_ages = self.popsize * self.nrecords.value
        # pdb.set_trace()
        return ctypes.byref(
            Out_c(
                self._inf[offset:].ctypes.data_as(POINTER(c_ubyte)),
                self._dis[offset:].ctypes.data_as(POINTER(c_ubyte)),
                self._lat[offset:].ctypes.data_as(POINTER(c_ubyte)),
                self.ages[offset_ages:].ctypes.data_as(POINTER(c_int)),
                self.nrecordsp,
            )
        )

    def write(self):
        """Write the output buffer to disk"""
        with open("infected_state.bin", "wb") as f:
            f.write(
                self._inf.tobytes()
            )
        with open("diseased_state.bin", "wb") as f:
            f.write(
                self._dis.tobytes()
            )
        with open("latent_state.bin", "wb") as f:
            f.write(
                self._lat.tobytes()
            )
        with open("ages.bin", "wb") as f:
            f.write(self.ages.astype(np.uint8).tobytes())


class Out_c(ctypes.Structure):
    _fields_ = [
        ("inf", ctypes.POINTER(ctypes.c_ubyte)),
        ("dis", ctypes.POINTER(ctypes.c_ubyte)),
        ("lat", ctypes.POINTER(ctypes.c_ubyte)),
        ("ages", ctypes.POINTER(ctypes.c_int)),
        ("nrecords", ctypes.POINTER(ctypes.c_int)),
    ]
