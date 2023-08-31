import ctypes
from ctypes import POINTER, c_int, c_double, c_ubyte
import numpy as np

# FIXME: Implement initial infection in this module
from .init import infected


class Population:
    """Represent a population of individuals.

    A Population object keeps track of population properties in a
    Structure of Arrays (SoA) fashion.  It maintains a collection of
    one-dimensional arrays of length equal to the number of
    individuals represented in the population. These arrays are:

    - The **latent** state of individuals (`lat`, boolean)
    - The **infected** state of individuals (`inf`, boolean)
    - The **diseased** state of individuals (`dis`, boolean)
    - The clock value for each individual (`clock`, integer)
    - The age value for each individual (`ages`, 8-bits unsigned integer)
    - The infection count for each individual (`count`, unsigned integer)
    - The bacterial load for each individual (`bact_load`, float)

    Initially, the arrays listed above are sorted according to
    individuals' increasing age value.  However, the a ``Population``
    instance does not garantuee that this order is maitained whenever
    one or more of these arrays are modified: is it up to users of the
    ``Population`` class to maintain this sorted order if necessary.

    ..note:

    The latent, infected and diseased state of individual is
    internally represented as an array of (N / 8) 8-bits unsigned
    integers (np.uint8), with one bit per individual. This is hidden
    to users of the ``Population`` class by providing access to these
    arrays through properties, packing or unpacking bits whenever an
    infection state array is accessed or set.

    :param ages: List of ages for each individuals in the population
    :type ages: array_like
    :param latent_base: The base latent period for each individuals in
        the population
    :type latent_base: array_like
    :param rng: A random generator instance
    :type rng: numpy.random.Generator
    :param lib: A handle to the C core library
    :type lib: ctypes.CDLL

    Examples:
    ~~~~~~~~~

    Instanciating a ``Population`` object:
    
    >>> import numpy, ctypes
    >>> from ntdmc_trachoma.state import Population
    >>> ages = [832, 156, 2340, 1664, 3328]
    >>> latent_base = [2] * 5
    >>> rng = numpy.random.default_rng()
    >>> lib = ctypes.CDLL('./libtrachoma.so')
    >>> pop = Population(ages, latent_base, rng, lib)


    Accessing ``Population`` properties:
    
    >>> pop.ages
    array([156, 832, 1664, 2340, 3328])
    >>> pop.lat
    array([True, False, False, False, False])
    >>> pop.inf
    array([True, False, False, True, False])
    >>> pop.dis
    array([False, True, False, True, False])
    """
    def __init__(self, ages, latent_base, rng, lib):

        self.size = len(ages)
        self.ages = np.sort(ages).astype(np.int32)
        self.clock = np.zeros(self.size, dtype=np.int32) - 1
        self.count = np.zeros(self.size, dtype=np.int32)
        self.bact_load = np.zeros(self.size, dtype=np.float64)

        #FIXME: I think a Population should be initialised by passing
        #ages and infected status only. This would avoid having to
        # pass the random generator reference.
        latent = infected(self.size, 0.3, rng)
        self._lat = np.packbits(latent)
        self._inf = self._lat.copy()
        self._dis = np.packbits(np.zeros(self.size, dtype=np.bool_))

        # FIXME: Population init shouldn't have to go down to C
        # just to init the latent time of initially infected people
        lib.setlatenttime.restype = ctypes.c_int
        lib.setlatenttime.argtypes = [ctypes.c_int] * 3
        self.clock[latent] = [
            lib.setlatenttime(base, 0, age)
            for age, base in zip(self.ages[latent], latent_base[latent])
        ]
        self.count[latent] = 1
        ninfected = np.sum(self.count)
        lib.get_load.restype = ctypes.c_double
        lib.get_load.argtypes = [ctypes.c_int]
        self.bact_load[latent] = [
            lib.get_load(1) for _ in range(ninfected)
        ]

    @property
    def _as_parameter_(self):
        # Tell ctypes what to actually pass as an argument when a
        # Population object is passed as a parameter to a foreign
        # function.
        return Pop_c(
            self.size,
            self._inf.ctypes.data_as(POINTER(c_ubyte)),
            self._dis.ctypes.data_as(POINTER(c_ubyte)),
            self._lat.ctypes.data_as(POINTER(c_ubyte)),
            self.clock.ctypes.data_as(POINTER(c_int)),
            self.ages.ctypes.data_as(POINTER(c_int)),
            self.count.ctypes.data_as(POINTER(c_int)),
            self.bact_load.ctypes.data_as(POINTER(c_double)),
        )

    @property
    def inf(self):
        return np.unpackbits(self._inf).astype(np.bool_)

    @inf.setter
    def inf(self, a):
        self._inf = np.packbits(a)

    @property
    def dis(self):
        return np.unpackbits(self._dis).astype(np.bool_)

    @dis.setter
    def dis(self, a):
        self._dis = np.packbits(a)

    @property
    def lat(self):
        return np.unpackbits(self._lat).astype(np.bool_)

    @lat.setter
    def lat(self, a):
        self._lat = np.packbits(a)


class Pop_c(ctypes.Structure):
    _fields_ = [
        ("n", c_int),
        ("inf", ctypes.POINTER(ctypes.c_ubyte)),
        ("dis", ctypes.POINTER(ctypes.c_ubyte)),
        ("lat", ctypes.POINTER(ctypes.c_ubyte)),
        ("clockm", ctypes.POINTER(ctypes.c_int)),
        ("ages", ctypes.POINTER(ctypes.c_int)),
        ("count", ctypes.POINTER(ctypes.c_int)),
        ("bactload", ctypes.POINTER(ctypes.c_double)),
    ]
