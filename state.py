import ctypes
from ctypes import POINTER, c_int, c_double, c_ubyte
import numpy as np

import periods
from infection import get_load
from init import infected


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


class Population:
    def __init__(self, ages, latent_base, rng):

        self.size = len(ages)
        self.ages = np.sort(ages).astype(np.int32)
        self.clock = np.zeros(self.size, dtype=np.int32) - 1
        self.count = np.zeros(self.size, dtype=np.int32)
        self.bact_load = np.zeros(self.size, dtype=np.float64)

        latent = infected(self.size, 0.3, rng)
        self._latent = np.packbits(latent)
        self._inf = self._latent.copy()
        self._dis = np.packbits(np.zeros(self.size, dtype=np.bool_))

        self.clock[latent] = periods.latent_time(latent_base[latent])
        self.count[latent] = 1
        ninfected = np.sum(self.count)
        self.bact_load[latent] = get_load(np.zeros(ninfected) + 1)

    @property
    def _as_parameter_(self):
        return Pop_c(
            self.size,
            self.inf.ctypes.data_as(POINTER(c_ubyte)),
            self.dis.ctypes.data_as(POINTER(c_ubyte)),
            self.latent.ctypes.data_as(POINTER(c_ubyte)),
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
