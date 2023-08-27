import numpy as np
import ctypes
from ctypes import POINTER, c_ubyte, c_int

import pdb

class Out_c(ctypes.Structure):
    _fields_ = [
        ("inf", ctypes.POINTER(ctypes.c_ubyte)),
        ("dis", ctypes.POINTER(ctypes.c_ubyte)),
        ("lat", ctypes.POINTER(ctypes.c_ubyte)),
        ("ages", ctypes.POINTER(ctypes.c_int)),
        ("nrecords", ctypes.POINTER(ctypes.c_int)),
    ]


class Output:
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
        return Out_c(
            self._inf[offset:].ctypes.data_as(POINTER(c_ubyte)),
            self._dis[offset:].ctypes.data_as(POINTER(c_ubyte)),
            self._lat[offset:].ctypes.data_as(POINTER(c_ubyte)),
            self.ages[offset_ages:].ctypes.data_as(POINTER(c_int)),
            self.nrecordsp,
        )

    def write(self):
        nrecords = len(self.ages) // self.popsize
        for r in range(nrecords):
            start = r * (self.popsize // 8)
            end = (r + 1) * (self.popsize // 8)
            lat = np.unpackbits(self._lat[start:end])
            inf = np.unpackbits(self._inf[start:end])
            dis = np.unpackbits(self._dis[start:end])
            ages = self.ages[r * self.popsize:(r + 1) * self.popsize]

            with open("out.txt", "a") as f:
                np.savetxt(f, np.stack((lat, inf, dis, ages)), fmt="%d")
