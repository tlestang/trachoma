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
        with open("infection_state.bin", "wb") as f:
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
