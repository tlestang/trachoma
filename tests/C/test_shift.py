import ctypes
from importlib import util
from ctypes import c_int
import numpy as np
import numpy.testing as nptest


LIBTRACHO_PATH = util.find_spec(
    "ntdmc_trachoma.libtrachoma"
).origin
lib = ctypes.CDLL(LIBTRACHO_PATH)


def test_shift():
    a = np.array(
        [
            1, 0, 1, 1, 1, 0, 0, 1,
            0, 0, 0, 1, 0, 0, 1, 1,
            0, 1, 1, 0, 0, 1, 0, 0,
            1, 0, 0, 0, 1, 0, 1, 1,
        ]
    )
    expected = np.array(
        [
            0, 0, 0, 1, 0, 1, 1, 1,
            0, 0, 1, 0, 0, 0, 1, 0,
            0, 1, 1, 0, 1, 1, 0, 0,
            1, 0, 0, 1, 0, 0, 0, 1,
        ]
    )
    a_packed = np.packbits(a)  # [185, 19, 100, 139]
    expected_packed = np.packbits(expected)

    shift = ctypes.c_int(3)
    size = ctypes.c_int(a.size // 8)
    lib.shift.restype = None
    lib.shift.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.uint8, ndim=1),
        ctypes.c_int,
        ctypes.c_int,
    ]
    lib.shift(a_packed, shift, size)
    nptest.assert_array_equal(expected_packed, a_packed)


def test_bgd_bitarray():
    a = np.array(
        [
            1, 0, 1, 1, 1, 0, 0, 1,
            0, 0, 0, 1, 0, 0, 1, 1,
        ]
    )
    expected = np.array(
        [
            0, 1, 0, 1, 1, 1, 0, 0,
            1, 0, 0, 0, 0, 0, 1, 1,
        ]
    )
    a_packed = np.packbits(a)  # [185, 19]
    expected_packed = np.packbits(expected)

    idx = ctypes.c_int(11)
    size = ctypes.c_int(a.size // 8)
    lib.bgd_death_bitarray.restype = None
    lib.bgd_death_bitarray.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.uint8, ndim=1),
        ctypes.c_int,
        ctypes.c_int,
    ]
    lib.bgd_death_bitarray(a_packed, idx, size)
    nptest.assert_array_equal(expected_packed, a_packed)


def test_rotate_bitarray():
    a = np.array(
        [
            1, 0, 1, 1, 1, 0, 0, 1,
            0, 0, 0, 1, 0, 0, 1, 1,
            0, 1, 1, 0, 0, 1, 0, 0,
            1, 0, 0, 0, 1, 0, 1, 1,
        ]
    )
    expected = np.array(
        [
            0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 1, 1, 1,
            0, 0, 1, 0, 0, 0, 1, 0,
            0, 1, 1, 0, 1, 1, 0, 0,
        ]
    )
    a_packed = np.packbits(a)  # [185, 19, 100, 139]
    expected_packed = np.packbits(expected)

    n = ctypes.c_int(11)
    size = ctypes.c_int(a.size // 8)
    lib.rotate_bitarray.restype = None
    lib.rotate_bitarray.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.uint8, ndim=1),
        ctypes.c_int,
        ctypes.c_int,
    ]
    lib.rotate_bitarray(a_packed, n, size)
    nptest.assert_array_equal(expected_packed, a_packed)


def test_rotate():
    TenIntegers = ctypes.c_int * 10
    a = TenIntegers(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    lib.rotate.restype = None
    lib.rotate.argtypes = [TenIntegers, ctypes.c_int, ctypes.c_int]
    lib.rotate(a, c_int(3), c_int(10), c_int(72))
    assert (
        [72, 72, 72, 1, 2, 3, 4, 5, 6, 7] ==
        [i for i in a]
    )

