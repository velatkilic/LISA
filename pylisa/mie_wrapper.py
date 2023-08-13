"""
ctypes wrapper for mie shared library
"""
from ctypes import cdll, Structure, c_float, c_int, POINTER, byref
import platform
import os
from pathlib import Path
import numpy as np


os_name = platform.system()
if os_name == "Linux":
    lib_path = os.path.join(Path(__file__).resolve().parent, "mie.so")
elif os_name == "Windows":
    lib_path = os.path.join(Path(__file__).resolve().parent, "mie.dll")
elif os_name == "Darwin":
    lib_path = os.path.join(Path(__file__).resolve().parent, "mie.dylib")
else:
    raise OSError(
        "Your OS is not supported but you can \
                  compile mie.c and change this line. \
                  Check out compilation instructions on \
                  Github: https://github.com/velatkilic/LISA"
    )

lib = cdll.LoadLibrary(lib_path)


class Complex(Structure):
    """
    Wrapper for the Complex struct in mie.c
    """

    _fields_ = [("r", c_float), ("i", c_float)]


cabs = lib.Cabs
cabs.argtypes = [
    Complex,
]
cabs.restype = c_float

cmulti = lib.Cmulti
cmulti.argtypes = [Complex, Complex]
cmulti.restype = Complex

cconj = lib.Cconj
cconj.argtypes = [
    Complex,
]
cconj.restype = Complex

_lib_bh_mie = lib.BHMie
_lib_bh_mie.argtypes = [
    POINTER(c_float),  # x
    POINTER(Complex),  # ref_rel
    POINTER(c_int),    # n_ang
    POINTER(Complex),  # S1
    POINTER(Complex),  # S2
    POINTER(c_float),  # q_ext
    POINTER(c_float),  # q_sca
    POINTER(c_float),  # q_back
    POINTER(c_float),  # g_aniso
]


def bh_mie(x, nref, nmed, n_ang=2):
    """Mie scattering calculation for spheres using Bohren Huffman implementation

    Args:
        x (float): Size parameter
        nref (complex or float): refractive index of the material
        nmed (float): refractive index of the medium
        n_ang (int, optional): Number of angles to calculate. Defaults to 2.

    Returns:
        q_ext (float): extinction efficiency
        q_sca (float): scattering efficiency
        q_back (float): back scattering efficiency
        g_aniso (float): anisotropy
        s1 (List[float]): scattering matrix component (as a function of angle)
        s2 (List[float]): scattering matrix component (as a function of angle)
    """
    n = n_ang * 2 - 1
    x = c_float(x)
    ref_rel = nref / nmed
    ref_rel = Complex(r=ref_rel.real, i=ref_rel.imag)
    n_ang = c_int(n_ang)
    s1 = (Complex * n)()
    s2 = (Complex * n)()
    q_ext = c_float(0)
    q_sca = c_float(0)
    q_back = c_float(0)
    g_aniso = c_float(0)
    _lib_bh_mie(
        byref(x),
        byref(ref_rel),
        byref(n_ang),
        s1,
        s2,
        byref(q_ext),
        byref(q_sca),
        byref(q_back),
        byref(g_aniso),
    )
    s1 = [num.r + 1j * num.i for num in s1]
    s2 = [num.r + 1j * num.i for num in s2]
    return q_ext.value, q_sca.value, q_back.value, g_aniso.value, s1, s2


def calc_qs(x_arr, nref, nmed):
    """_summary_

    Args:
        x_arr (np.array[float]): Size parameter array of size N
        nref (complex or float): refractive index of the sphere
        nmed (float): refractive index of the medium

    Returns:
        qs (np.array[float]): {extinction, scattering and back scattering} efficiencies (N,3)
    """
    qs = np.zeros((len(x_arr), 3), dtype=np.float32)
    n_ang = 2
    n = n_ang * 2 - 1
    ref_rel = nref / nmed
    ref_rel = Complex(r=ref_rel.real, i=ref_rel.imag)
    n_ang = c_int(n_ang)
    s1 = (Complex * n)()
    s2 = (Complex * n)()
    q_ext = c_float(0)
    q_sca = c_float(0)
    q_back = c_float(0)
    g_aniso = c_float(0)
    for i in range(len(x_arr)):
        x = c_float(x_arr[i])
        _lib_bh_mie(
            byref(x),
            byref(ref_rel),
            byref(n_ang),
            s1,
            s2,
            byref(q_ext),
            byref(q_sca),
            byref(q_back),
            byref(g_aniso),
        )
        qs[i, 0] = q_ext.value
        qs[i, 1] = q_sca.value
        qs[i, 2] = q_back.value
    return qs


if __name__ == "__main__":
    PI = 3.14159265
    # Reproduce Bohren & Huffman Appendix A data
    ref_med = 1.0
    ref_re = 1.55
    ref_im = 0.0
    sphere_radius = 0.525
    wavelength = 0.6328
    size_parameter = 2 * PI * sphere_radius * ref_med / wavelength
    print(f"Size parameter: {size_parameter: .3f}")

    n_ang = 11
    n = n_ang * 2 - 1

    x = c_float(size_parameter)
    ref_rel = Complex(r=ref_re / ref_med, i=ref_im / ref_med)
    n_ang = c_int(n_ang)
    s1 = (Complex * n)()
    s2 = (Complex * n)()
    q_ext = c_float(0)
    q_sca = c_float(0)
    q_back = c_float(0)
    g_aniso = c_float(0)
    _lib_bh_mie(
        byref(x),
        byref(ref_rel),
        byref(n_ang),
        s1,
        s2,
        byref(q_ext),
        byref(q_sca),
        byref(q_back),
        byref(g_aniso),
    )

    print(f"{q_sca=} {q_ext=} {q_back=} {g_aniso=}")
    ang = np.linspace(0, 180, n)
    s11nor = 0.5 * (cabs(s1[0]) ** 2 + cabs(s2[0]) ** 2)
    for j in range(n):
        temp1 = cabs(s1[j])
        temp2 = cabs(s2[j])
        s11 = 0.5 * (temp1**2 + temp2**2)
        s12 = 0.5 * (temp1**2 - temp2**2)

        pol = -s12 / s11
        ctemp = cmulti(s2[j], cconj(s1[j]))
        s33 = ctemp.r
        s33 = s33 / s11
        s34 = ctemp.i
        s34 = s34 / s11
        s11 = s11 / s11nor
        print(f"{ang[j]:.2f}, {s11:14.6g}, {pol:14.6g}, {s33:14.6g}, {s34:14.6g}")
