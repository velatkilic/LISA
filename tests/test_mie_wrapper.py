import unittest
import os
from pathlib import Path
from pylisa.mie_wrapper import bh_mie, calc_qs

import numpy as np
import matplotlib.pyplot as plt

class TestMieWrapperBHAppendixA(unittest.TestCase):
    def setUp(self):
        PI = 3.14159265
        n_med = 1.
        n_ref = 1.55
        sphere_radius = 0.525
        wavelength = .6328

        self.size_parameter = 2 * PI * sphere_radius * n_med / wavelength
        self.q_ext, self.q_sca, self.q_back, self.g_aniso, self.s1, self.s2 = bh_mie(self.size_parameter, n_ref, n_med, n_ang=11)

    def test_efficiencies_q_are_float(self):
        self.assertIsInstance(self.q_ext, float)
        self.assertIsInstance(self.q_sca, float)
        self.assertIsInstance(self.q_back, float)
    
    def test_s_params_are_complex(self):
        for s in self.s1:
            self.assertIsInstance(s, complex)
        
        for s in self.s2:
            self.assertIsInstance(s, complex)
    
    def test_size_parameter(self):
        self.assertAlmostEqual(self.size_parameter, 5.213, 3)
    
    def test_q_sca_has_correct_value(self):
        self.assertAlmostEqual(self.q_sca, 3.1054, 4)
    
    def test_q_ext_has_correct_value(self):
        self.assertAlmostEqual(self.q_ext, 3.1054, 4)
    
    def test_q_back_has_correct_value(self):
        self.assertAlmostEqual(self.q_back, 2.9253, 4)
    
    def test_s_and_p_params(self):
        # TODO: transfer over the table from Appendix A
        pass

class TestMieWrapperImagIndex(unittest.TestCase):
    def setUp(self):
        PI = 3.14159265
        n_med = 1.
        n_ref = 1.55 + 1j
        sphere_radius = 0.525
        wavelength = .6328

        size_parameter = 2 * PI * sphere_radius * n_med / wavelength
        self.q_ext, self.q_sca, self.q_back, self.g_aniso, self.s1, self.s2 = bh_mie(size_parameter, n_ref, n_med, n_ang=5)
    
    def test_q_sca_has_correct_value(self):
        self.assertAlmostEqual(self.q_sca, 1.3413, 2)
    
    def test_q_ext_has_correct_value(self):
        self.assertAlmostEqual(self.q_ext, 2.5895, 2)
    
    def test_q_back_has_correct_value(self):
        self.assertAlmostEqual(self.q_back, 0.14703, 2)

class TestMieWrapperXArray(unittest.TestCase):
    def test_generate_plot(self):
        nref = 1.55
        nmed = 1.
        x_arr = np.linspace(0.1,30,100)
        qs = calc_qs(x_arr, nref, nmed)
        q_sca_mie = qs[:,1]
        xi = (nref**2 - 1) / (nref**2 + 2)
        q_sca_rayleigh = (8./3.)* x_arr**4 * xi**2
        q_sca_geom = np.ones((len(x_arr, ))) * 2. # see the "extinction paradox"
        plt.figure()
        plt.loglog(x_arr, q_sca_mie, label="Mie")
        plt.loglog(x_arr, q_sca_rayleigh, label="Rayleigh")
        plt.loglog(x_arr, q_sca_geom, label = "Geometric")
        plt.ylim([1e-5, 10])
        plt.xlabel("Size parameter, ka")
        plt.ylabel("Scattering Efficiency, Qsca")
        plt.legend()
        plt.savefig(os.path.join(Path(__file__).resolve().parent, "mie_plot.png"))
        plt.close()

if __name__ == "__main__":
    unittest.main()