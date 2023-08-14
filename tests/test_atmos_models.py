import unittest
import os
from pathlib import Path
from pylisa.mie_wrapper import calc_qs
from pylisa.atmos_models import alpha_beta, n_mp_rain

import numpy as np
import matplotlib.pyplot as plt

class TestAtmosModels(unittest.TestCase):
    def test_compare_alpha_beta_to_approx_theory(self):
        m = 1.328
        lam = 905. # nm
        d = np.logspace(-5, 1, 2000) # diameter range for Mie, (mm)
        x = 1e6 * np.pi * 0.5 * d / lam # 1e6 since D in mm and lam in nm
        qs = calc_qs(x, m, 1.)
        q_ext = qs[:,0]
        q_back = qs[:,2]
        N = 100
        rain_rate = np.linspace(1,20,N)

        # From: Ulbrich and Atlas, "Extinction of visible and
        # Infrared Radiation in rain...", eq 8
        k1 = 1.45*(rain_rate**.64)
        k2 = np.zeros(k1.shape)

        alpha = np.zeros(k1.shape)
        beta = np.zeros(k1.shape)

        for i in range(N):
            nd_rain = n_mp_rain(d, rain_rate[i])
            
            alpha[i], beta[i] = alpha_beta(nd_rain, d, q_ext, q_back)
            
            k2[i] = alpha[i] * 1000 * (10/np.log(10))
        
        plt.figure()
        plt.plot(rain_rate, k1, lw=2)
        plt.plot(rain_rate, k2, lw=2)
        plt.ylabel('Extinction Coefficient, dB/km')
        plt.xlabel('Rain rate, mm/hr')
        plt.legend(('Asymptotic Model','Mie Model'))
        plt.savefig(os.path.join(Path(__file__).resolve().parent, "extinction_mie_vs_asymp.png"))
        plt.close()

if __name__ == "__main__":
    unittest.main()