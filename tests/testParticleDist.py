import pylisa
import numpy as np

import unittest

class TestParticleDist(unittest.TestCase):
	def test_MP_default_constructor(self):
		mp = pylisa.MarshallPalmerRain()
		self.assertNotEqual(mp,None)
		print("Test Marshall Palmer rain model default constructor passed")


	def test_MP_Nmodel(self):
		rain = pylisa.MarshallPalmerRain()
		Nd = rain.N_model(0.01, 30)
		self.assertEqual(Nd, 7841.0256295643185)
		print("Test Marshall Palmer rain N model passed")

	def test_MP_Ntotal(self):
		rain = pylisa.MarshallPalmerRain()
		Nt = rain.N_total(30, 0.050e-3)
		self.assertEqual(Nt, 3985.2723697019424)
		print("Test Marshall Palmer rain N total passed")

	def test_MP_Nsample(self):
		rain = pylisa.MarshallPalmerRain()
		dsam = rain.N_sample(30, 0.050e-3, 1000)
		avg  = np.mean(dsam)
		self.assertAlmostEqual(avg, 0.5, 1)
		print("Test Marshall Palmer rain Nsample passed")

	def test_ParticleDist_default_constructor(self):
		pd = pylisa.ParticleDist()
		self.assertNotEqual(pd,None)
		print("Test ParticleDist default constructor passed")

if __name__=="__main__":
	unittest.main()
