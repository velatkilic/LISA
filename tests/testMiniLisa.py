import pylisa
import numpy as np

import unittest

class TestMiniLisa(unittest.TestCase):

	def test_MiniLisa_default_constructor(self):
		augmentor = pylisa.MiniLisa()
		self.assertNotEqual(augmentor, None)
		print("Test MiniLIsa default constructor passed")

	def test_MiniLisa_constructor(self):
		lidar = pylisa.Lidar()
		water = pylisa.Water()
		rain  = pylisa.MarshallPalmerRain()
		augmentor = pylisa.MiniLisa(lidar, water, rain)
		self.assertNotEqual(augmentor, None)
		print("Test MiniLisa constructor passed")

	def test_MiniLisa_rmin(self):
		augmentor = pylisa.MiniLisa() 
		pc = [[.1,.1,.1,0.5]]
		pcnew = augmentor.augment(pc)
		self.assertEqual(pcnew,[[0,0,0,0]])
		print("Test MiniLisa with range smaller than rmin passed")

	def test_MiniLisa_rmin_np(self):
		augmentor = pylisa.MiniLisa() 
		pc = np.array([[.1,.1,.1,0.5]])
		pcnew = augmentor.augment(pc)
		self.assertEqual(pcnew,[[0,0,0,0]])
		print("Test MiniLisa with range smaller than rmin passed (numpy array input)")

if __name__=="__main__":
	unittest.main()
