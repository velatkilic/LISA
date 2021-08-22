import lisa
import numpy as np

import unittest

class TestLisa(unittest.TestCase):

	def test_Lisa_default_constructor(self):
		augmentor = lisa.Lisa()
		self.assertNotEqual(augmentor, None)
		print("Test Lisa default constructor passed")

	def test_Lisa_constructor(self):
		lidar = lisa.Lidar()
		water = lisa.Water()
		rain  = lisa.MarshallPalmerRain()
		augmentor = lisa.Lisa(lidar, water, rain)
		self.assertNotEqual(augmentor, None)
		print("Test Lisa constructor passed")

	def test_Lisa_rmin(self):
		augmentor = lisa.Lisa() 
		pc = [[.1,.1,.1,0.5]]
		pcnew = augmentor.augment(pc)
		self.assertEqual(pcnew,[[0,0,0,0,0]])
		print("Test Lisa with range smaller than rmin passed")

	def test_Lisa_rmin_np(self):
		augmentor = lisa.Lisa() 
		pc = np.array([[.1,.1,.1,0.5]])
		pcnew = augmentor.augment(pc)
		self.assertEqual(pcnew,[[0,0,0,0,0]])
		print("Test Lisa with range smaller than rmin passed (numpy array input)")

if __name__=="__main__":
	unittest.main()
