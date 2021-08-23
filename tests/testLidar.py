import pylisa
import numpy as np

import unittest

class TestLidar(unittest.TestCase):
	def test_lidar_default_constructor(self):
		lidar = pylisa.Lidar()
		self.assertNotEqual(lidar, None)
		print("Test Lidar default constructor passed")

	def test_lidar_setget_uncer(self):
		lidar = pylisa.Lidar()
		lidar.set_ran_uncer(1.0)
		self.assertEqual(lidar.get_ran_uncer(), 1.0)
		print("Test Lidar set/get for range uncertainy passed")

	def test_lidar_setget_max_range(self):
		lidar = pylisa.Lidar()
		lidar.set_max_range(1.0)
		self.assertEqual(lidar.get_max_range(), 1.0)
		print("Test Lidar set/get for max range passed")

	def test_lidar_setget_min_range(self):
		lidar = pylisa.Lidar()
		lidar.set_min_range(1.0)
		self.assertEqual(lidar.get_min_range(), 1.0)
		print("Test Lidar set/get for min range passed")

	def test_lidar_setget_wavelength(self):
		lidar = pylisa.Lidar()
		lidar.set_wavelength(1.0)
		self.assertEqual(lidar.get_wavelength(), 1.0)
		print("Test Lidar set/get for wavelength passed")

	def test_lidar_setget_b_div(self):
		lidar = pylisa.Lidar()
		lidar.set_b_div(1.0)
		self.assertEqual(lidar.get_b_div(), 1.0)
		print("Test Lidar set/get for beam divergence passed")

if __name__=="__main__":
	unittest.main()
