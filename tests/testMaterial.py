import pylisa
import numpy as np

import unittest

class TestMaterial(unittest.TestCase):
	def test_water_default_constructor(self):
		water = pylisa.Water()
		self.assertNotEqual(water, None)
		print("Test water default constructor passed")

	def test_water_ref_ind(self):
		water = pylisa.Water()
		self.assertAlmostEqual(water.ref_ind(0.5876), 1.3335, 4)
		print("Test water refractive index passed")

	def test_material_constructor(self):
		material = pylisa.Material()
		self.assertNotEqual(material, None)
		print("Test material default constructor passed")


if __name__=="__main__":
	unittest.main()
