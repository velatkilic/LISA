import pylisa
import numpy as np

import unittest

class TestUtilsMethods(unittest.TestCase):
	def test_trapz(self):
		x = np.linspace(0,1,101)
		f = 2*x
		y = pylisa.trapz(f,x)
		self.assertEqual(y,1)
		print("Test trapz passed")

	def test_logspace(self):
		x  = pylisa.logspace(0,1,10)
		x2 = np.logspace(0,1,10)
		self.assertEqual(np.array_equal(x,x2), True)
		print("Test logspace passed")

	def test_calc_qext_with_real_index(self):
		xm   = 2.0 * 3.14159265359 * .525 / .6328;
		qext = pylisa.calc_qext(xm, 1.55);
		self.assertAlmostEqual(qext, 3.10543, 5)
		print("Test calc_qext with real index passed")

	def test_calc_qext_with_imag_index(self):
		qext = pylisa.calc_qext(1.0, 1.5 + 1.0j);
		self.assertAlmostEqual(qext, 2.336, 3)
		print("Test calc_qext with imaginary index passed")


if __name__=="__main__":
	unittest.main()
