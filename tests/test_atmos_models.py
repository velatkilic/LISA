import unittest
from pylisa.atmos_models import LISA

class TestLisa(unittest.TestCase):
    def test_lisa_constructor(self):
        lisa = LISA()
        self.assertIsNotNone(lisa)

if __name__ == "__main__":
    unittest.main()