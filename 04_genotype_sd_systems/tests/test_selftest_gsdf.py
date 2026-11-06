import os, sys, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import selftest_gsdf as M


class T(unittest.TestCase):
    def test_selftest_gsdf(self):
        M.run()


if __name__ == '__main__':
    unittest.main()
