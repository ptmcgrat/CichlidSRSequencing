import os, sys, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import selftest_bc as M


class T(unittest.TestCase):
    def test_selftest_bc(self):
        M.main()


if __name__ == '__main__':
    unittest.main()
