import os, sys, unittest
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import selftest_breakpoints as M


class T(unittest.TestCase):
    def test_selftest_breakpoints(self):
        M.run()


if __name__ == '__main__':
    unittest.main()
