"""
All steps S1 -> S5 on synthetic data (the --selftest scenarios).
"""

import unittest

import helpers  # noqa: F401  (puts the package on sys.path)
from barcodex.selftest import run_scenario


class TestAllSteps(unittest.TestCase):

    def test_pos_scenario(self):
        self.assertTrue(run_scenario("pos"))

    def test_flank_scenario(self):
        self.assertTrue(run_scenario("flank"))


if __name__ == "__main__":
    unittest.main()
