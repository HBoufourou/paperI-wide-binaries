import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from data_source import load_catalog  # noqa: E402
from moteur_population import triple_vtilde  # noqa: E402


class CatalogueSafetyTests(unittest.TestCase):
    def test_newton_filename_is_refused_before_reading(self):
        with self.assertRaisesRegex(ValueError, "Newton_\\*"):
            load_catalog(ROOT / "Newton_dr3_MSMS_d200pc_5.csv", download=False)


class TripleGravityTests(unittest.TestCase):
    def setUp(self):
        self.sample = {
            "vt_outer_x": np.array([0.3, -0.2]),
            "vt_outer_y": np.array([0.4, 0.1]),
            "vt_inner_x": np.array([0.05, 0.02]),
            "vt_inner_y": np.array([-0.03, 0.04]),
        }

    def test_gamma_one_reconstructs_vector_sum(self):
        expected = np.hypot(
            self.sample["vt_outer_x"] + self.sample["vt_inner_x"],
            self.sample["vt_outer_y"] + self.sample["vt_inner_y"],
        )
        np.testing.assert_allclose(triple_vtilde(self.sample, 1.0), expected)

    def test_gamma_scales_only_outer_orbit(self):
        gamma = 1.4
        expected = np.hypot(
            np.sqrt(gamma) * self.sample["vt_outer_x"] + self.sample["vt_inner_x"],
            np.sqrt(gamma) * self.sample["vt_outer_y"] + self.sample["vt_inner_y"],
        )
        np.testing.assert_allclose(triple_vtilde(self.sample, gamma), expected)


if __name__ == "__main__":
    unittest.main()
