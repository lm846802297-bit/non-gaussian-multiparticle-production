import json
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "phase4-multiconfig-001"

class MultiConfigurationSDTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = {p.stem: json.loads(p.read_text()) for p in RESULTS.glob("*.json")}

    def test_all_four_configurations_have_small_residual(self):
        for name in ("A_reference", "B_spatial_translation", "C_changed_separations_times", "D_early_pulse_side"):
            self.assertLess(self.data[name]["SD_relative_residual"], 7e-6)

    def test_spatial_translation_invariance_of_correlators_and_source(self):
        a, b = self.data["A_reference"], self.data["B_spatial_translation"]
        for key in ("F4c", "M33_projected", "connected_source", "source_term"):
            x, y = np.asarray(a[key]), np.asarray(b[key])
            self.assertLess(np.linalg.norm(x-y)/np.linalg.norm(x), 1e-10)

    def test_time_stencil_refinement(self):
        coarse = self.data["A_reference"]["SD_relative_residual"]
        fine = self.data["A_N8_h0p0015625"]["SD_relative_residual"]
        self.assertLess(fine, coarse/3.0)

    def test_fock_cutoff_refinement(self):
        a, b = self.data["A_reference"], self.data["A_N10_h0p003125"]
        for key in ("F4c", "M33_projected", "connected_source"):
            x, y = np.asarray(a[key]), np.asarray(b[key])
            self.assertLess(np.linalg.norm(x-y)/np.linalg.norm(y), 1e-8)

    def test_source_and_derivative_terms_cancel(self):
        for name in ("A_reference", "B_spatial_translation", "C_changed_separations_times", "D_early_pulse_side"):
            r = self.data[name]
            x = np.asarray(r["box_plus_mass_F4c"])
            y = np.asarray(r["source_term"])
            residual = np.asarray(r["SD_residual"])
            self.assertLess(np.linalg.norm(x+y-residual), 1e-18)

if __name__ == "__main__":
    unittest.main(verbosity=2)