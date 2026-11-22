import json
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results"
P = R / "phase6-ctp-sd-001"
K = R / "phase6-kcut-scan-001"


def load(path):
    return json.loads(path.read_text())


class CTPAllBranchSDTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a8 = load(P / "all16-n8.json")
        cls.a10 = load(P / "all16-n10.json")
        cls.c8 = load(P / "config-C-all16-n8.json")
        cls.free = load(P / "free-vacuum-all16.json")

    def test_all_16_components_have_small_sd_residual(self):
        for data in (self.a8, self.c8, self.a10):
            self.assertEqual(data["checks"]["number_of_components"], 16)
            for comp in data["components"].values():
                self.assertLess(comp["SD_relative_residual"], 7e-6)

    def test_ctp_largest_time_identity(self):
        for data in (self.a8, self.c8, self.a10):
            self.assertLess(data["checks"]["max_largest_time_F4c_residual"], 1e-14)

    def test_complement_branch_conjugation(self):
        for data in (self.a8, self.c8, self.a10, self.free):
            self.assertLess(data["checks"]["max_F4c_complement_conjugation_error"], 1e-14)
            self.assertLess(data["checks"]["max_residual_complement_conjugation_error"], 1e-14)

    def test_all_forward_matches_independent_time_ordered_engine(self):
        old = load(R / "phase4-sd-001" / "spectral-n8-h0.003125-a0.2.json")
        new = self.a8["components"]["++++"]
        for key_old, key_new in (("F4c", "F4c"), ("M33_projected", "M33"),
                                 ("connected_source", "connected_composite_source"),
                                 ("SD_residual", "SD_residual")):
            x = np.asarray(old[key_old], dtype=float)
            y = np.asarray(new[key_new], dtype=float)
            self.assertLess(np.linalg.norm(x-y) / max(np.linalg.norm(x), 1e-300), 2e-10)

    def test_all_branch_fock_cutoff_convergence(self):
        a, b = self.a8["components"], self.a10["components"]
        for key in ("F4c", "M33", "connected_composite_source", "source_term"):
            changes = []
            for label in a:
                x = np.asarray(a[label][key]); y = np.asarray(b[label][key])
                changes.append(np.linalg.norm(x-y)/max(np.linalg.norm(y), 1e-300))
            self.assertLess(max(changes), 2e-8)

    def test_second_spacetime_configuration(self):
        self.assertEqual(self.c8["checks"]["number_of_components"], 16)
        self.assertLess(self.c8["checks"]["max_SD_relative_residual"], 7e-6)
        self.assertLess(self.c8["checks"]["max_largest_time_F4c_residual"], 1e-14)

    def test_free_vacuum_connected_four_point(self):
        for comp in self.free["components"].values():
            self.assertLess(np.linalg.norm(comp["F4c"]), 3e-16)
            self.assertLess(comp["SD_residual_abs"], 4e-12)

    def test_kcut_scan_results_are_archived_and_finite(self):
        for Kcut, Nmax in ((1, 5), (2, 5), (3, 5)):
            data = load(K / f"K{Kcut}-N{Nmax}.json")
            self.assertEqual(data["parameters"]["Kcut"], Kcut)
            self.assertTrue(np.all(np.isfinite(data["F4c"])))
            self.assertTrue(np.all(np.isfinite(data["M33_projected"])))
            self.assertTrue(np.all(np.isfinite(data["connected_source"])))
            self.assertLess(data["SD_relative_residual"], 3e-5)


if __name__ == "__main__":
    unittest.main(verbosity=2)