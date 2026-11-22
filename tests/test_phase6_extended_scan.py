import json
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results" / "phase6-kcut-extension-001"

def load(name):
    return json.loads((R / name).read_text())

def z(obj, key):
    return complex(*obj[key])

class Phase6ExtendedScanTests(unittest.TestCase):
    def test_prefix_cache_matches_original_reference_exactly(self):
        old = load("K1-N5-reference.json")
        new = load("K1-N5-cached-tight.json")
        for key in ("F4c", "M33_projected", "connected_source", "box_plus_mass_F4c", "SD_residual"):
            self.assertEqual(old[key], new[key], key)

    def test_k3_fock_cutoff_change_is_small_for_source(self):
        n5 = load("K3-N5-reference.json")
        n6 = load("K3-N6-h0p00625.json")
        f4_change = abs(z(n6, "F4c") - z(n5, "F4c")) / abs(z(n5, "F4c"))
        source_change = abs(z(n6, "connected_source") - z(n5, "connected_source")) / abs(z(n5, "connected_source"))
        self.assertLess(f4_change, 0.002)
        self.assertLess(source_change, 0.0001)

    def test_k3_time_stencil_refinement_reduces_residual(self):
        coarse = load("K3-N6-h0p00625.json")
        fine = load("K3-N6-h0p003125.json")
        self.assertLess(fine["SD_relative_residual"], coarse["SD_relative_residual"] / 3.0)

    def test_kcut_four_point_scan_is_not_monotonic(self):
        vals = [abs(z(load(f"K{k}-N5-matched-loose.json"), "F4c")) for k in (1, 2, 3)]
        vals.append(abs(z(load("K4-N5-h0p00625.json"), "F4c")))
        self.assertGreater(vals[1], vals[0])
        self.assertGreater(vals[2], vals[1])
        self.assertLess(vals[3], vals[2])

    def test_k4_nmax4_is_not_adequate_for_sd_residual(self):
        n4 = load("K4-N4-h0p00625.json")
        n5 = load("K4-N5-h0p00625.json")
        self.assertGreater(n4["SD_relative_residual"], 1e-3)
        self.assertLess(n5["SD_relative_residual"], 2e-5)
        # The h refinement at Nmax=5 is not yet monotonic; preserve this as an
        # explicit unresolved numerical limitation rather than a PASS criterion.
        n5_fine = load("K4-N5-h0p003125.json")
        self.assertGreater(n5_fine["SD_relative_residual"], n5["SD_relative_residual"])

    def test_all_reported_states_preserve_norm_and_hermiticity(self):
        names = [
            "K1-N5-matched-loose.json", "K2-N5-matched-loose.json",
            "K3-N5-matched-loose.json", "K3-N6-h0p00625.json",
            "K3-N6-h0p003125.json", "K4-N4-h0p00625.json",
            "K4-N5-h0p00625.json", "K4-N5-h0p003125.json"
        ]
        for name in names:
            r = load(name)
            self.assertLess(abs(r["vacuum_norm_out"] - 1.0), 1e-10, name)
            self.assertLess(r["interaction_hermiticity_error"], 1e-13, name)
            self.assertTrue(math.isfinite(r["SD_relative_residual"]), name)

if __name__ == "__main__":
    unittest.main(verbosity=2)