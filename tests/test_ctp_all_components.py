import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ctp_all_components import ctp_components


class CTPAllComponentsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.events = [(0.8, 0.0), (0.3, 0.5), (-0.2, 1.0), (-0.7, 1.5)]

    def test_all_sixteen_branch_components_exist(self):
        r = ctp_components(self.events, kcut=1, nmax=6, lam=0.1)
        self.assertEqual(r["checks"]["number_of_components"], 16)
        self.assertEqual(len(r["components"]), 16)

    def test_complement_branch_components_are_conjugates(self):
        r = ctp_components(self.events, kcut=1, nmax=6, lam=0.1)
        self.assertLess(r["checks"]["max_complement_conjugation_error"], 2e-9)

    def test_free_vacuum_connected_components_vanish(self):
        r = ctp_components(self.events, kcut=1, nmax=6, lam=0.0)
        for item in r["components"].values():
            z = complex(*item["F4_connected"])
            self.assertLess(abs(z), 2e-9)

    def test_interaction_generates_nonzero_connected_component(self):
        r = ctp_components(self.events, kcut=1, nmax=6, lam=0.1)
        values = [abs(complex(*x["F4_connected"])) for x in r["components"].values()]
        self.assertGreater(max(values), 1e-8)

    def test_largest_time_branch_independence(self):
        r = ctp_components(self.events, kcut=1, nmax=6, lam=0.1)
        c = r["components"]
        # Event 0 has the latest physical time. Its + and - components agree
        # for every fixed assignment of the other three branch indices.
        errors = []
        for mask in range(8):
            label_minus = "".join("+" if (mask >> i) & 1 else "-" for i in range(4))
            label_plus = "".join("+" if ((mask | 1) >> i) & 1 else "-" for i in range(4))
            zm = complex(*c[label_minus]["F4_connected"])
            zp = complex(*c[label_plus]["F4_connected"])
            errors.append(abs(zm-zp))
        self.assertLess(max(errors), 2e-9)

    def test_equal_time_components_are_branch_independent(self):
        events = [(0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 1.0)]
        r = ctp_components(events, kcut=1, nmax=6, lam=0.1)
        values = [complex(*x["F4_connected"]) for x in r["components"].values()]
        self.assertLess(max(abs(z-values[0]) for z in values), 2e-9)


if __name__ == "__main__":
    unittest.main()
