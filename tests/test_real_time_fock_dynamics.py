import math
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from real_time_fock_dynamics import (
    build_interaction, evolve, leading_n2, mode_observables,
    equal_time_connected_four_field,
)


class RealTimeFockDynamicsTests(unittest.TestCase):
    def test_projected_interaction_is_hermitian(self):
        result = build_interaction(kcut=2, nmax=6, mass=1.0, length=20.0)
        self.assertLess(result[6], 1e-13)

    def test_vacuum_to_four_particle_vertex_normalization(self):
        modes, basis, index, V, energies, ws, herm, nmon = build_interaction(
            kcut=2, nmax=4, mass=1.0, length=20.0
        )
        vacuum = index[tuple(0 for _ in modes)]
        target = (0, 1, 2, 1, 0)  # momenta -1,0,0,1,2 -> occupation {0,1,2}
        row = index[target]
        # Direct normalized-state coefficient: L prod u_p / sqrt(prod n_p!).
        direct = 20.0
        occ = target
        for n, w in zip(occ, ws):
            if n:
                direct *= (1.0 / math.sqrt(2.0 * w * 20.0)) ** n
                direct /= math.sqrt(math.factorial(n))
        self.assertAlmostEqual(V[row, vacuum], direct, places=14)

    def test_repeated_and_distinct_occupation_vertex_factors(self):
        kcut, nmax, mass, length = 3, 4, 1.0, 20.0
        modes, basis, index, V, energies, ws, herm, nmon = build_interaction(
            kcut, nmax, mass, length
        )
        vacuum = index[tuple(0 for _ in modes)]
        cases = [
            {0: 4},
            {-1: 2, 1: 2},
            {-1: 1, 0: 2, 1: 1},
            {-1: 3, 3: 1},
            {-2: 1, -1: 1, 1: 1, 2: 1},
        ]
        for occupations in cases:
            target = tuple(occupations.get(p, 0) for p in modes)
            row = index[target]
            direct = length
            for p, n in occupations.items():
                ip = modes.index(p)
                direct *= (1.0 / math.sqrt(2.0 * ws[ip] * length)) ** n
                direct /= math.sqrt(math.factorial(n))
            self.assertAlmostEqual(V[row, vacuum], direct, places=14)

    def test_norm_and_zero_coupling(self):
        r = evolve(2, 6, 1.0, 20.0, 0.0, 1.0, rtol=1e-9, atol=1e-11)
        self.assertAlmostEqual(r["norm"], 1.0, places=10)
        o = mode_observables(r["modes"], r["basis"], r["psi"], 1, 2)
        self.assertLess(abs(o["N2"]), 1e-24)

    def test_weak_coupling_matches_leading_perturbation(self):
        r = evolve(2, 8, 1.0, 20.0, 0.01, 1.0, rtol=1e-9, atol=1e-11)
        o = mode_observables(r["modes"], r["basis"], r["psi"], 1, 2)
        lo = leading_n2(2, 1.0, 20.0, 0.01, 1.0, 1, 2)
        self.assertLess(abs(o["N2"] / lo - 1.0), 0.01)

    def test_coincident_channel_is_factorial_moment(self):
        r = evolve(2, 8, 1.0, 20.0, 0.1, 1.0, rtol=1e-9, atol=1e-11)
        o = mode_observables(r["modes"], r["basis"], r["psi"], 1, 1)
        self.assertGreater(o["N2"], 0.0)
        self.assertGreater(o["kappa4"], 0.0)

    def test_free_vacuum_connected_four_field_vanishes(self):
        r = evolve(2, 6, 1.0, 20.0, 0.0, 1.0, save_times=3, rtol=1e-9, atol=1e-11)
        c = equal_time_connected_four_field(2, 6, 1.0, 20.0, r["states"][:, 1],
                                            (0.0, 0.0, 0.0, 2.0))
        self.assertLess(abs(c["G4_connected"][0]), 1e-13)
        self.assertLess(abs(c["G4_connected"][1]), 1e-13)

    def test_interaction_generates_equal_time_connected_four_field(self):
        r = evolve(2, 8, 1.0, 20.0, 0.2, 1.0, save_times=3, rtol=1e-9, atol=1e-11)
        c = equal_time_connected_four_field(2, 8, 1.0, 20.0, r["states"][:, 1],
                                            (0.0, 0.0, 0.0, 2.0))
        self.assertGreater(abs(c["G4_connected"][0]), 1e-6)


if __name__ == "__main__":
    unittest.main()