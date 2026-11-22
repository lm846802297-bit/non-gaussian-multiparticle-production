import math
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from ctp_all_components import build_hamiltonian, field_matrix
from ctp_four_point_sd_residual import composite_phi3_matrix, evaluate_at


class CTPFourPointSDResidualTests(unittest.TestCase):
    def test_normal_ordered_phi3_matches_phi_cubed_subtraction_on_interior_columns(self):
        modes, basis, index, V, energies, omega, h_herm = build_hamiltonian(1, 6, 1.0, 20.0)
        phi = field_matrix(modes, basis, index, 1.0, 20.0, 0.0)
        a, a_herm = composite_phi3_matrix(modes, basis, index, 1.0, 20.0, 0.0)
        delta = sum(1.0 / (2.0 * w * 20.0) for w in omega)
        direct = phi @ phi @ phi - 3.0 * delta * phi
        interior = [i for i, state in enumerate(basis) if sum(state) <= 3]
        difference = (a - direct)[:, interior]
        error = np.max(np.abs(difference.data)) if difference.nnz else 0.0
        self.assertLess(error, 2e-14)
        self.assertLess(a_herm, 2e-14)

    def test_free_gaussian_vacuum_has_wick_composite_source(self):
        kcut, nmax, mass, length = 1, 4, 1.0, 20.0
        modes, basis, index, V, energies, omega, h_herm = build_hamiltonian(kcut, nmax, mass, length)
        psi0 = np.zeros(len(basis), dtype=np.complex128)
        psi0[index[tuple(0 for _ in modes)]] = 1.0
        template = [(1.0, 0.0), (0.5, 1.0), (-0.5, 2.0), (-1.0, 3.0)]
        f4c, k3c, comp_herm = evaluate_at(1.0, 0.0, template, psi0, -12.0, 12.0,
                                          modes, basis, index, energies, V,
                                          mass, length, 0.0, 1.0)
        self.assertLess(abs(f4c), 2e-9)
        # In a Gaussian vacuum the composite cumulant is not zero: all three
        # constituents of :phi^3: can contract with the three external fields.
        # Normal ordering removes internal contractions, leaving 3! F12 F13 F14.
        propagators = []
        for tj, xj in template[1:]:
            value = 0.0j
            for p in modes:
                momentum = 2.0 * math.pi * p / length
                w = math.sqrt(momentum * momentum + mass * mass)
                value += np.exp(-1j * w * (1.0 - tj) - 1j * momentum * xj) / (2.0 * w * length)
            propagators.append(value)
        expected = 6.0 * propagators[0] * propagators[1] * propagators[2]
        self.assertLess(abs(k3c - expected), 2e-12)
        self.assertLess(comp_herm, 2e-14)


if __name__ == "__main__":
    unittest.main()
