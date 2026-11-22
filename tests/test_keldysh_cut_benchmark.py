import unittest

from src.keldysh_cut_benchmark import (
    lambda_tilde,
    n2_cut_sum,
    n2_fock_sum,
    simpson_fourier,
)


class KeldyshCutBenchmarkTests(unittest.TestCase):
    def test_cut_sum_matches_unordered_fock_sum(self):
        for nk, nq in ((1, 2), (1, -1), (1, 1)):
            for cutoff in (4, 8, 16, 24):
                with self.subTest(nk=nk, nq=nq, cutoff=cutoff):
                    cut = n2_cut_sum(nk, nq, cutoff)
                    fock = n2_fock_sum(nk, nq, cutoff)
                    self.assertAlmostEqual(cut, fock, delta=5e-14 * abs(cut))
                    self.assertGreaterEqual(cut, 0.0)

    def test_coincident_mode_uses_factorial_moment(self):
        # The k=q observable is <n_k(n_k-1)>, not <n_k^2>.
        cut = n2_cut_sum(1, 1, 24)
        fock = n2_fock_sum(1, 1, 24)
        self.assertAlmostEqual(cut, fock, delta=5e-14 * abs(cut))

    def test_analytic_sech_squared_transform(self):
        for energy in (2.0, 4.0, 6.0, 8.0):
            with self.subTest(energy=energy):
                analytic = lambda_tilde(energy)
                numeric = abs(simpson_fourier(energy))
                self.assertLess(abs(numeric-analytic)/analytic, 1e-10)


if __name__ == "__main__":
    unittest.main()
