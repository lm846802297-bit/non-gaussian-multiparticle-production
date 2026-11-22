import json
import math
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from sd_residual import make_hamiltonian, projected_normal_ordered_cubic, field_matrix

RESULTS = ROOT / "results" / "phase4-sd-001"


def momentum_matrix(modes, basis, index, mass, length, x):
    rows, cols, vals = [], [], []
    for j, state in enumerate(basis):
        for i, nmode in enumerate(modes):
            p = 2.0 * math.pi * nmode / length
            w = math.sqrt(p*p + mass*mass)
            c = math.sqrt(w/(2.0*length))
            if state[i]:
                ns = list(state); ns[i] -= 1
                row = index.get(tuple(ns))
                if row is not None:
                    rows.append(row); cols.append(j)
                    vals.append(-1j*c*np.exp(1j*p*x)*math.sqrt(state[i]))
            ns = list(state); ns[i] += 1
            row = index.get(tuple(ns))
            if row is not None:
                rows.append(row); cols.append(j)
                vals.append(1j*c*np.exp(-1j*p*x)*math.sqrt(state[i]+1))
    return coo_matrix((vals,(rows,cols)),shape=(len(basis),len(basis)),dtype=complex).tocsr()


class SDResidualTests(unittest.TestCase):
    def test_projected_cubic_matches_hamiltonian_commutator_away_from_fock_boundary(self):
        K, N, m, L, x = 1, 6, 1.0, 20.0, 0.37
        modes, basis, index, V, E, ws, herm = make_hamiltonian(K,N,m,L)
        pi = momentum_matrix(modes,basis,index,m,L,x)
        cubic = projected_normal_ordered_cubic(modes,basis,index,m,L,x)
        # For states with at most N-1 particles, no intermediate field
        # insertion reaches the Fock boundary. The projected canonical EOM is
        # -6 i [V, pi(x)] = P_K :phi^3(x):.
        comm_cubic = -6j*(V@pi-pi@V)
        ids = [i for i,s in enumerate(basis) if sum(s) <= N-1]
        diff = (comm_cubic-cubic)[ids,:][:,ids]
        ref = cubic[ids,:][:,ids]
        rel = np.linalg.norm(diff.data)/np.linalg.norm(ref.data)
        self.assertLess(rel, 2e-13)
        self.assertLess(herm, 1e-13)


    def test_normal_ordered_cubic_matches_fourier_projected_phi_cube_minus_tadpole(self):
        K, N, m, L, x0 = 1, 6, 1.0, 20.0, 0.37
        modes, basis, index, V, E, ws, herm = make_hamiltonian(K,N,m,L)
        delta = sum(1.0/(2.0*w*L) for w in ws)
        ngrid = 8*K+1
        cubes=[]
        for j in range(ngrid):
            x=L*j/ngrid
            phi=field_matrix(modes,basis,index,m,L,x)
            cubes.append(phi@phi@phi)
        coeffs=[]
        for q in range(-K,K+1):
            acc=None
            for j,cube in enumerate(cubes):
                term=(np.exp(-2j*np.pi*q*j/ngrid)/ngrid)*cube
                acc=term if acc is None else acc+term
            coeffs.append(acc)
        reconstructed=sum(coeffs[i]*np.exp(2j*np.pi*(-K+i)*x0/L)
                           for i in range(len(coeffs)))
        phi0=field_matrix(modes,basis,index,m,L,x0)
        explicit=reconstructed-3.0*delta*phi0
        normal=projected_normal_ordered_cubic(modes,basis,index,m,L,x0)
        ids=[i for i,state in enumerate(basis) if sum(state)<=N-3]
        diff=(explicit-normal)[ids,:][:,ids]
        ref=normal[ids,:][:,ids]
        rel=np.linalg.norm(diff.data)/np.linalg.norm(ref.data)
        self.assertLess(rel,2e-12)

    def test_sd_residual_converges_with_time_stencil(self):
        coarse = json.loads((RESULTS/'spectral-n8-h0.00625-a0.2.json').read_text())
        fine = json.loads((RESULTS/'spectral-n8-h0.003125-a0.2.json').read_text())
        self.assertLess(fine['SD_relative_residual'], coarse['SD_relative_residual']/3)
        self.assertLess(fine['SD_relative_residual'], 7e-6)

    def test_spatial_stencil_scan_is_stable(self):
        a05 = json.loads((RESULTS/'spectral-n8-h0.00625-a0.05.json').read_text())
        a10 = json.loads((RESULTS/'spectral-n8-h0.00625-a0.1.json').read_text())
        a20 = json.loads((RESULTS/'spectral-n8-h0.00625-a0.2.json').read_text())
        for item in (a05,a10,a20):
            fd=np.array(item['dxx_F4c_fd'])
            sp=np.array(item['dxx_F4c_spectral'])
            rel=np.linalg.norm(fd-sp)/np.linalg.norm(sp)
            self.assertLess(rel, 4e-4)
        err05=np.linalg.norm(np.array(a05['dxx_F4c_fd'])-np.array(a05['dxx_F4c_spectral']))
        err10=np.linalg.norm(np.array(a10['dxx_F4c_fd'])-np.array(a10['dxx_F4c_spectral']))
        err20=np.linalg.norm(np.array(a20['dxx_F4c_fd'])-np.array(a20['dxx_F4c_spectral']))
        self.assertLess(err05,err10)
        self.assertLess(err10,err20)

    def test_fock_cutoff_convergence_of_source_and_four_point(self):
        n8 = json.loads((RESULTS/'spectral-n8-h0.00625-a0.2.json').read_text())
        n10 = json.loads((RESULTS/'spectral-n10-h0.00625-a0.2.json').read_text())
        for key in ('F4c','M33_projected','connected_source'):
            x=np.array(n8[key]); y=np.array(n10[key])
            self.assertLess(np.linalg.norm(x-y)/np.linalg.norm(y), 2e-8)

    def test_free_vacuum_connected_four_point_and_sd_residual(self):
        free=json.loads((RESULTS/'free-vacuum.json').read_text())
        self.assertLess(np.linalg.norm(free['F4c']), 1e-14)
        self.assertLess(free['SD_residual_abs'], 1e-12)
        self.assertAlmostEqual(free['vacuum_norm_out'],1.0,places=12)


if __name__ == '__main__':
    unittest.main(verbosity=2)