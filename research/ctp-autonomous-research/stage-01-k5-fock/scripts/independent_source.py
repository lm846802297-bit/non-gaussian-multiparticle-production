"""Independent finite-Fock Wick identity, using sparse field products."""
import argparse, json, math, platform, sys, time
from pathlib import Path
import numpy as np
import scipy
from scipy.sparse import coo_matrix
from sd_residual import make_basis, projected_normal_ordered_cubic


def setup(K,N):
    modes,basis,_=make_basis(K,N+1)
    b=np.asarray(basis,dtype=np.int64)
    m=np.asarray(modes,dtype=np.int64)
    return b,b@m,np.asarray(basis,dtype=np.int64).sum(axis=1),m


def wick_matrix(b, p, occupations, modes, K, N, mass, L, x):
    # Build phi(0) by only populating raising transitions then taking adjoint;
    # independent of the repository's cubic normal-order expansion.
    dim = b.shape[0]
    d = int(np.searchsorted(occupations, N, side='right'))
    omega = np.sqrt(mass**2 + (2*np.pi*modes/L)**2)
    u = 1/np.sqrt(2*L*omega)
    place = np.asarray([(N+2)**i for i in range(len(modes))], dtype=np.int64)
    code = b@place
    order = np.argsort(code)
    valid = np.where(occupations <= N)[0]
    assert len(valid) == d and np.array_equal(valid,np.arange(d))
    rows, cols, vals = [], [], []
    for q in range(len(modes)):
        target = code[valid]+place[q]
        row = order[np.searchsorted(code[order],target)]
        assert np.array_equal(code[row],target)
        rows.extend(row); cols.extend(valid); vals.extend(u[q]*np.sqrt(b[valid,q]+1))
    a = coo_matrix((vals,(rows,cols)),shape=(dim,dim)).tocsr()
    f = a+a.T
    b1 = f[:,:d]
    # Exact for these endpoints: a 3-step ±1 path cannot exceed N+1.
    unw = (b1.T@(f@b1)-3*np.sum(u*u)*f[:d,:d]).tocoo()
    transfer = p[unw.col]-p[unw.row]
    keep = np.abs(transfer)<=K
    projected = coo_matrix((unw.data[keep]*np.exp(2j*np.pi*x*transfer[keep]/L),(unw.row[keep],unw.col[keep])),shape=(d,d)).tocsr()
    projected.eliminate_zeros()
    return projected


def summarize(m,tol=1e-14):
    return {'shape':m.shape,'nnz':m.nnz,'nnz_over_tolerance':int(np.count_nonzero(np.abs(m.data)>tol)),
            'max_abs':float(np.max(np.abs(m.data),initial=0)),
            'frobenius':float(np.linalg.norm(m.data)), 'tolerance':tol}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--K',type=int,default=5); p.add_argument('--N',type=int,default=6)
    p.add_argument('--x',type=float,default=0.0); p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    tic=time.monotonic()
    b, mom, occupations, modes = setup(args.K,args.N)
    w=wick_matrix(b,mom,occupations,modes,args.K,args.N,1.,20.,args.x)
    basis=list(map(tuple,b[:w.shape[0]]))
    ref=projected_normal_ordered_cubic(tuple(modes),basis,{s:j for j,s in enumerate(basis)},1.,20.,args.x)
    d=w-ref
    result={'comparison':summarize(d),'reference':summarize(ref),
            'independent':summarize(w),'independent_antihermitian':summarize(w-w.getH()),
            'vacuum_to_1particle_max_abs':float(np.max(np.abs(w[:len(modes)+1,0].toarray()))),
            'relative_frobenius':float(np.linalg.norm(d.data)/np.linalg.norm(ref.data)),
            'K':args.K,'Nmax':args.N,'m':1.0,'L':20.0,'x':args.x,
            'intermediate_basis_dimension':b.shape[0],
            'source_commit':'eab1f5b3198d03eafc3f162299c9038956beaa0b',
            'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),
            'elapsed_s':time.monotonic()-tic}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,indent=2,sort_keys=True))
    if not result['comparison']['nnz_over_tolerance']==0: raise SystemExit(1)
if __name__=='__main__':main()
