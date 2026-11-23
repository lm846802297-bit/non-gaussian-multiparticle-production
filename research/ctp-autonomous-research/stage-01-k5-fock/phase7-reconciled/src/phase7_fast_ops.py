#!/usr/bin/env python3
"""Faster exact sparse construction for the same finite-Fock operators.

This module preserves the monomial sums and projections of the Phase-VI
source, but enumerates only annihilation multisets supported by each ket.
"""
from __future__ import annotations
import itertools, math
from collections import Counter, defaultdict
import numpy as np
from scipy.sparse import coo_matrix
from sd_residual import make_basis


def _multiset_data(modes, u, max_order, include_factorials):
    nm=len(modes)
    by_order_momentum={}
    all_by_order={}
    for n in range(max_order+1):
        groups=defaultdict(list)
        rows=[]
        for tup in itertools.combinations_with_replacement(range(nm),n):
            c=Counter(tup)
            occ=tuple(c.get(i,0) for i in range(nm))
            mom=sum(modes[i]*v for i,v in enumerate(occ))
            if include_factorials:
                coeff=math.prod(u[i]**v/math.factorial(v) for i,v in enumerate(occ))
            else:
                coeff=math.prod(u[i]**v for i,v in enumerate(occ))
            item=(occ,coeff)
            groups[mom].append(item)
            rows.append((occ,mom,coeff))
        by_order_momentum[n]=groups
        all_by_order[n]=rows
    return by_order_momentum,all_by_order


def _apply_monomial(state, av, cv, coeff, index, col, rows, cols, vals):
    new=tuple(state[i]-av[i]+cv[i] for i in range(len(state)))
    row=index.get(new)
    if row is None:
        return
    amp=coeff
    for i in range(len(state)):
        a=av[i]; c=cv[i]
        if a:
            amp *= math.sqrt(math.factorial(state[i])/math.factorial(state[i]-a))
        if c:
            amp *= math.sqrt(math.factorial(state[i]-a+c)/math.factorial(state[i]-a))
    if amp:
        rows.append(row); cols.append(col); vals.append(amp)


def make_hamiltonian_fast(kcut,nmax,mass,length):
    modes,basis,index=make_basis(kcut,nmax); nm=len(modes)
    ws=[math.sqrt((2*math.pi*p/length)**2+mass**2) for p in modes]
    u=[1/math.sqrt(2*w*length) for w in ws]
    _, create_by_order=_multiset_data(modes,u,4,True)
    create_groups={n:defaultdict(list) for n in range(5)}
    for n in range(5):
        for cv,pc,fc in create_by_order[n]:
            create_groups[n][pc].append((cv,fc))
    ann_by_order={n:[(av,pa,fa) for av,pa,fa in create_by_order[n]] for n in range(5)}
    rows=[]; cols=[]; vals=[]
    for j,state in enumerate(basis):
        for na in range(5):
            nc=4-na
            for av,pa,fa in ann_by_order[na]:
                if any(state[i]<av[i] for i in range(nm)):
                    continue
                for cv,fc in create_groups[nc].get(pa,()):
                    _apply_monomial(state,av,cv,length*fc*fa,index,j,rows,cols,vals)
    V=coo_matrix((vals,(rows,cols)),shape=(len(basis),len(basis)),dtype=float).tocsr(); V.sum_duplicates()
    diff=V-V.T
    herm=float(np.max(np.abs(diff.data))) if diff.nnz else 0.0
    energies=np.array([sum(n*w for n,w in zip(s,ws)) for s in basis])
    return modes,basis,index,V,energies,ws,herm


def projected_normal_ordered_cubic_fast(modes,basis,index,mass,L,x):
    nm=len(modes)
    ws=[math.sqrt((2*math.pi*p/L)**2+mass**2) for p in modes]
    u=[1/math.sqrt(2*w*L) for w in ws]
    K=max(abs(p) for p in modes)
    _, all_multisets=_multiset_data(modes,u,3,False)
    create_by_order={n:defaultdict(list) for n in range(4)}
    ann_by_order={n:[] for n in range(4)}
    for n in range(4):
        for occ,mom,coeff in all_multisets[n]:
            create_by_order[n][mom].append((occ,coeff))
            ann_by_order[n].append((occ,mom,coeff))
    rows=[]; cols=[]; vals=[]
    for j,state in enumerate(basis):
        for na in range(4):
            nc=3-na
            for av,pa,fa in ann_by_order[na]:
                if any(state[i]<av[i] for i in range(nm)):
                    continue
                for pc in range(pa-K,pa+K+1):
                    Q=pa-pc
                    phase=np.exp(1j*(2*math.pi/L)*Q*x)
                    for cv,fc in create_by_order[nc].get(pc,()):
                        mult=math.factorial(3)/(math.prod(math.factorial(v) for v in cv)*math.prod(math.factorial(v) for v in av))
                        _apply_monomial(state,av,cv,mult*fc*fa*phase,index,j,rows,cols,vals)
    O=coo_matrix((vals,(rows,cols)),shape=(len(basis),len(basis)),dtype=np.complex128).tocsr(); O.sum_duplicates()
    return O
