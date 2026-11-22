"""Three-loop 2PI + full bilocal Bethe–Salpeter response in finite momentum modes.

Uses all total pair momenta Q=k+l of the left pair. Each Q is closed under
the kernel derivative, even though background G is diagonal in momentum.
"""
from collections import Counter
from itertools import product
import numpy as np
from born_bubble import ctp_mode, LABELS

def solve(lam,h,events,K=1,L=20.,m=1.,T=1.,tmin=-4.,tol=2e-13,maxiter=150):
    stop=max(t for t,_ in events)
    steps=round((stop-tmin)/h)
    assert abs(steps*h-stop+tmin)<1e-9
    tt=np.linspace(tmin,stop,steps+1)
    it=np.array([round((t-tmin)/h) for t,_ in events])
    assert np.max(np.abs(tt[it]-[t for t,_ in events]))<1e-9
    nb=len(tt)
    w=np.full(nb,h); w[[0,-1]]*=.5
    eta_g=np.concatenate((w,-w))*np.tile(lam/(L*np.cosh(tt/T)**2),2)
    eta_gg=np.outer(eta_g,eta_g)
    modes=list(range(-K,K+1)); freqs=np.hypot(2*np.pi*np.array(modes)/L,m)
    g0_ij=[]
    for a in (0,1):
        g0_ij.append([ctp_mode(tt[:,None],tt[None,:],a,b,freqs).transpose((2,0,1)) for b in (0,1)])
    g0=np.concatenate([np.concatenate(row,axis=2) for row in g0_ij],axis=1)
    # Precompute ordered momentum triples, grouped only to avoid duplicate multiplies.
    triples={k:Counter(tuple(sorted((a,b,c))) for a,b,c in product(modes,repeat=3) if (a+b+c)==k)
             for k in modes}
    twolines={P:Counter(tuple(sorted((a,b))) for a,b in product(modes,repeat=2) if a+b==P)
              for P in range(-2*K,2*K+1)}
    counterterms=(.5/freqs).sum()
    g=g0.copy()
    for iteration in range(maxiter):
        hartree=-.5j*eta_g*(np.diagonal(g,axis1=1,axis2=2).sum(axis=0)-counterterms)
        nonlocal_terms=[]
        for k in modes:
            loop=sum(mult*g[a+K]*g[b+K]*g[c+K] for (a,b,c),mult in triples[k].items())
            nonlocal_terms.append(-(1/6)*eta_gg*loop)
        gnew=g0+g0@(hartree[None,:,None]*g+np.asarray(nonlocal_terms)@g)
        diff=float(np.max(np.abs(gnew-g)));g=gnew
        if diff<tol:break
    else:raise RuntimeError(f'Dyson did not converge at lam={lam}, h={h}, err={diff}.')
    loops={P:sum(v*g[a+K]*g[b+K] for (a,b),v in ab.items()) for P,ab in twolines.items()}
    ph=np.exp(2j*np.pi*np.asarray([x for _,x in events])[:,None]*np.array(modes)/L)
    by_Q={};block_its={};max_bs=0.;max_forbidden=0
    for Q in range(-2*K,2*K+1):
        pairs=[(k,l) for k,l in product(modes,repeat=2) if k+l==Q]
        # X and its kernel derivative never require r+s != Q; explicit assert.
        routes=[]
        for k,l in pairs:
            dst=[]
            for P in loops:
                r,s=k-P,l+P
                if r in modes and s in modes:
                    max_forbidden=max(max_forbidden,abs(r+s-Q))
                    dst.append((P,pairs.index((r,s))))
            routes.append(dst)
        gl=np.asarray([g[k+K] for k,l in pairs]);gr=np.asarray([g[l+K] for k,l in pairs])
        phases=np.asarray([ph[0,k+K]*ph[1,l+K] for k,l in pairs])
        components={};block_its[Q]={}
        for a,b in product((0,1),repeat=2):
            e3,e4=it[2]+a*nb,it[3]+b*nb
            x0=np.stack([np.outer(g[k+K,:,e3],g[l+K,e4,:])*np.conj(ph[2,k+K]*ph[3,l+K])
                       +np.outer(g[k+K,:,e4],g[l+K,e3,:])*np.conj(ph[3,k+K]*ph[2,l+K]) for k,l in pairs])
            c=np.zeros_like(x0)
            for it_bs in range(maxiter):
                xx=x0+c
                loc=-.5j*eta_g*np.diagonal(xx,axis1=1,axis2=2).sum(axis=0)
                nl=np.stack([-.5*eta_gg*sum(loops[P]*xx[ridx] for P,ridx in dst) for dst in routes])
                nxt=gl@(loc[None,:,None]*gr+nl@gr)
                err=float(np.max(np.abs(nxt-c)));c=nxt
                if err<tol:break
            else:raise RuntimeError(f'BS did not converge at lam={lam}, h={h}, Q={Q}, err={err}.')
            max_bs=max(max_bs,err)
            block_its[Q][f'{a}{b}']=it_bs+1
            for c1,c2 in product((0,1),repeat=2):
                tag=''.join('+' if p==0 else '-' for p in (c1,c2,a,b))
                components[tag]=(phases@c[:,it[0]+c1*nb,it[1]+c2*nb])/L**2
        by_Q[Q]=components
    assert max_forbidden==0
    comp={tag:sum(v[tag] for v in by_Q.values()) for tag in LABELS}
    return {'F4c':comp,'F4c_by_Q':by_Q,'step':h,'nt_phys':nb,'dyson_iterations':iteration+1,
            'dyson_residual':diff,'bs_residual':max_bs,'bs_iterations':block_its,
            'G0_delta':float(np.max(np.abs(g-g0))),'momentum_route_error':max_forbidden}
