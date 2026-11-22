"""Four-loop 4PI stationary vertex and propagator, K=0 CTP pilot.

Use the symmetric normalized pair basis. V includes all contour weights
so convolutions have ordinary finite-index sums. The sunset self-energy
uses one bare and one full vertex as in Carrington et al. (arXiv:1310.4352),
Eqs. (10)–(13), after continuation to G=<T_C q q> and g=lambda/L.
"""
import itertools
import numpy as np
from born_bubble import ctp_mode
from two_pi_k0 import EV

def solve(lam,h,L=20.,m=1.,T=1.,tmin=-4.,tmax=1.,tol=2e-13,maxiter=50):
    n=round((tmax-tmin)/h)
    assert abs(n*h-tmax+tmin)<1e-10
    tt=np.linspace(tmin,tmax,n+1)
    idx=np.array([round((t-tmin)/h) for t in EV])
    assert np.max(np.abs(tt[idx]-EV))<1e-10
    w=np.full(n+1,h); w[[0,-1]]*=.5
    eta_g=np.concatenate((w,-w))*np.tile(lam/(L*np.cosh(tt/T)**2),2)
    bare=-1j*eta_g
    g0=np.block([[ctp_mode(tt[:,None],tt[None,:],a,b,np.array([m]))[...,0] for b in (0,1)] for a in (0,1)])
    nt=g0.shape[0]
    i,j=np.triu_indices(nt)
    scale=np.sqrt(1+(i!=j))
    norm=np.outer(scale,scale)
    pair=np.full((nt,nt),-1,int);pair[i,j]=pair[j,i]=np.arange(len(i))
    diag=pair[np.arange(nt),np.arange(nt)]
    # All recombinations use the same pair normalization:
    it=np.ix_(i,i); iu=np.ix_(i,j); jt=np.ix_(j,j); ju=np.ix_(j,i)
    ik,jl,il,jk=pair[it],pair[jt],pair[iu],pair[ju]
    v0=np.zeros((len(i),len(i)),complex);v0[diag,diag]=bare
    g=g0.copy();v=v0.copy()
    for iteration in range(maxiter):
        # On a normalized symmetric pair, G² matrix has the two symmetrized products.
        pairgreen=(g[it]*g[jt]+g[iu]*g[ju])*(.5*norm)
        chan=(v@pairgreen@v)*(.5/norm)
        vnew=v0+(chan+chan[ik,jl]+chan[il,jk])*norm
        # Connected self energy from one weighted bare vertex and a weighted full one.
        gproj=g[:,i]*g[:,j]*scale
        conv=gproj@vnew.T
        sig=(1/6)*bare[:,None]*np.einsum('xyc,xc,yc->xy',conv[:,pair],g,1/scale[pair],optimize=False)
        lo=bare*.5*(np.diag(g)-.5/m)
        gnew=g0+g0@(lo[:,None]*g+sig@g)
        resv=float(np.max(np.abs(vnew-v)));resg=float(np.max(np.abs(gnew-g)))
        g,v=gnew,vnew
        if max(resv,resg)<tol:break
    else:raise RuntimeError(f'4PI vertex={resv:g} and Dyson={resg:g} updates did not converge.')
    c4={}
    for bs in itertools.product((0,1),repeat=4):
        a,b,c,d=idx+np.array(bs)*(n+1)
        def proj(a,b):return .5*scale*(g[a,i]*g[b,j]+g[a,j]*g[b,i])
        c4[''.join('+' if b==0 else '-' for b in bs)]=proj(a,b)@v@proj(c,d)/L**2
    return dict(iter=iteration+1,v_res=resv,g_res=resg,nt_phys=n+1,ct_pair_dimension=len(i),G0_delta=float(np.max(np.abs(g-g0))),F4c=c4)
