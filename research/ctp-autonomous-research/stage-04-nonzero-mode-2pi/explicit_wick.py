"""Independent coordinate-space Wick sums on a discrete time/space contour."""
import itertools
import numpy as np

def wick(events, h, K=1,L=20.,m=1.,T=1.,tmin=-4.,selected=None):
    t=np.linspace(tmin,max(z[0] for z in events),round((max(z[0] for z in events)-tmin)/h)+1)
    p=2*np.pi*np.arange(-K,K+1)/L; omega=np.hypot(p,m)
    nx=4*K+1
    s=np.arange(nx)*L/nx
    ww=np.full(len(t),h);ww[[0,-1]]*=.5
    ns=np.tile(s,2*len(t)); nt=np.repeat(np.tile(t,2),nx)
    eta=np.repeat([1,-1],len(t)*nx)
    w=(eta*np.tile(np.repeat(ww,nx),2)*(L/nx))/np.cosh(nt/T)**2
    def c(t,s,x,y,a,b):
        d=t-s
        f=np.where(a==b,np.abs(d)*a,np.where(a==-1,d,-d))
        return (np.exp(-1j*f[...,None]*omega+1j*(x-y)[...,None]*p)/(2*L*omega)).sum(axis=-1)
    g=c(nt[:,None],nt[None,:],ns[:,None],ns[None,:],eta[:,None],eta[None,:])
    g2=g*g
    ee=[[c(ti,nt,xi,ns,a,eta) for a in (1,-1)] for ti,xi in events]
    out={}
    for bs in itertools.product((0,1),repeat=4):
        tag=''.join('+' if a==0 else '-' for a in bs)
        if selected is not None and tag not in selected:continue
        f=[ee[j][b] for j,b in enumerate(bs)]
        c1=-1j*np.dot(w,np.prod(f,axis=0))
        c2=-.5*sum((w*f[0]*f[j])@g2@(w*f[k]*f[l]) for j,k,l in ((1,2,3),(2,1,3),(3,1,2)))
        out[tag]=c1,c2
    return out

def wick_by_q(events,h,K=1,L=20.,m=1.,T=1.,tmin=-4.,selected=None):
    """Project total left pair momentum by translating the left two points."""
    shifts=np.arange(4*K+1)*L/(4*K+1)
    f=[wick([(ti,xi+(d if j<2 else 0.)) for j,(ti,xi) in enumerate(events)],h,K,L,m,T,tmin,selected) for d in shifts]
    return {Q:{tag:tuple(np.asarray([g[tag] for g in f]).T @ np.exp(-2j*np.pi*Q*shifts/L)/len(shifts)) for tag in f[0]}
            for Q in range(-2*K,2*K+1)}
