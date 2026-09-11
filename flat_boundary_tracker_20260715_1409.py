import numpy as np
from scipy.linalg import expm

BAND = 0.02

def make_step(N, dt):
    A = np.array([[0.0,1.0],[0.0,-N]]); B = np.array([0.0,1.0])
    M = np.zeros((3,3)); M[:2,:2]=A; M[:2,2]=B
    E = expm(M*dt)
    return E[:2,:2].copy(), E[:2,2].copy()

def sim_vec(k, tau, N, dt=5e-4, Tend=10.0):
    """Reduced loop L=k(tau s+1)e^{-s}/[s(1+s/N)], vectorized over (k,tau) pairs.
    Returns Ts (sub-grid refined), dymin (exact derivative, parabola-refined), y' at band entry."""
    k = np.atleast_1d(np.asarray(k,float)); tau = np.atleast_1d(np.asarray(tau,float))
    G = len(k)
    Ad, Bd = make_step(N, dt)
    n = int(round(Tend/dt)); nd = int(round(1.0/dt))
    c1, c2 = k*N, k*N*tau
    ub = np.zeros((n+1,G)); vb = np.zeros((n+1,G))  # u and u' histories
    x = np.zeros((2,G))
    for i in range(n):
        y = ub[i-nd] if i>=nd else np.zeros(G)
        e = 1.0 - y
        ub[i] = c1*x[0] + c2*x[1]
        # u' = C(Ax+Be) = c1*x2 + c2*(e - N x2)
        vb[i] = c1*x[1] + c2*(e - N*x[1])
        x = Ad@x + Bd[:,None]*e
    ub[n] = c1*x[0] + c2*x[1]
    # y(t)=u(t-1), y'(t)=v(t-1); build on grid indices i>=nd
    y = ub[:n+1-nd]; v = vb[:n+1-nd]      # y[i] corresponds to t=(i+nd)*dt
    t0 = nd*dt
    Ts = np.full(G, np.nan); dyTs = np.full(G, np.nan); dymin = np.full(G, np.nan)
    for g in range(G):
        yy = y[:,g]; vv = v[:,g]
        out = np.abs(1.0-yy) > BAND
        if not out.any() or out[-1]:
            continue
        j = np.max(np.where(out))          # last sample outside band
        # refine crossing y=0.98 between j and j+1 (monotone from below)
        y0,y1 = yy[j], yy[j+1]
        frac = (0.98-y0)/(y1-y0) if y1!=y0 else 0.0
        Ts[g] = t0 + (j+frac)*dt
        dyTs[g] = vv[j] + frac*(vv[j+1]-vv[j])
        # min y' after motion starts, parabola-refined
        i0 = 5
        m = i0 + int(np.argmin(vv[i0:]))
        if 0 < m < len(vv)-1:
            a,b,c = vv[m-1],vv[m],vv[m+1]
            denom = (a-2*b+c)
            dymin[g] = b - (a-c)**2/(8*denom) if abs(denom)>0 else b
        else:
            dymin[g] = vv[m]
    return Ts, dymin, dyTs

def kstar_bisect(tau_grid, N, klo=0.2, khi=2.5, iters=34, dt=5e-4):
    """Boundary k*(tau): max k with dymin >= 0, bisected vectorized over tau grid."""
    G = len(tau_grid)
    lo = np.full(G, klo); hi = np.full(G, khi)
    # ensure hi infeasible, lo feasible
    for _ in range(iters):
        mid = 0.5*(lo+hi)
        Ts, dm, _ = sim_vec(mid, tau_grid, N, dt=dt)
        feas = (dm >= 0.0) & ~np.isnan(Ts)
        lo = np.where(feas, mid, lo); hi = np.where(feas, hi, mid)
    Ts, dm, dyTs = sim_vec(lo, tau_grid, N, dt=dt)
    return lo, Ts, dm, dyTs

def solve_N(N, tau_lo=0.20, tau_hi=0.55, coarse=24, dt=5e-4):
    taus = np.linspace(tau_lo, tau_hi, coarse)
    ks, Ts, dm, dyTs = kstar_bisect(taus, N, dt=dt)
    j = int(np.nanargmin(Ts))
    # local refine around best tau
    w = (tau_hi-tau_lo)/(coarse-1)
    taus2 = np.linspace(max(tau_lo,taus[j]-1.5*w), min(tau_hi,taus[j]+1.5*w), 21)
    ks2, Ts2, dm2, dyTs2 = kstar_bisect(taus2, N, dt=dt)
    j2 = int(np.nanargmin(Ts2))
    return dict(N=N, tau=taus2[j2], k=ks2[j2], Ts=Ts2[j2], dymin=dm2[j2], dyTs=dyTs2[j2])

if __name__ == '__main__':
    import sys, json, time
    Ns = [float(x) for x in sys.argv[1:]]
    out = []
    for N in Ns:
        t0=time.time()
        r = solve_N(N)
        r['sec'] = round(time.time()-t0,1)
        out.append(r)
        print(f"N={N:6.1f}  Ts*={r['Ts']:.6f}  k*={r['k']:.6f}  tau*={r['tau']:.4f}  "
              f"miny'={r['dymin']:.2e}  y'(Ts)={r['dyTs']:.4f}  [{r['sec']}s]", flush=True)
    json.dump(out, open(f'flat_{int(Ns[0])}_{int(Ns[-1])}.json','w'))
