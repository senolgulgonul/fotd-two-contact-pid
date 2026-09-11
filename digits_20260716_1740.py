import numpy as np
from scipy.optimize import fsolve
exec(open('trace2_20260716_1740.py').read().split('# ---- Anchor')[0])

# Step-size convergence of the corner and of the 1%-band optimum (paper Section 8).
def mins(k, tau, N, dt):
    t, y, v = simulate(k, tau, N, dt=dt, Tend=6.5)
    return (v[(t>=2.5)&(t<=3.16)].min(), v[(t>=3.18)&(t<=3.99)].min(), v[(t>=4.02)&(t<=4.99)].min())
def res3(x, dt): return list(mins(x[0],x[1],x[2],dt))

for dt in [1e-4, 5e-5]:
    x = fsolve(res3, [0.61923,0.35326,21.35], args=(dt,), xtol=1e-11, epsfcn=1e-8)
    print(f"corner dt={dt:.0e}: k={x[0]:.6f} tau={x[1]:.6f} N={x[2]:.4f}")

k,tau,N = fsolve(res3, [0.61923,0.35326,21.35], args=(5e-5,), xtol=1e-11, epsfcn=1e-8)
t,y,v = simulate(k,tau,N,dt=2.5e-5,Tend=9.0)
def wmin(lo,hi):
    m=(t>=lo)&(t<=hi); j=np.argmin(v[m]); return t[m][j], v[m][j]
t0,w0=wmin(2.5,3.16); t3,w3=wmin(3.18,3.99); t4,w4=wmin(4.02,4.99)
a2=np.where(np.abs(y-1)>0.02)[0]; a5=np.where(np.abs(y-1)>0.05)[0]; a1=np.where(np.abs(y-1)>0.01)[0]
print(f"corner contacts: t0={t0:.4f} t1={t3:.4f} t2={t4:.4f} (residuals {w0:+.1e},{w3:+.1e},{w4:+.1e})")
print(f"corner Ts: 1%={t[a1[-1]+1]:.4f} 2%={t[a2[-1]+1]:.4f} 5%={t[a5[-1]+1]:.4f}")

def res2(x, N, dt):
    m = mins(x[0],x[1],N,dt); return [m[1], m[2]]
for dt in [1e-4, 5e-5]:
    best=None; x=np.array([0.6087,0.3533])
    for N in np.arange(15.2,16.21,0.2):
        x=fsolve(res2,x,args=(float(N),dt),xtol=1e-11,epsfcn=1e-8)
        t,y,v=simulate(x[0],x[1],float(N),dt=dt,Tend=10.0)
        a=np.where(np.abs(y-1)>0.01)[0]; ts=t[a[-1]+1]
        if best is None or ts<best[0]: best=(ts,N,x.copy())
    print(f"1% optimum dt={dt:.0e}: N={best[1]:.2f} k={best[2][0]:.5f} tau={best[2][1]:.5f} Ts1={best[0]:.5f}")
