import numpy as np
from scipy.optimize import fsolve
exec(open('trace2_20260716_1740.py').read().split('# ---- Anchor')[0])

# Three-contact corner. The third contact is the boundary-layer trough of the interval-3
# expression at t = 3 + O(1/N); it falls between naive windows [.,2.99] and [3.02,.], so the
# trough window here spans the interval boundary: [2.5, 3.16].
def mins(k, tau, N, dt=5e-5):
    t, y, v = simulate(k, tau, N, dt=dt, Tend=6.5)
    trough = v[(t>=2.5)&(t<=3.16)].min()          # layer trough at t=3+
    c3 = v[(t>=3.18)&(t<=3.99)].min()             # interior interval-3 contact (~3.7)
    c4 = v[(t>=4.02)&(t<=4.99)].min()             # interval-4 contact
    return trough, c3, c4

def res2(x, N, dt=5e-5):
    tr, c3, c4 = mins(x[0], x[1], N, dt); return [c3, c4]

x = np.array([0.6172, 0.35255])
print(f"{'N':>6} {'k':>9} {'tau':>9} {'trough':>10} {'c3':>10} {'c4':>10}")
sols={}
for N in [20.0,20.5,21.0,21.2,21.4,21.6,21.8,22.0]:
    x = fsolve(res2, x, args=(N,), xtol=1e-11, epsfcn=1e-8)
    tr, c3, c4 = mins(x[0], x[1], N); sols[N]=x.copy()
    print(f"{N:>6.2f} {x[0]:>9.6f} {x[1]:>9.6f} {tr:>10.6f} {c3:>10.2e} {c4:>10.2e}")

def res3(x, dt=5e-5):
    tr, c3, c4 = mins(x[0], x[1], x[2], dt); return [tr, c3, c4]
Ns=sorted(sols); trs=[mins(*sols[N],N)[0] for N in Ns]; x0=None
for i in range(len(Ns)-1):
    if trs[i]>0>=trs[i+1]:
        w=trs[i]/(trs[i]-trs[i+1]); N0=Ns[i]+(Ns[i+1]-Ns[i])*w
        x0=np.array([*(sols[Ns[i]]*(1-w)+sols[Ns[i+1]]*w), N0]); break
if x0 is None: x0=np.array([0.619,0.3533,21.3])
x = fsolve(res3, x0, xtol=1e-11, epsfcn=1e-8)
tr,c3,c4 = mins(*x)
t,y,v = simulate(*x, dt=2.5e-5, Tend=9.0)
m=(t>=2.5)&(t<=3.16); j=np.argmin(v[m]); ttr=t[m][j]
a2=np.where(np.abs(y-1.0)>0.02)[0]; a5=np.where(np.abs(y-1.0)>0.05)[0]
print(f"\nSTRICT CORNER: k={x[0]:.6f} tau={x[1]:.6f} N={x[2]:.4f}")
print(f"  residuals: trough={tr:+.2e} c3={c3:+.2e} c4={c4:+.2e}; trough at t={ttr:.5f}")
print(f"  Ts(2%)={t[a2[-1]+1]:.4f}  Ts(5%)={t[a5[-1]+1]:.4f}  OS={max(0,(y.max()-1)*100):.4f}%  minv={v[(t>0.05)&(t<8.5)].min():+.2e}")
