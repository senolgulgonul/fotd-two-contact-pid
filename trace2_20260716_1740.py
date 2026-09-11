import numpy as np
from scipy.optimize import fsolve

# Reduced open loop L(s)=k(tau s+1)e^{-s}/(s(1+s/N)) = e^{-s} * [ k/s + k(tau-1/N)/(1+s/N) ]
# States: q'=e, f'=N(e-f); u = k q + k(tau-1/N) f; y(t) = u(t-1); e = 1 - y.
# This partial-fraction realization avoids differentiating the setpoint step (no delta term).
# Anchor: (k,tau,N)=(0.610489,0.352510,16.351) must give Ts(2%)=3.0985 and two grazing
# slope minima in [3,4) and [4,5)  (paper eq. 8).

def simulate(k, tau, N, dt=1e-4, Tend=9.0):
    nd = int(round(1.0/dt)); steps = int(round(Tend/dt))
    ubuf = np.zeros(steps)
    q = 0.0; f = 0.0
    y = np.zeros(steps); v = np.zeros(steps)
    c2 = k*(tau - 1.0/N)
    for i in range(steps):
        yd = ubuf[i-nd] if i >= nd else 0.0
        y[i] = yd
        e = 1.0 - yd
        u = k*q + c2*f
        ubuf[i] = u
        q += dt*e
        f += dt*N*(e - f)
    v[1:-1] = (y[2:] - y[:-2])/(2*dt); v[0]=0.0; v[-1]=v[-2]
    t = np.arange(steps)*dt
    return t, y, v

def win_min(t, v, lo, hi):
    m = (t >= lo) & (t <= hi)
    j = np.argmin(v[m]); tt = t[m]
    return tt[j], v[m][j], (j == 0 or j == m.sum()-1)

def settling(t, y, band):
    above = np.where(np.abs(y - 1.0) > band)[0]
    return t[above[-1]+1] if len(above) and above[-1]+1 < len(t) else t[-1]

# ---- Anchor validation ----
k0, tau0, N0 = 0.610489, 0.352510, 16.351
t, y, v = simulate(k0, tau0, N0, dt=5e-5)
print("ANCHOR (k*,tau*,N*):")
for lo,hi in [(2.0,3.0),(3.0,4.0),(4.0,5.0),(5.0,6.0)]:
    tm, vm, edge = win_min(t, v, lo, hi)
    print(f"  [{lo},{hi}): min v={vm:.6f} at t={tm:.4f}{' (edge)' if edge else ''}")
print(f"  Ts(1%)={settling(t,y,0.01):.4f}  Ts(2%)={settling(t,y,0.02):.4f} (paper 3.0985)  Ts(5%)={settling(t,y,0.05):.4f}")
print(f"  overshoot: {max(0.0,(y.max()-1.0)*100):.4f}%   min v (t>0.05): {v[(t>0.05)&(t<8.5)].min():.6f}")

# ---- Contact solve and curve trace (C1-C2 curve: interior contacts in [3,4) and [4,5)) ----
def residuals(x, N, dt=2e-4):
    k, tau = x
    t, y, v = simulate(k, tau, N, dt=dt, Tend=6.5)
    _, r3, _ = win_min(t, v, 3.02, 3.99)
    _, r4, _ = win_min(t, v, 4.02, 4.99)
    return [r3, r4]

Ngrid = [10,11,12,13,14,15,16,16.351,17,18,19,20,20.5,21,21.5,22,23]
sol = {}
x = np.array([k0, tau0])
for N in [n for n in Ngrid if n >= 16.0]:
    x = fsolve(residuals, x, args=(float(N),), xtol=1e-10, epsfcn=1e-8); sol[N] = x.copy()
x = np.array([k0, tau0])
for N in sorted([n for n in Ngrid if n < 16.0], reverse=True):
    x = fsolve(residuals, x, args=(float(N),), xtol=1e-10, epsfcn=1e-8); sol[N] = x.copy()

print(f"\n{'N':>7} {'k':>9} {'tau':>9} {'t3':>7} {'t4':>7} {'Ts1%':>8} {'Ts2%':>8} {'Ts5%':>8} {'w2':>9} {'w5':>9}")
for N in sorted(sol):
    k, tau = sol[N]
    t, y, v = simulate(k, tau, float(N), dt=5e-5, Tend=10.0)
    t3, r3, _ = win_min(t, v, 3.02, 3.99)
    t4, r4, _ = win_min(t, v, 4.02, 4.99)
    _, w2, e2 = win_min(t, v, 2.05, 2.99)   # NOTE: this window misses the t=3+ layer trough; see corner2.py
    _, w5, e5 = win_min(t, v, 5.02, 5.99)
    ts1, ts2, ts5 = settling(t,y,0.01), settling(t,y,0.02), settling(t,y,0.05)
    print(f"{N:>7.3f} {k:>9.6f} {tau:>9.6f} {t3:>7.3f} {t4:>7.3f} {ts1:>8.4f} {ts2:>8.4f} {ts5:>8.4f} {w2:>9.5f} {w5:>9.5f}")
