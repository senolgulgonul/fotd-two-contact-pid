import numpy as np
from spectrum import Delta, dDelta, newton, k, tau, N
import sys
sys.path.insert(0, '/home/claude/flatness')
from flat import make_step

# collect roots again (quiet)
roots = []
for re0 in np.arange(-8, 0.1, 0.5):
    for im0 in np.arange(0, 80.1, 0.7):
        r = newton(re0 + 1j*im0)
        if abs(Delta(r)) < 1e-10 and r.real > -12 and -0.5 < r.imag < 90:
            if not any(abs(r - q) < 1e-6 for q in roots):
                roots.append(r)
roots = sorted(roots, key=lambda z: -z.real)
res = {p: k*(tau*p+1)/(p*dDelta(p)) for p in roots}

def yprime_modes(t, mode_set):
    """y'(t) = sum p c e^{p(t-1)} + conj for complex."""
    v = np.zeros_like(t)
    for p in mode_set:
        c = res[p]
        term = (p*c*np.exp(p*(t-1)))
        v += 2*term.real if p.imag > 1e-8 else term.real
    return v

# exact simulation for ground truth
Ad, Bd = make_step(N, 2.5e-4)
dt=2.5e-4; n=int(10/dt); nd=int(1/dt)
c1,c2 = k*N, k*N*tau
ub=np.zeros(n+1); vb=np.zeros(n+1); x=np.zeros(2)
for i in range(n):
    y = ub[i-nd] if i>=nd else 0.0
    e = 1.0-y
    ub[i]=c1*x[0]+c2*x[1]; vb[i]=c1*x[1]+c2*(e-N*x[1])
    x = Ad@x + Bd*e
v_sim = vb[:n+1-nd]; t = 1.0+np.arange(len(v_sim))*dt

reals = [p for p in roots if abs(p.imag)<1e-8]
pairs = sorted([p for p in roots if p.imag>1e-8], key=lambda z: z.imag)
smooth = reals + pairs[:1]          # real pole + slow pair (omega=1.58)
ladder = pairs[1:]                  # the 2pi/L ladder

tt = np.array([3.75, 4.26, 5.10, 6.02, 8.09, 8.95, 7.0, 7.5])
v_sm  = yprime_modes(tt, smooth)
v_all = yprime_modes(tt, smooth+ladder)
v_lad = yprime_modes(tt, ladder)
idx = ((tt-1.0)/dt).astype(int)
print(f"{'t':>6} {'y_sim':>11} {'smooth':>11} {'ladder':>11} {'smooth+ladder':>13} {'err':>9}")
for j,tj in enumerate(tt):
    print(f"{tj:6.2f} {v_sim[idx[j]]:11.3e} {v_sm[j]:11.3e} {v_lad[j]:11.3e} {v_all[j]:13.3e} {v_all[j]-v_sim[idx[j]]:9.1e}")

# ladder pulse-train picture: min of ladder over each unit interval vs smooth floor there
print("\nper-interval: floor (smooth) at ladder-pulse minimum vs pulse depth")
for n0 in [3,4,5,6,7,8]:
    seg = (t>=n0)&(t<n0+1)
    ts = t[seg]
    lad = yprime_modes(ts, ladder); sm = yprime_modes(ts, smooth)
    i = np.argmin(lad)
    print(f"interval [{n0},{n0+1}): pulse min at t={ts[i]:.3f}, ladder={lad[i]:+.3e}, smooth={sm[i]:+.3e}, sum={lad[i]+sm[i]:+.3e}")
