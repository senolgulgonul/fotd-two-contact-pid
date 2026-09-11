import numpy as np

k, tau, N = 0.613, 0.352, 17.4

def Delta(s):  return s*(1+s/N) + k*(tau*s+1)*np.exp(-s)
def dDelta(s): return 1 + 2*s/N + k*np.exp(-s)*(tau - tau*s - 1)

def newton(s0, iters=60):
    s = complex(s0)
    for _ in range(iters):
        d = Delta(s)/dDelta(s)
        s -= d
        if abs(d) < 1e-14: break
    return s

# scan initial guesses over the dominant strip
roots = []
for re0 in np.arange(-8, 0.1, 0.5):
    for im0 in np.arange(0, 40.1, 0.7):
        r = newton(re0 + 1j*im0)
        if abs(Delta(r)) < 1e-10 and r.real > -12 and -0.5 < r.imag < 45:
            if not any(abs(r - q) < 1e-6 for q in roots):
                roots.append(r)
roots = sorted(roots, key=lambda z: -z.real)

# residues of y(t) = 1 + sum c_j e^{p_j (t-1)}: c_j = k(tau p+1)/(p Delta'(p))  [e^{-p} folded into shift]
print(f"{'p (pole)':>28} {'|c| (residue of y)':>18} {'mode of y_prime |p c|':>20}")
for p in roots[:8]:
    c = k*(tau*p+1)/(p*dDelta(p))
    print(f"{p.real:12.6f} {p.imag:+12.6f}j {abs(c):18.6f} {abs(p*c):20.6f}")

# hypothesis test: a vs b of dominant pair
reals = [p for p in roots if abs(p.imag) < 1e-8]
pairs = [p for p in roots if p.imag > 1e-8]
a = -max(p.real for p in reals) if reals else None
b1 = -max(p.real for p in pairs); w1 = [p for p in pairs if -p.real==b1][0].imag
print(f"\nreal pole a = {a:.6f}" if a else "\nno real pole in strip")
print(f"dominant pair b = {b1:.6f}, omega = {w1:.6f}  (echo spacing 2pi/omega = {2*np.pi/w1:.4f})")
if a: print(f"decay-degeneracy test: a - b = {a-b1:+.6f}  (relative {(a-b1)/a:+.3%})")
