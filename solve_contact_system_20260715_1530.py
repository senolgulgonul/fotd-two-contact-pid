import numpy as np, pickle, sympy as sp
from scipy.optimize import fsolve, brentq

d = pickle.load(open('exprs.pkl','rb'))
k, tau, N, s = sp.symbols('k tau N sigma', positive=True)
P3, L3 = sp.sympify(d['P3']), sp.sympify(d['L3'])
P4, L4 = sp.sympify(d['P4']), sp.sympify(d['L4'])

# contact 1: double root of quadratic P3 = A s^2 + B s + C
A = sp.expand(P3).coeff(s, 2); B = sp.expand(P3).coeff(s, 1); C = sp.expand(P3).coeff(s, 0)
disc = sp.expand(B**2 - 4*A*C)
sig1 = -B/(2*A)
print('P3 coefficients (quadratic in sigma):')
for nm, ex in [('A',A),('B',B),('C',C)]:
    print(f'  {nm} =', sp.simplify(ex))

Q  = P4 + L4*sp.exp(-N*s)
Qp = sp.diff(Q, s)

f_disc = sp.lambdify((k,tau,N), disc, 'numpy')
f_sig1 = sp.lambdify((k,tau,N), sig1, 'numpy')
f_Q    = sp.lambdify((k,tau,N,s), Q, 'numpy')
f_Qp   = sp.lambdify((k,tau,N,s), Qp, 'numpy')

# y on [3,4): u_2 = Upoly2 + Ulayer2 * X. Rebuild from stored derivative? Integrate P3+L3:
# y3(sigma) = y3(0) + int_0^sigma y'. y3(0) = y at t=3 = u_2(0): get from P-chain instead:
# easier: y3' = P3 + L3 e^{-Ns}; integrate symbolically with constant y3(0) = c30(k,tau,N).
# y at t=3 equals u_2 evaluated at sigma=1 of interval 2... we did not store U; recover via
# continuity: y is continuous, y3(0) = y_2(1)|_{X->0}. Rebuild y_2' similarly? Simpler:
# integrate the ODE chain again quickly: reuse steps machinery.
exec(open('steps.py').read().split('# validation')[0])  # rebuilds U dict silently
Up2, Ul2 = U[2]
y3 = sp.expand(Up2 + Ul2*sp.exp(-N*s))     # y on [3,4), sigma local
f_y3 = sp.lambdify((k,tau,N,s), y3, 'numpy')

def curve_point(Nv, x0=(0.612, 0.352, 0.26)):
    F = lambda x: [f_disc(x[0], x[1], Nv),
                   f_Q(x[0], x[1], Nv, x[2]),
                   f_Qp(x[0], x[1], Nv, x[2])]
    sol, info, ier, msg = fsolve(F, x0, full_output=True)
    kv, tv, s2 = sol
    s1 = f_sig1(kv, tv, Nv)
    Ts = 3.0 + brentq(lambda x: f_y3(kv, tv, Nv, x) - 0.98, 0.0, 0.9)
    return kv, tv, s2, s1, Ts, ier

print(f"\n{'N':>6} {'k':>10} {'tau':>9} {'t1=3+s1':>9} {'t2=4+s2':>9} {'Ts':>10}")
x0 = (0.612, 0.352, 0.26)
for Nv in [12, 15, 17.4, 17.6, 20, 22, 25]:
    kv, tv, s2, s1, Ts, ier = curve_point(Nv, x0)
    x0 = (kv, tv, s2)
    print(f"{Nv:6.1f} {kv:10.6f} {tv:9.5f} {3+s1:9.4f} {4+s2:9.4f} {Ts:10.6f}  {'ok' if ier==1 else 'FAIL'}")
print("\ntracker reference (dt=2.5e-4, biased ~+0.0016 in Ts): N=17.6 -> k=0.612821 tau=0.3522 Ts=3.10086")
