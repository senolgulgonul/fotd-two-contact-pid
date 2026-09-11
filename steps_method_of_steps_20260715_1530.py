import sympy as sp

k, tau, N, s = sp.symbols('k tau N sigma', positive=True)
X = sp.Symbol('X', positive=True)   # placeholder for exp(-N*sigma) within the current interval

def hatW(E):
    """particular solution of v' + N v = E(sigma), E polynomial: sum (-1)^j E^(j) / N^(j+1)"""
    out, term, j = sp.S(0), E, 0
    while term != 0:
        out += (-1)**j * term / N**(j+1)
        term = sp.diff(term, s); j += 1
    return sp.expand(out)

def int_poly_exp(Phi):
    """int_0^sigma Phi(xi) e^{-N xi} dxi = Cphi - G(sigma) X, G = sum Phi^(j)/N^(j+1)"""
    G, term, j = sp.S(0), Phi, 0
    while term != 0:
        G += term / N**(j+1)
        term = sp.diff(term, s); j += 1
    Cphi = G.subs(s, 0)
    return sp.expand(Cphi), sp.expand(G)

def step_interval(E, Phi, I0, w0):
    """Given e(sigma) = E + Phi*X on this interval and entry states, return
    u(sigma) = Upoly + Ulayer*X on this interval and exit states (I1, w1)."""
    intE = sp.integrate(E, (s, 0, s))
    Cphi, G = int_poly_exp(Phi)
    Ipoly = sp.expand(I0 + intE + Cphi)
    Ilayer = sp.expand(-G)
    W = hatW(E)
    Psi = sp.integrate(Phi, (s, 0, s))
    wpoly = W
    wlayer = sp.expand(w0 - W.subs(s, 0) + Psi)
    Upoly  = sp.expand(k*Ipoly  + k*(tau*N-1)*wpoly)
    Ulayer = sp.expand(k*Ilayer + k*(tau*N-1)*wlayer)
    I1 = (Ipoly + Ilayer*X).subs([(s,1),(X,0)])       # drop e^{-N}
    w1 = (wpoly + wlayer*X).subs([(s,1),(X,0)])
    return Upoly, Ulayer, sp.expand(I1), sp.expand(w1)

# interval 0: e = 1, states zero
E, Phi, I0, w0 = sp.S(1), sp.S(0), sp.S(0), sp.S(0)
U = {}
for n in range(0, 4):
    Upoly, Ulayer, I0, w0 = step_interval(E, Phi, I0, w0)
    U[n] = (Upoly, Ulayer)
    # next interval's e: y_{n+1}(sigma) = u_n(sigma), e = 1 - y
    E   = sp.expand(1 - Upoly)
    Phi = sp.expand(-Ulayer)

# y' on interval [m, m+1) is d/dsigma of u_{m-1}: y'_m = Upoly' + (Ulayer' - N*Ulayer) X
def yprime(m):
    Up, Ul = U[m-1]
    return sp.expand(sp.diff(Up, s)), sp.expand(sp.diff(Ul, s) - N*Ul)

P3, L3 = yprime(3)
P4, L4 = yprime(4)
print('P3 degree:', sp.degree(P3, s), ' L3 degree:', sp.degree(L3, s))
print('P4 degree:', sp.degree(P4, s), ' L4 degree:', sp.degree(L4, s))

sp.pickle = None
import pickle
pickle.dump({'P3':sp.srepr(P3),'L3':sp.srepr(L3),'P4':sp.srepr(P4),'L4':sp.srepr(L4)}, open('exprs.pkl','wb'))

# validation at the reference optimum
import numpy as np, sys
sys.path.insert(0, '/home/claude/flatness')
from flat import make_step
kv, tv, Nv = 0.612821, 0.3522, 17.6
f3 = sp.lambdify(s, (P3 + L3*sp.exp(-N*s)).subs([(k,kv),(tau,tv),(N,Nv)]), 'numpy')
f4 = sp.lambdify(s, (P4 + L4*sp.exp(-N*s)).subs([(k,kv),(tau,tv),(N,Nv)]), 'numpy')

dt=2e-5
Ad, Bd = make_step(Nv, dt)
n=int(5.2/dt); nd=int(1/dt); c1,c2=kv*Nv,kv*Nv*tv
ub=np.zeros(n+1); vb=np.zeros(n+1); x=np.zeros(2)
for i in range(n):
    y = ub[i-nd] if i>=nd else 0.0
    e = 1.0-y
    ub[i]=c1*x[0]+c2*x[1]; vb[i]=c1*x[1]+c2*(e-Nv*x[1])
    x = Ad@x + Bd*e
v=vb[:n+1-nd]; t=1.0+np.arange(len(v))*dt
for m,f in [(3,f3),(4,f4)]:
    seg=(t>=m+0.01)&(t<m+0.99)
    err = np.max(np.abs(f(t[seg]-m) - v[seg]))
    print(f'interval [{m},{m+1}): max |symbolic - simulation| = {err:.2e}')
print(f'layer coefficient at contact 2 (sigma=0.26): L4 = {float(L4.subs([(k,kv),(tau,tv),(N,Nv),(s,0.26)])):.5f}, '
      f'layer term = {float((L4*sp.exp(-N*s)).subs([(k,kv),(tau,tv),(N,Nv),(s,0.26)])):.2e}')
