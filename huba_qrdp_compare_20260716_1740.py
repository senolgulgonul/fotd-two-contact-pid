import numpy as np
from math import sqrt, exp

# Huba & Vrancic 2018 (IFAC PID'18) QRDP filtered-PID tuning for FOTD plant.
# Plant S(s) = K e^{-L s}/(T s + 1), normalize K=1, L=Td=1, a=1/T, A = a*Td = L/T.
# Controller (their eq 16): C(s) = Kc (1 + Ti s + Ti TD s^2)/(Ti s)  [ideal PID]
#   => Kp = Kc, Ki = Kc/Ti, Kd = Kc*TD
# Prefilter (their eq 16): Fp(s) = (c Ti TD s^2 + b Ti s + 1)/(Ti TD s^2 + Ti s + 1)
# Quadruple real dominant pole formulas (their eq 17,18), corrected reading:

def qrdp(A):
    S = sqrt(A*A + 12.0)
    Ko  = 0.5*(S*(A+12.0) - (A*A + 2*A + 36.0))*exp((S - A - 6.0)/2.0)   # = Kco*Ks*Td
    tio = 2.0*(36.0 + 2*A + A*A - (A+12.0)*S) / (A**3 + 12*A*A + 36*A + 288.0 - (A*A + 12*A + 84.0)*S)
    tDo = (2.0 - S)/(A*A + 2*A + 36.0 - (A+12.0)*S)
    num = A**3 + 12*A*A + 36*A + 288.0 - (A*A + 12*A + 84.0)*S
    bo  = 2.0*num / ((A*A + 2*A + 36.0 - S*(A+12.0))*(A + 6.0 - S))
    co  = num / ((S - 2.0)*(A + 6.0 - S))
    so  = -(6.0 + A - S)/2.0
    Kc, Ti, TD = Ko, tio, tDo      # Ks=Td=1
    return dict(A=A, S=S, so=so, Kc=Kc, Ti=Ti, TD=TD, b=bo, c=co,
                Kp=Kc, Ki=Kc/Ti, Kd=Kc*TD)

# ---- Validation anchor: a=0 (IPDT) must give IAEs = 2.1547, IAEi=4.7626 (Huba 2018 eq 19) ----
t0 = qrdp(0.0)
IAEs_a0 = t0['Ti']*(1 - t0['b'])          # a=0 => IAEs = Ti(1-b)
IAEi_a0 = t0['Ti']/t0['Kc']               # = 1/Ki
print("VALIDATION a=0 (IPDT):")
print(f"  Kc={t0['Kc']:.4f} Ti={t0['Ti']:.4f} TD={t0['TD']:.4f} b={t0['b']:.4f}")
print(f"  IAEs={IAEs_a0:.4f} (Huba 2.1547)   IAEi={IAEi_a0:.4f} (Huba 4.7626)")
print()

# ---- Time-domain simulation of the 2DOF loop on FOTD via method-of-steps ODE (delay buffer) ----
# States: plant y ; controller integrator qi ; derivative filter state (fast) ; prefilter states.
# Realize C(s)=Kp+Ki/s+Kd s with derivative filter 1/(1+s/N) (N large ~ Huba ideal Tf->0).
# Prefilter Fp(s)=(c Ti TD s^2 + b Ti s +1)/(Ti TD s^2 + Ti s +1): 2nd-order, realize in controllable canonical form.

def simulate(A, N=200.0, dt=2e-4, Tend=40.0):
    p = qrdp(A)
    Kp, Ki, Kd = p['Kp'], p['Ki'], p['Kd']
    Ti, TD, b, c = p['Ti'], p['TD'], p['b'], p['c']
    Tinv = A                      # 1/T = A (L=1)
    # Prefilter Fp(s)=Nf(s)/Df(s), Df = TiTD s^2 + Ti s + 1, Nf = c TiTD s^2 + b Ti s + 1
    a2, a1, a0 = Ti*TD, Ti, 1.0
    n2, n1, n0 = c*Ti*TD, b*Ti, 1.0
    # controllable canonical: state z (2), input w(step=1). x1'=x2 ; x2' = (-a0 x1 - a1 x2 + w)/a2
    # output wf = (n0 - n2*a0/a2) x1 + (n1 - n2*a1/a2) x2 + (n2/a2) w
    L = 1.0
    ndelay = int(round(L/dt))
    steps = int(round(Tend/dt))
    ubuf = np.zeros(ndelay+1)     # circular buffer for u(t-L)
    bi = 0
    y = 0.0                       # plant output
    qi = 0.0                      # integrator state of controller (integral of error into Ki/s)
    df = 0.0                      # derivative-filter state: df approximates d/dt of (Kd * ef) filtered
    z1 = z2 = 0.0                 # prefilter states
    w = 1.0
    ys = np.empty(steps)
    ts = np.empty(steps)
    # derivative filter: realize Kd*s/(1+s/N) acting on error e = wf - y.
    # state xd, xd' = -N xd + e ; output ud = Kd*N*(e*... ) -> use: G=Kd*N*s/(s+N)=Kd*N*(1 - N/(s+N))
    # implement xf' = -N xf + N e ; ud = Kd*N*(e - xf)   (since s/(s+N) e = e - N/(s+N)e)
    xf = 0.0
    IAE = 0.0
    for k in range(steps):
        # prefilter output wf
        x2dot = (-a0*z1 - a1*z2 + w)/a2
        wf = (n0 - n2*a0/a2)*z1 + (n1 - n2*a1/a2)*z2 + (n2/a2)*w
        e = wf - y
        # controller u = Kp e + Ki qi + ud
        ud = Kd*N*(e - xf)
        u = Kp*e + Ki*qi + ud
        # delayed input
        ud_delayed = ubuf[bi]
        # plant: y' = A*(u(t-L) - y)
        ydot = Tinv*(ud_delayed - y)
        # record
        ys[k] = y; ts[k] = k*dt
        IAE += abs(w - y)*dt
        # integrate (explicit Euler; dt tiny)
        y  += dt*ydot
        qi += dt*e
        xf += dt*(-N*xf + N*e)
        z1 += dt*z2
        z2 += dt*x2dot
        # push current controller output into delay buffer
        ubuf[bi] = u
        bi = (bi+1) % (ndelay+1)
    # settling time (2% band) and overshoot
    yf = ys[-1]
    band = 0.02*abs(yf)
    Ts = ts[-1]
    for k in range(steps-1, -1, -1):
        if abs(ys[k]-yf) > band:
            Ts = ts[k+1] if k+1 < steps else ts[k]
            break
    OS = max(0.0, (ys.max()-yf)/yf*100.0)
    return dict(A=A, TL=1.0/A if A>0 else np.inf, Ts=Ts, OS=OS, IAE=IAE,
                Kp=Kp, Ki=Ki, Kd=Kd, N=N, yfinal=yf)

# Validate simulation on a=0 IPDT: compare IAE to 2.1547 (approx; our sim uses filter+Euler)
print("Simulation self-check (a=0 IPDT, expect IAEs approx 2.15):")
r0 = simulate(1e-6, N=200, dt=2e-4, Tend=40)
print(f"  sim IAEs={r0['IAE']:.3f}  Ts(2%)={r0['Ts']:.3f}  OS={r0['OS']:.2f}%")
print()

print("QRDP FPID on FOTD, T/L sweep (N=200 ~ ideal Huba Tf->0):")
print(f"{'T/L':>5} {'A=L/T':>6} {'Kp':>7} {'Ki':>7} {'Kd':>7} {'Ts/L':>7} {'OS%':>6} {'IAEs':>7}")
for A in [2.0, 1.0, 0.5, 0.2]:
    r = simulate(A, N=200, dt=2e-4, Tend=40)
    print(f"{1.0/A:>5.1f} {A:>6.2f} {r['Kp']:>7.3f} {r['Ki']:>7.3f} {r['Kd']:>7.3f} {r['Ts']:>7.3f} {r['OS']:>6.2f} {r['IAE']:>7.3f}")
