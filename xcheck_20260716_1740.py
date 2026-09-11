import numpy as np, control as ct
from math import sqrt, exp

def qrdp(A):
    S=sqrt(A*A+12.0)
    Ko=0.5*(S*(A+12.0)-(A*A+2*A+36.0))*exp((S-A-6.0)/2.0)
    tio=2.0*(36.0+2*A+A*A-(A+12.0)*S)/(A**3+12*A*A+36*A+288.0-(A*A+12*A+84.0)*S)
    tDo=(2.0-S)/(A*A+2*A+36.0-(A+12.0)*S)
    num=A**3+12*A*A+36*A+288.0-(A*A+12*A+84.0)*S
    bo=2.0*num/((A*A+2*A+36.0-S*(A+12.0))*(A+6.0-S))
    co=num/((S-2.0)*(A+6.0-S))
    Kc=Ko/A; Ti=tio; TD=tDo
    return Kc,Ti,TD,bo,co

for A in [2.0,1.0,0.5,0.2]:
    Kc,Ti,TD,b,c = qrdp(A)
    T=1.0/A
    s=ct.tf('s')
    # controller ideal PID with fast filter N
    N=200.0
    C = (Kc*(1+Ti*s+Ti*TD*s**2)/(Ti*s)) * (1/(1+s/N))
    Fp = (c*Ti*TD*s**2+b*Ti*s+1)/(Ti*TD*s**2+Ti*s+1)
    # plant with high-order Pade delay
    npade=40
    numd,dend = ct.pade(1.0,npade)
    Pdelay=ct.tf(numd,dend)
    P = (1/(T*s+1))*Pdelay
    Fcl = ct.feedback(C*P,1)
    Fs = ct.minreal(Fp*Fcl, verbose=False)
    tt=np.linspace(0,60,600000)
    tout,yout=ct.step_response(Fs,T=tt)
    yf=yout[-1]
    OS=max(0.0,(yout.max()-yf)/yf*100)
    band=0.02*abs(yf); Ts=tt[-1]
    for k in range(len(yout)-1,-1,-1):
        if abs(yout[k]-yf)>band:
            Ts=tout[k+1] if k+1<len(tout) else tout[k]; break
    # monotonicity: any decrease?
    dif=np.diff(yout); mono = (dif>=-1e-6).all()
    print(f"T/L={T:>4.1f} A={A:>4.2f}  Kp={Kc:.3f} Ki={Kc/Ti:.3f} Kd={Kc*TD:.3f}  Ts/L={Ts:6.3f} OS={OS:5.2f}% mono={mono}")