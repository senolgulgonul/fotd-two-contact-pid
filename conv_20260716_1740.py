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
    return Ko/A,tio,tDo,bo,co,-(6.0+A-S)/2.0

# prefilter double-root check at A=1
Kc,Ti,TD,b,c,so=qrdp(1.0)
# numerator c*Ti*TD s^2 + b*Ti s + 1 ; roots:
coef=[c*Ti*TD, b*Ti, 1.0]
roots=np.roots(coef)
print(f"A=1: so={so:.5f}  prefilter-num roots={roots}  (should be near so, double)")
print()

def euler(A,N,dt,Tend):
    Kc,Ti,TD,b,c,so=qrdp(A); Tinv=A
    a2,a1,a0=Ti*TD,Ti,1.0; n2,n1,n0=c*Ti*TD,b*Ti,1.0
    nd=int(round(1.0/dt)); steps=int(round(Tend/dt))
    ubuf=np.zeros(nd+1); bi=0
    y=qi=xf=z1=z2=0.0; w=1.0
    ys=np.empty(steps); tarr=np.empty(steps); IAE=0.0
    for k in range(steps):
        x2dot=(-a0*z1-a1*z2+w)/a2
        wf=(n0-n2*a0/a2)*z1+(n1-n2*a1/a2)*z2+(n2/a2)*w
        e=wf-y
        ud=Kc*TD*N*(e-xf)  # Kd=Kc*TD
        u=Kc*e+ (Kc/Ti)*qi + ud
        ydot=Tinv*(ubuf[bi]-y)
        ys[k]=y; tarr[k]=k*dt; IAE+=abs(w-y)*dt
        y+=dt*ydot; qi+=dt*e; xf+=dt*(-N*xf+N*e); z1+=dt*z2; z2+=dt*x2dot
        ubuf[bi]=u; bi=(bi+1)%(nd+1)
    yf=ys[-1]; band=0.02*abs(yf); Ts=tarr[-1]
    for k in range(steps-1,-1,-1):
        if abs(ys[k]-yf)>band: Ts=tarr[min(k+1,steps-1)]; break
    OS=max(0.0,(ys.max()-yf)/yf*100)
    return Ts,OS,IAE,yf

for dt in [2e-4,1e-4,5e-5]:
    Ts,OS,IAE,yf=euler(1.0,200,dt,60)
    print(f"Euler dt={dt:.0e}: Ts/L={Ts:.3f} OS={OS:.2f}% IAE={IAE:.3f} yf={yf:.4f}")

# stable low-order Pade lsim
for npade in [8,12,16]:
    Kc,Ti,TD,b,c,so=qrdp(1.0); T=1.0; N=200.0; s=ct.tf('s')
    C=(Kc*(1+Ti*s+Ti*TD*s**2)/(Ti*s))*(1/(1+s/N))
    Fp=(c*Ti*TD*s**2+b*Ti*s+1)/(Ti*TD*s**2+Ti*s+1)
    numd,dend=ct.pade(1.0,npade); P=(1/(T*s+1))*ct.tf(numd,dend)
    Fcl=ct.feedback(C*P,1); Fs=Fp*Fcl
    tt=np.linspace(0,60,240000)
    try:
        to,yo=ct.step_response(Fs,T=tt)
        yf=yo[-1]; band=0.02*abs(yf); Ts=tt[-1]
        for k in range(len(yo)-1,-1,-1):
            if abs(yo[k]-yf)>band: Ts=to[min(k+1,len(to)-1)]; break
        OS=max(0.0,(yo.max()-yf)/yf*100)
        print(f"Pade-{npade}: Ts/L={Ts:.3f} OS={OS:.2f}% yf={yf:.4f}")
    except Exception as ex:
        print(f"Pade-{npade}: FAIL {ex}")