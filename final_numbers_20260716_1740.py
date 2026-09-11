import numpy as np
from math import sqrt, exp
def qrdp(A):
    S=sqrt(A*A+12.0)
    Ko=0.5*(S*(A+12.0)-(A*A+2*A+36.0))*exp((S-A-6.0)/2.0)
    tio=2.0*(36.0+2*A+A*A-(A+12.0)*S)/(A**3+12*A*A+36*A+288.0-(A*A+12*A+84.0)*S)
    tDo=(2.0-S)/(A*A+2*A+36.0-(A+12.0)*S)
    so=-(6.0+A-S)/2.0
    Kc=Ko/A;Ti=tio;TD=tDo;c=1.0/(Ti*TD*so*so);b=-2.0/(Ti*so)
    return Kc,Ti,TD,b,c,so
def euler(A,N,dt,Tend,rule=None):
    if rule is None:
        Kc,Ti,TD,b,c,so=qrdp(A);Tinv=A
        a2,a1,a0=Ti*TD,Ti,1.0;n2,n1,n0=c*Ti*TD,b*Ti,1.0
        Kp,Ki,Kd=Kc,Kc/Ti,Kc*TD
    else:
        Kp,Ki,Kd,N=rule;Tinv=A
        a2,a1,a0=1,0,0;n2,n1,n0=0,0,1  # prefilter=1 (1DOF)
    nd=int(round(1.0/dt));steps=int(round(Tend/dt))
    ubuf=np.zeros(nd+1);bi=0;y=qi=xf=z1=z2=0.0;w=1.0
    ys=np.empty(steps);IAE=0.0
    for k in range(steps):
        if rule is None:
            x2dot=(-a0*z1-a1*z2+w)/a2
            wf=(n0-n2*a0/a2)*z1+(n1-n2*a1/a2)*z2+(n2/a2)*w
        else:
            wf=w;x2dot=0.0
        e=wf-y;ud=Kd*N*(e-xf) if rule else Kd*N*(e-xf)
        u=Kp*e+Ki*qi+ud
        ydot=Tinv*(ubuf[bi]-y)
        ys[k]=y;IAE+=abs(w-y)*dt
        y+=dt*ydot;qi+=dt*e;xf+=dt*(-N*xf+N*e)
        if rule is None: z1+=dt*z2;z2+=dt*x2dot
        ubuf[bi]=u;bi=(bi+1)%(nd+1)
    yf=ys[-1];band=0.02*abs(yf);Ts=(steps-1)*dt
    for k in range(steps-1,-1,-1):
        if abs(ys[k]-yf)>band: Ts=min(k+1,steps-1)*dt;break
    OS=max(0.0,(ys.max()-yf)/yf*100)
    return Ts,OS,IAE

print("HUBA QRDP FPID (2DOF, ideal):")
print(f"{'T/L':>5}{'Ts/L':>8}{'OS%':>7}{'IAEs':>8}")
for A in [2.0,1.0/0.7,1.0,0.5,0.2]:
    Ts,OS,IAE=euler(A,200,1e-4,60)
    print(f"{1/A:>5.2f}{Ts:>8.3f}{OS:>7.2f}{IAE:>8.3f}")

print("\nTWO-CONTACT RULE (1DOF, N=16.35):")
print(f"{'T/L':>5}{'Ts/L':>8}{'OS%':>7}{'IAEs':>8}")
for TL in [0.5,0.7,1.0,2.0,5.0]:
    T=TL;Ki=0.6105;Kp=0.6105*(T+0.3525);Kd=0.2152*T;N=16.35
    Ts,OS,IAE=euler(1.0/T if T>0 else 1,None,1e-4,60,rule=(Kp,Ki,Kd,N)) if False else euler(1.0/T,None,1e-4,60,rule=(Kp,Ki,Kd,N))
    print(f"{TL:>5.1f}{Ts:>8.3f}{OS:>7.2f}{IAE:>8.3f}")