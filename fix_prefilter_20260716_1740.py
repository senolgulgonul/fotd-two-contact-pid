import numpy as np, control as ct
from math import sqrt, exp

def qrdp(A):
    S=sqrt(A*A+12.0)
    Ko=0.5*(S*(A+12.0)-(A*A+2*A+36.0))*exp((S-A-6.0)/2.0)
    tio=2.0*(36.0+2*A+A*A-(A+12.0)*S)/(A**3+12*A*A+36*A+288.0-(A*A+12*A+84.0)*S)
    tDo=(2.0-S)/(A*A+2*A+36.0-(A+12.0)*S)
    so=-(6.0+A-S)/2.0
    Kc=Ko/A; Ti=tio; TD=tDo
    # prefilter cancels DOUBLE dominant pole at so: Nf(s)=c*Ti*TD*(s-so)^2
    c = 1.0/(Ti*TD*so*so)
    b = -2.0/(Ti*so)
    return Kc,Ti,TD,b,c,so

# check double root now
for A in [1.0]:
    Kc,Ti,TD,b,c,so=qrdp(A)
    coef=[c*Ti*TD, b*Ti, 1.0]
    print(f"A={A}: so={so:.5f} derived b={b:.4f} c={c:.4f} prefilter roots={np.roots(coef)}")
    # compare to eq18 b (validated at a=0 for IAEs)
print()

def euler(A,N,dt,Tend):
    Kc,Ti,TD,b,c,so=qrdp(A); Tinv=A
    a2,a1,a0=Ti*TD,Ti,1.0; n2,n1,n0=c*Ti*TD,b*Ti,1.0
    nd=int(round(1.0/dt)); steps=int(round(Tend/dt))
    ubuf=np.zeros(nd+1); bi=0
    y=qi=xf=z1=z2=0.0; w=1.0
    ys=np.empty(steps); tarr=np.empty(steps); IAE=0.0; ymono=True; ymin_after_peak=None
    for k in range(steps):
        x2dot=(-a0*z1-a1*z2+w)/a2
        wf=(n0-n2*a0/a2)*z1+(n1-n2*a1/a2)*z2+(n2/a2)*w
        e=wf-y
        ud=Kc*TD*N*(e-xf)
        u=Kc*e+(Kc/Ti)*qi+ud
        ydot=Tinv*(ubuf[bi]-y)
        ys[k]=y; tarr[k]=k*dt; IAE+=abs(w-y)*dt
        y+=dt*ydot; qi+=dt*e; xf+=dt*(-N*xf+N*e); z1+=dt*z2; z2+=dt*x2dot
        ubuf[bi]=u; bi=(bi+1)%(nd+1)
    yf=ys[-1]; band=0.02*abs(yf); Ts=tarr[-1]
    for k in range(steps-1,-1,-1):
        if abs(ys[k]-yf)>band: Ts=tarr[min(k+1,steps-1)]; break
    OS=max(0.0,(ys.max()-yf)/yf*100)
    # monotonic check on output
    dif=np.diff(ys); mono=(dif>=-1e-5).all()
    return Ts,OS,IAE,yf,mono

print("With prefilter cancelling double pole at so:")
print(f"{'T/L':>5} {'Kp':>7} {'Ki':>7} {'Kd':>7} {'Ts/L':>7} {'OS%':>6} {'mono':>6} {'IAEs':>7}")
res={}
for A in [2.0,1.0,0.5,0.2]:
    Kc,Ti,TD,b,c,so=qrdp(A)
    Ts,OS,IAE,yf,mono=euler(A,200,1e-4,60)
    res[1.0/A]=(Ts,OS)
    print(f"{1.0/A:>5.1f} {Kc:>7.3f} {Kc/Ti:>7.3f} {Kc*TD:>7.3f} {Ts:>7.3f} {OS:>6.2f} {str(mono):>6} {IAE:>7.3f}")

# cross-check T/L=1 with Pade
Kc,Ti,TD,b,c,so=qrdp(1.0); T=1.0; N=200.0; s=ct.tf('s')
C=(Kc*(1+Ti*s+Ti*TD*s**2)/(Ti*s))*(1/(1+s/N))
Fp=(c*Ti*TD*s**2+b*Ti*s+1)/(Ti*TD*s**2+Ti*s+1)
numd,dend=ct.pade(1.0,12); P=(1/(T*s+1))*ct.tf(numd,dend)
Fs=Fp*ct.feedback(C*P,1); tt=np.linspace(0,60,240000)
to,yo=ct.step_response(Fs,T=tt); yf=yo[-1]; band=0.02; Ts=tt[-1]
for k in range(len(yo)-1,-1,-1):
    if abs(yo[k]-yf)>band: Ts=to[min(k+1,len(to)-1)]; break
print(f"\nPade-12 xcheck T/L=1: Ts/L={Ts:.3f} OS={max(0,(yo.max()-yf)/yf*100):.2f}%")