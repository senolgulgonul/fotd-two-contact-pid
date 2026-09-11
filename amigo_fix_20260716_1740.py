import numpy as np, control as ct
s=ct.tf('s')
def settle(Fs,Tend=40.0,npts=800000):
    tt=np.linspace(0,Tend,npts); to,yo=ct.step_response(Fs,T=tt)
    yf=yo[-1]; band=0.02*abs(yf); Ts=tt[-1]
    for k in range(len(yo)-1,-1,-1):
        if abs(yo[k]-yf)>band: Ts=to[min(k+1,len(to)-1)]; break
    return Ts,max(0.0,(yo.max()-yf)/yf*100),yf

# AMIGO PID (Astrom & Hagglund, Advanced PID Control 2006), FOTD K=1, L=1, T=TL.
# Setpoint weights: beta = 1 if L > T else 0; gamma = 0. Derivative filter Tf = Td/10.
def amigo(TL):
    K=1.0;L=1.0;T=TL
    Kc=(1/K)*(0.2+0.45*T/L)
    Ti=L*(0.4*L+0.8*T)/(L+0.1*T)
    Td=0.5*L*T/(0.3*L+T)
    b=1.0 if L>T else 0.0
    return Kc,Ti,Td,b

for TL in [0.5,1.0,2.0,5.0]:
    Kc,Ti,Td,b=amigo(TL); Tf=Td/10.0; g=0.0
    Cff=Kc*b+(Kc/Ti)/s+Kc*Td*s*g/(1+Tf*s)     # reference path
    Cfb=Kc+(Kc/Ti)/s+Kc*Td*s/(1+Tf*s)         # feedback path
    numd,dend=ct.pade(1.0,14); P=(1/(TL*s+1))*ct.tf(numd,dend)
    Fs=ct.minreal(P*Cff/(1+P*Cfb),verbose=False)
    Ts,OS,yf=settle(Fs)
    print(f"AMIGO T/L={TL} (beta={b:.0f}): Ts/L={Ts:.3f} OS={OS:.2f}% yf={yf:.4f}")

# CHR-0% series vs parallel, full loops, several T/L
print()
for TL in [0.5,1.0,2.0,5.0]:
    T=TL; Kc=0.6*T; Ti=T; Td=0.5
    numd,dend=ct.pade(1.0,14); P=(1/(T*s+1))*ct.tf(numd,dend)
    Cser=Kc*(1+1/(Ti*s))*(1+Td*s)          # series/interacting
    Cpar=Kc*(1+1/(Ti*s)+Td*s)              # parallel/ideal
    Ts1,OS1,_=settle(ct.feedback(Cser*P,1))
    Ts2,OS2,_=settle(ct.feedback(Cpar*P,1))
    print(f"T/L={TL}: series Ts/L={Ts1:.3f} OS={OS1:.2f}%   parallel Ts/L={Ts2:.3f} OS={OS2:.2f}%")
