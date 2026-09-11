import numpy as np, control as ct
# User's two-contact rule reduced loop: L(s)=k(tau s+1)e^{-s}/(s(1+s/N)), k=0.6105, tau=0.3525, N=16.35
# Closed loop Y/W = L/(1+L). Verify Ts(2%)=3.0985, OS=0, and plant-invariance across T/L.
k,tau,N=0.610489,0.352510,16.351
s=ct.tf('s')
Gr=k*(tau*s+1)/(s*(1+s/N))
for npade in [10,14,18]:
    numd,dend=ct.pade(1.0,npade); D=ct.tf(numd,dend)
    L=Gr*D
    Fs=ct.feedback(L,1)
    tt=np.linspace(0,8,400000)
    to,yo=ct.step_response(Fs,T=tt)
    yf=yo[-1]; band=0.02*abs(yf); Ts=tt[-1]
    for kk in range(len(yo)-1,-1,-1):
        if abs(yo[kk]-yf)>band: Ts=to[min(kk+1,len(to)-1)]; break
    OS=max(0.0,(yo.max()-yf)/yf*100)
    dif=np.diff(yo); mono=(dif>=-1e-6).all()
    print(f"Pade-{npade}: Ts/L={Ts:.4f} (paper 3.0985) OS={OS:.3f}% mono={mono}")

# plant-invariance: simulate full plant K e^{-Ls}/(Ts+1) with un-collapsed controller for T/L in {0.5,1,5}
print("\nPlant-invariance check (full loop, un-collapsed controller):")
for TL in [0.5,1.0,5.0]:
    T=TL; L_=1.0
    Ki=0.6105; Kp=0.6105*(T+0.3525)/1.0; Kd=0.2152*T; Nf=16.35
    Cnum=[Kd,Kp,Ki]; Cden=[1/Nf,1,0]  # (Kd s^2+Kp s+Ki)/(s(1+s/N))=(...)/( (1/N)s^2+s )
    C=ct.tf(Cnum,Cden)
    numd,dend=ct.pade(1.0,14); P=(1.0/(T*s+1))*ct.tf(numd,dend)
    Fs=ct.feedback(C*P,1); tt=np.linspace(0,8,400000)
    to,yo=ct.step_response(Fs,T=tt); yf=yo[-1]; band=0.02*abs(yf); Ts=tt[-1]
    for kk in range(len(yo)-1,-1,-1):
        if abs(yo[kk]-yf)>band: Ts=to[min(kk+1,len(to)-1)]; break
    OS=max(0.0,(yo.max()-yf)/yf*100)
    print(f"  T/L={TL}: Ts/L={Ts:.4f} OS={OS:.3f}%")