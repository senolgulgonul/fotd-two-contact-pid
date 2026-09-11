import numpy as np, control as ct
s=ct.tf('s')

def settle(Fs, Tend=30.0, npts=600000):
    tt=np.linspace(0,Tend,npts)
    to,yo=ct.step_response(Fs,T=tt)
    yf=yo[-1]; band=0.02*abs(yf); Ts=tt[-1]
    for k in range(len(yo)-1,-1,-1):
        if abs(yo[k]-yf)>band: Ts=to[min(k+1,len(to)-1)]; break
    OS=max(0.0,(yo.max()-yf)/yf*100)
    IAE=np.trapezoid(np.abs(1-yo),to)
    return Ts,OS,IAE,yf

# --- CHR 0% servo PID: Kc=0.6T/(KL), Ti=T, Td=0.5L; ideal derivative, 1DOF ---
# Ti=T cancels the plant pole exactly (SERIES form): reduced loop 0.6*(0.5 L s + 1) e^{-Ls}/(L s), L=1
for npade in [12,16]:
    numd,dend=ct.pade(1.0,npade); D=ct.tf(numd,dend)
    Lred=0.6*(0.5*s+1)/s * D
    Fs=ct.feedback(Lred,1)
    Ts,OS,IAE,yf=settle(Fs)
    print(f"CHR-0% PID (reduced, Pade-{npade}): Ts/L={Ts:.4f} OS={OS:.3f}% IAE={IAE:.3f} yf={yf:.4f}")

# --- AMIGO PID, 2DOF, derivative filter Tf=Td/10 (beta rule corrected in amigo_fix.py) ---
def amigo(TL):
    K=1.0; L=1.0; T=TL
    Kc=(1/K)*(0.2+0.45*T/L)
    Ti=L*(0.4*L+0.8*T)/(L+0.1*T)
    Td=0.5*L*T/(0.3*L+T)
    return Kc,Ti,Td
