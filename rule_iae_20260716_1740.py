import numpy as np, control as ct
s=ct.tf('s')
# Rule setpoint IAEs, whole-filter controller, plant-invariant so use reduced loop
k,tau,N=0.610489,0.352510,16.351
Gr=k*(tau*s+1)/(s*(1+s/N))
numd,dend=ct.pade(1.0,16); L=Gr*ct.tf(numd,dend)
Fs=ct.feedback(L,1)
tt=np.linspace(0,20,400000); to,yo=ct.step_response(Fs,T=tt)
IAEs=np.trapz(np.abs(1-yo),tt)
print(f"Two-contact rule setpoint IAEs = {IAEs:.3f} L  (Huba QRDP: 1.60-2.05 L)")
# also confirm Ts once more high-res
yf=yo[-1];band=0.02;Ts=tt[-1]
for kk in range(len(yo)-1,-1,-1):
    if abs(yo[kk]-yf)>band: Ts=to[min(kk+1,len(to)-1)];break
print(f"rule Ts/L={Ts:.4f}")