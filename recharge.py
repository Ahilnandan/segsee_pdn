import numpy as np
VDD=0.75; L=507e-12; R=79.6e-3
# ---- Bank sized from the earlier result: 100 mA step tapered over 5 ns -> q = 250 pC
Cb=0.5e-9; Vfull=1.45; Vempty=0.95; q=Cb*(Vfull-Vempty)
print(f"bank: Cb={Cb*1e9:.1f} nF, swing {Vempty}->{Vfull} V, charge/event q={q*1e12:.0f} pC\n")
print("(a) recharge from the 1.8 V FIVR input rail via a current-limited switch (separate rail):")
for Ir in [2e-3,5e-3,10e-3,25e-3]:
    print(f"    I_r = {Ir*1e3:4.0f} mA  -> recharge time {q/Ir*1e9:6.1f} ns = cycles @1 GHz")
print("\n(b) 2:1 switched-cap doubler from Vdd (no FIVR). Output tends to 2*Vdd=1.50 V;")
print("    slow-switching-limit output resistance R_out = 1/(C_fly*f_sw); bank charges with tau = R_out*Cb")
for Cfly,fsw in [(0.1e-9,500e6),(0.1e-9,250e6),(0.05e-9,250e6),(0.25e-9,500e6)]:
    Rout=1/(Cfly*fsw); tau=Rout*Cb
    t90=tau*np.log((2*VDD-Vempty)/(2*VDD-Vfull))       # Vempty -> Vfull toward 1.5 V asymptote
    i0=(2*VDD-Vempty)/Rout                               # initial pump output current
    print(f"    C_fly {Cfly*1e9:4.2f} nF @ {fsw/1e6:4.0f} MHz: R_out {Rout:5.1f} ohm, tau {tau*1e9:5.1f} ns, "
          f"recharge {t90*1e9:6.1f} ns; initial draw from Vdd ~{2*i0*1e3:4.0f} mA")
# ---- Self-inflicted droop: pump input current is a load on Vdd itself. Trapezoidal (converged) solver.
dt=0.1e-9; t=np.arange(0,300e-9,dt)
def droop(C, iload):
    A=np.array([[-R/L,-1/L],[1/C,0]]); I2=np.eye(2); M1=np.linalg.inv(I2-dt/2*A); M2=I2+dt/2*A
    s=np.array([iload[0],VDD-R*iload[0]]); v0=s[1]; vmin=v0
    for k in range(len(t)-1):
        s=M1@(M2@s+dt/2*np.array([2*VDD/L,-(iload[k]+iload[k+1])/C])); vmin=min(vmin,s[1])
    return (v0-vmin)*1e3
base=0.05; t0=50e-9
Rout=20.0; tau=Rout*Cb; i0=(2*VDD-Vempty)/Rout
hard = base + 2*i0*np.exp(-(t-t0)/tau)*(t>=t0)                     # pump just switches on
Ilim=10e-3; ramp=np.clip((t-t0)/20e-9,0,1)*(t>=t0)                 # soft start: 20 ns ramp, 10 mA cap
soft = base + np.minimum(2*i0*np.exp(-(t-t0)/tau), 2*Ilim)*ramp
print("\nself-droop on Vdd from the pump's own input current (1 nF passive left on die):")
print(f"    pump hard start ({2*i0*1e3:.0f} mA in, tau {tau*1e9:.0f} ns): {droop(1e-9,hard):5.1f} mV")
print(f"    soft start, 20 mA-in cap, 20 ns ramp:  {droop(1e-9,soft):5.1f} mV "
      f"(recharge then takes ~{q/Ilim*1e9:.0f} ns)")
print(f"    reference: the 100 mA load step it exists to cover, same 1 nF: {droop(1e-9, base+0.1*(t>=t0)):5.1f} mV")
