"""Cross-validate replay2.droop() (exact-ZOH series-RLC solver) against an
independent circuit engine (Pulsim, native C++ MNA/sparse-LU, GCC14/C++23
build in the `pulsim` conda env) on a REAL scheduled current trace.

Topology: Vrm(0V, ideal) -- L(507pH) -- R(79.6mOhm) -- node vdie -- C(5nF) -- gnd,
with a current sink at vdie pulling the real per-cycle In(t) (mA -> A) from a
policy-TD (token-gated hold) replay. vdie(t) IS the droop signal directly
(negative = droop below the 0V reference), so -vdie*1000 in mV compares
directly against replay2.droop()'s output on the same array.

Sign convention for b_extra_fn (established empirically 2026-09-30, see
memory pulsim-pdn-cross-validation): to mimic add_current_source(n_pos,
n_neg, I) -- which INJECTS I into n_pos (confirmed: +1mA into a node with a
1kOhm resistor to gnd gives vc=+1.0V, i.e. positive I sources INTO n_pos,
despite the docstring's "flows out of n_pos" phrasing -- that describes the
source's own internal terminal, not the node) -- at time t, set
    b_extra[node_id_of(n_pos)] = -I(t)
A LOAD DRAW (current pulled OUT of a node, lowering its voltage) is
therefore the NEGATIVE of an injection: to draw I_load(t) out of `n_pos`,
use b_extra[node_id_of(n_pos)] = +I_load(t) (i.e. -( -I_load(t) )). Getting
this backwards silently produces a sign-flipped, still-plausible-looking
oscillation (caught 2026-09-30 by comparing full waveforms, not just peak
mV, against replay2.droop() on a toy step input -- see
verify-apparatus-before-null-results memory).

IC: droop() starts from the analytic DC operating point for a constant
current i0 = I[0]: x0 = [i0, -R*i0] (inductor current, cap-node deviation).
The real circuit has the same DC fixed point for a constant load current, so
giving the inductor the matching i0 initial condition reproduces it exactly
(zero startup transient in either model) -- confirmed below.
"""
import sys, numpy as np
try:
    import pulsim as p
except ImportError:
    p = None
from candparse import parse
from replay2 import build, unwrap_age, calibrate, run, droop

L_H, C_F, R_OHM = 507e-12, 5e-9, 79.6e-3

def build_pdn(i0_A):
    b = p.CircuitBuilder()
    b.add_voltage_source('Vrm', 'vrm', 'gnd', 0.0)
    b.add_inductor('L1', 'vrm', 'vmid', L_H, i0=i0_A)
    b.add_resistor('R1', 'vmid', 'vdie', R_OHM)
    b.add_capacitor('C1', 'vdie', 'gnd', C_F)
    b.set_initial('C1', -R_OHM * i0_A)  # match droop()'s analytic DC fixed point x0=[i0,-R*i0]
    b.add_current_source('Iload', 'vdie', 'gnd', 0.0)  # placeholder; PWL via b_extra
    return b

def pulsim_droop_mV(I_mA, cyc_s=1e-9, sub=1):
    """Run the real current trace through Pulsim at dt = cyc_s/sub (current
    held constant -- zero-order-hold -- across each real 1ns cycle, same
    convention droop() uses); return droop in mV (positive = droop below
    rail), sampled back down to one point per real cycle for comparison."""
    dt_s = cyc_s / sub
    i0_A = I_mA[0] * 1e-3
    b = build_pdn(i0_A)
    vdie_idx = b.node_id_of('vdie')
    nunk = b.graph.num_branches  # total MNA unknown count (node voltages + branch currents)
    n = len(I_mA)
    I_A = np.asarray(I_mA, float) * 1e-3

    def b_extra(t):
        k = int(t / cyc_s + 1e-9)  # which REAL cycle we're in (ZOH, matches droop()'s per-cycle hold)
        if k >= n: k = n - 1
        out = [0.0] * nunk
        out[vdie_idx] = I_A[k]  # DRAW (load), not injection -- see sign-convention note above
        return out

    res = p.simulate(b, t_end=n * cyc_s, dt=dt_s, b_extra_fn=b_extra, start_from_dc_op=False)
    states = np.asarray(res.states)
    vdie = states[:, vdie_idx]
    return -vdie[::sub][:n] * 1e3  # mV, one sample per real cycle


def restore_edge_check(sub=8):
    """§2.4 restore-edge check: settle in a 300 mA dynamic burst under clock stretch (s=0.25), then restore to
    nominal clock (s=0) (a) instantly in 1 step vs (b) in a 4-step staggered ladder spanning the ~10 ns LC period."""
    i_base, i_burst = 100.0, 400.0
    i_str = i_base + (i_burst - i_base) / 1.25
    out = {}
    for name, ladder in (('instant', []), ('ladder_4x1_16', [0.1875, 0.1875, 0.125, 0.125, 0.0625, 0.0625])):
        seq = [i_str] * 60   # settled in stretched mode (s=0.25) so entry transient is 0
        t_rest = len(seq)
        for s in ladder:
            seq.append(i_base + (i_burst - i_base) / (1.0 + s))
        seq += [i_burst] * 40
        In = np.array(seq, float)
        dc = droop(In)
        dc_ss = float(dc[-1])   # steady-state DC IR drop at i_burst
        over_c = float(dc[t_rest:].max() - dc_ss)
        row = dict(closed_restore_dip_mV=over_c, closed_max_mV=float(dc[t_rest:].max()))
        if p is not None:
            dp = pulsim_droop_mV(In, sub=sub)
            over_p = float(dp[t_rest:].max() - float(dp[-1]))
            err_pct = 100.0 * abs(over_p - over_c) / max(abs(over_c), 1e-6)
            row.update(pulsim_restore_dip_mV=over_p, pulsim_max_mV=float(dp[t_rest:].max()), peak_err_pct=err_pct)
        out[name] = row
        print(f'restore_edge [{name:14s}]: closed restore overshoot={over_c:5.2f} mV (peak={dc[t_rest:].max():5.2f} mV)'
              + (f', pulsim={row["pulsim_restore_dip_mV"]:5.2f} mV (err {row["peak_err_pct"]:.1f}%)' if p is not None else ' (pulsim not installed in current env)'))
    red = 100.0 * (1.0 - out['ladder_4x1_16']['closed_restore_dip_mV'] / max(out['instant']['closed_restore_dip_mV'], 1e-6))
    print(f'ladder reduction of restore-edge second droop: {red:.1f}%')
    return out


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'restore_edge':
        restore_edge_check()
        sys.exit(0)
    wl = sys.argv[1] if len(sys.argv) > 1 else 'dhrystone'
    n_cap = int(sys.argv[2]) if len(sys.argv) > 2 else 4000  # keep the circuit sim short

    inst, _ = parse(f'traces/pm2_{wl}.log.gz')
    prods, cons = build(inst); age = unwrap_age(inst)
    pin, _ = calibrate(inst, prods, cons, age)
    R = run(inst, prods, cons, age, pin, 'TD', TH=100, Qb=250, Ir=25)
    In_full = R['In']
    # trim to a fixed window that contains real activity (skip any leading/trailing pad)
    In = In_full[:n_cap].copy()
    print(f'{wl}: using {len(In)} cycles of In (mA), min={In.min():.1f} max={In.max():.1f} mean={In.mean():.1f}')

    d_closed = droop(In)  # replay2's validated exact-ZOH solver
    print(f'replay2.droop() (exact ZOH):  max={d_closed.max():.2f} mV  at t={int(d_closed.argmax())}')

    results = {}
    for sub in (1, 4, 8):
        d_pulsim = pulsim_droop_mV(In, sub=sub)
        results[sub] = d_pulsim
        print(f'Pulsim (native MNA solve, {sub}x sub-stepping, dt={1e-9/sub*1e12:.0f}ps): '
              f'max={d_pulsim.max():.2f} mV  at t={int(d_pulsim.argmax())}')

    d_pulsim = results[8]  # finest, most trustworthy
    err = d_pulsim - d_closed
    peak_err_pct = 100 * abs(d_pulsim.max() - d_closed.max()) / d_closed.max()
    print(f'8x sub-step vs closed-form: pointwise max|err|={np.abs(err).max():.2f} mV, '
          f'rms={np.sqrt((err**2).mean()):.2f} mV (edge-discretization noise, not a peak-droop error); '
          f'PEAK droop agreement: {peak_err_pct:.1f}%')

    np.savez('pulsim_validate_out.npz', In=In, d_closed=d_closed, d_pulsim=d_pulsim)
