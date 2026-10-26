// PDN time-domain kernel: linear network x_{k+1} = P x_k + Q (u_k + u_{k+1}) (trapezoidal on E x' = A x + B u)
// plus per-region devices updated explicitly from the previous step's state:
//   - small on-die buck (alternate parallel path, 1.8 V rail, inductor L): current-target control. A trigger (sensed
//     droop after comparator delay and/or the issue-stream predictor) sets the target to Ipk; the target then tapers
//     to 0 over t_taper so the package path takes the current back smoothly. Switches follow the target with a
//     hysteresis band and a minimum dwell (bounded switching frequency). Mode 2 = PWM-regulating buck (DVFS / FIVR).
//   - active decap (Meng/Saleh style): capacitor Ca between die node and a bottom-plate node driven by an amplifier,
//     Vbp* = Vmid - G * (Vdie - LPF(Vdie)), first-order amplifier lag tau_a, clamped to [0.05, 1.75] V. Ca and the
//     amplifier output resistance are LINEAR network elements; only the bottom-plate drive voltage is computed here.
//   - off-chip VR: drive voltage of the board inductor branch, integral loop on the mean die voltage (crossover fc_vr)
// Inputs u: [0] VR voltage, [1..R] load currents drawn (A), [R+1..2R] buck currents injected (A), [2R+1..3R] Vbp drives.
#include <cmath>
#include <vector>
#include <algorithm>
extern "C" {

struct Dev {
  int buck_mode; double Vin, L, Rs, Ipk, band, th_b; int dly_b, use_pred_b, dwell; double t_taper;
  double fsw, fc_b, Esw;
  int ad_on; double G, tau_a, f_lp, Pq; int ad_br;     // ad_br: state index of the bottom-plate drive branch current
  double Cr; int no_int;                                // mode 2: region capacitance for loop gains; P-only (parallel AC regulator)
};

int run(int ncyc, int sub, int nx, int R, const double* P, const double* Q, const int* vidx,
        const float* load_mA, const double* leak_mA, const double* shape, const unsigned char* pred,
        const Dev* dev, const double* vref_cyc, double Vnom, double fc_vr, double slew_vr, double dt, double vr_v0,
        double* x0, float* vmin_out, float* vavg_out, double* stats)
{
  const int nu = 1 + 3 * R;
  std::vector<double> x(x0, x0 + nx), xn(nx), u(nu, 0.0), un(nu, 0.0);
  std::vector<double> iL(R, 0.0), itgt(R, 0.0), integ(R, 0.0), vbp(R, 0.9), vlp(R, Vnom), vsf(R, Vnom);
  std::vector<int> sw(R, 0), dwell(R, 0), pwm_cnt(R, 0);
  std::vector<std::vector<double>> vhist(R, std::vector<double>(256, Vnom));
  double vr = vr_v0, vr_int = vr_v0 - Vnom, vref_vr = vref_cyc ? vref_cyc[0] : Vnom;
  double E_vr = 0, E_load = 0, E_buck_in = 0, E_buck_out = 0, E_buck_sw = 0, E_ad = 0, E_adq = 0;
  long nbuck_on = 0, nsw = 0;
  for (int r = 0; r < R; ++r) { u[1 + 2 * R + r] = 0.9; un[1 + 2 * R + r] = 0.9; }
  u[0] = vr;
  for (int c = 0; c < ncyc; ++c) {
    double vref = vref_cyc ? vref_cyc[c] : Vnom;
    double vmin[16], vsum[16];
    for (int r = 0; r < R; ++r) { vmin[r] = 1e9; vsum[r] = 0; }
    for (int s = 0; s < sub; ++s) {
      long k = (long)c * sub + s;
      double vd[16];
      for (int r = 0; r < R; ++r) { vd[r] = x[vidx[r]]; vsf[r] += (vd[r] - vsf[r]) * dt / 1e-9; vhist[r][k & 255] = vsf[r]; }   // 1 ns sense filter
      double vm = 0; for (int r = 0; r < R; ++r) vm += vd[r]; vm /= R;
      { double dv = vref - vref_vr, lim = slew_vr * dt; vref_vr += std::max(-lim, std::min(lim, dv)); }   // VR DVFS slew limit
      vr_int += dt * 2 * M_PI * fc_vr * (vref_vr - vm);
      vr = vref_vr + vr_int;
      for (int r = 0; r < R; ++r) {
        const Dev& d = dev[r];
        double inj = 0.0;
        if (d.buck_mode == 1) {
          double vs = vhist[r][(k - d.dly_b) & 255];
          bool trig = (vs < vref - d.th_b) || (d.use_pred_b && pred[(long)c * R + r]);
          if (trig) itgt[r] = d.Ipk;
          else itgt[r] = std::max(0.0, itgt[r] - d.Ipk * dt / d.t_taper);
          if (dwell[r] > 0) --dwell[r];
          else {
            int want = sw[r];
            if (iL[r] < itgt[r] - d.band) want = 1;
            else if (iL[r] > itgt[r] + d.band || itgt[r] <= 0) want = 0;
            if (want != sw[r]) { sw[r] = want; dwell[r] = d.dwell; E_buck_sw += d.Esw; ++nsw; }
          }
          double vsw = sw[r] ? d.Vin : 0.0;
          if (iL[r] <= 0 && !sw[r]) iL[r] = 0.0;
          else { iL[r] += (vsw - vd[r] - iL[r] * d.Rs) / d.L * dt; if (iL[r] < 0) iL[r] = 0.0; }
          if (sw[r]) { E_buck_in += d.Vin * iL[r] * dt; ++nbuck_on; }
          inj = iL[r];
        } else if (d.buck_mode == 2) {
          // averaged current-mode buck: outer PI voltage loop (crossover fc_b on the region capacitance Cr), inner current
          // tracking limited by the inductor slopes (Vin - v)/L up and v/L down; switching loss charged at fsw
          double err = vref - vd[r];
          double kp = 2 * M_PI * d.fc_b * d.Cr, ki = kp * 2 * M_PI * d.fc_b / 4;
          if (!d.no_int) integ[r] += ki * err * dt;
          double icmd = kp * err + integ[r];                 // synchronous buck in CCM: can source and sink
          double up = (d.Vin - vd[r]) / d.L * dt, dn = vd[r] / d.L * dt;
          iL[r] += std::max(-dn, std::min(up, icmd - iL[r]));
          double duty = std::min(1.0, std::max(0.0, (vd[r] + iL[r] * d.Rs) / d.Vin));
          E_buck_in += d.Vin * iL[r] * duty * dt;
          E_buck_sw += 2 * d.Esw * d.fsw * dt * (iL[r] > 1e-4 || !d.no_int ? 1.0 : 0.0);
          inj = iL[r];
        }
        un[1 + R + r] = inj;
        if (d.buck_mode) E_buck_out += vd[r] * inj * dt;
        if (d.ad_on) {
          vlp[r] += (vd[r] - vlp[r]) * dt * 2 * M_PI * d.f_lp;
          double tgt = 0.9 - d.G * (vd[r] - vlp[r]);
          vbp[r] += (tgt - vbp[r]) * dt / d.tau_a;
          vbp[r] = std::min(1.75, std::max(0.05, vbp[r]));
          un[1 + 2 * R + r] = vbp[r];
          double ia = x[d.ad_br];
          E_ad += 1.8 * std::fabs(ia) * dt; E_adq += d.Pq * dt;
        } else un[1 + 2 * R + r] = 0.9;
      }
      for (int r = 0; r < R; ++r) {
        double Ic = load_mA[(long)c * R + r];
        double dyn = std::max(0.0, Ic - leak_mA[r]);
        un[1 + r] = 1e-3 * (leak_mA[r] + dyn * shape[s]);
        E_load += vd[r] * un[1 + r] * dt;
      }
      un[0] = vr;
      for (int i = 0; i < nx; ++i) {
        double acc = 0; const double* Pi = P + (long)i * nx;
        for (int j = 0; j < nx; ++j) acc += Pi[j] * x[j];
        const double* Qi = Q + (long)i * nu;
        for (int j = 0; j < nu; ++j) acc += Qi[j] * (u[j] + un[j]);
        xn[i] = acc;
      }
      std::swap(x, xn); std::swap(u, un);
      for (int r = 0; r < R; ++r) { double v = x[vidx[r]]; if (v < vmin[r]) vmin[r] = v; vsum[r] += v; }
      E_vr += vr * x[vidx[R]] * dt;
    }
    for (int r = 0; r < R; ++r) { vmin_out[(long)c * R + r] = (float)vmin[r]; if (vavg_out) vavg_out[(long)c * R + r] = (float)(vsum[r] / sub); }
  }
  for (int i = 0; i < nx; ++i) x0[i] = x[i];
  stats[0] = E_vr; stats[1] = E_load; stats[2] = E_buck_in; stats[3] = E_buck_out; stats[4] = E_buck_sw;
  stats[5] = E_ad; stats[6] = E_adq; stats[7] = (double)nbuck_on; stats[8] = (double)nsw;
  return 0;
}
}
