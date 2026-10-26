// pdn4: PDN time-domain kernel with role-switchable regulators.
// Linear network x_{k+1} = P x_k + Q (u_k + u_{k+1}) (trapezoidal on E x' = A x + B u), devices explicit per step.
// Inputs u: [0] off-chip buck inductor current into board node vo (A)
//           [1..R] load currents drawn at die nodes, [R+1..2R] on-die injections into die nodes (buck + LDO)
//           [2R+1..3R] active-decap bottom-plate drive voltages, [3R+1..4R] series-LDO currents drawn from rail nodes
// Devices
//   off-chip buck (averaged, L1 from Vin1): role per cycle: 0 = voltage setter (current-mode PI on remote-sensed
//     die/rail voltage, crossover fc_off, DVFS slew limit), 1 = current helper (command = integral of the on-die
//     regulators' output current -> takes over their DC, plus DVFS feed-forward C_board * dVref/dt).
//     Bumpless transfer: the command integrator is preloaded with the present inductor current on every switch.
//   on-die regulator per region, role per cycle: 1 = transient injector (hysteretic current target, switch level,
//     sensed-droop or predictor trigger, taper), 2 = voltage regulator (averaged current-mode, PI at fc_on,
//     inductor slope limits, N interleaved phases; ripple handled analytically in Python), 0 = off.
//     Predictive boost: in role 2 the reference is raised by boost_on while the issue-queue predictor flags droop.
//   LDO per region: parallel (sources current from a separate rail V_ldo through a pass device, gm loop toward
//     vref - offset (+ predictive boost), first-order bandwidth) or series (bumps feed a rail node; the LDO is the
//     only path from rail node to die).
//   active decap (Meng/Saleh): Ca between die node and an amplifier-driven bottom plate (linear in the network).
#include <cmath>
#include <vector>
#include <algorithm>
extern "C" {

struct Glob {
  double Vin1, L1, R1, fc_off, Cb, slew, f_share, ff_C;   // off-chip stage; Cb = board-side capacitance for loop gain
  int sense_rail;                                         // 1: off-chip senses rail nodes (series LDO), else die nodes
  double rail_off;                                        // off-chip regulates to vref + rail_off
};
struct Dev {
  // on-die regulator
  int has_reg; double Vin, L, Rs, Ipk, band, th_b; int dly_b, dwell; double t_taper;
  int phases; double fsw, fc_on, Cr, Esw, boost_on;
  // LDO
  int ldo; double V_ldo, R_pass, gm, tau_ldo, off_ldo, boost_ldo, Iq_ldo; int ldo_series;
  // active decap
  int ad_on; double G, tau_a, f_lp, Pq; int ad_br;
  int cp; double I_ch, R_cp, eta_cp, Vmul;          // on-die charge pump filling the LDO's reservoir node from Vdd
};

int run(int ncyc, int sub, int nx, int R, const double* P, const double* Q, const int* vidx, const int* ridx, int vo_idx,
        const float* load_mA, const double* leak_mA, const double* shape, const unsigned char* pred,
        const unsigned char* role_off, const unsigned char* role_on, const Glob* g, const Dev* dev,
        const double* vref_cyc, const double* vref_off_cyc, double Vnom, double dt, double iL1_0,
        double* x0, float* vmin_out, float* vavg_out, float* ion_out, double* stats)
{
  const int nu = 2 + 4 * R;
  std::vector<double> x(x0, x0 + nx), xn(nx), u(nu, 0.0), un(nu, 0.0);
  std::vector<double> iL(R, 0.0), itgt(R, 0.0), integ(R, 0.0), vbp(R, 0.9), vlp(R, Vnom), vsf(R, Vnom), ildo(R, 0.0);
  std::vector<int> sw(R, 0), dwell(R, 0), last_role(R, -1);
  std::vector<std::vector<double>> vhist(R, std::vector<double>(256, Vnom));
  double iL1 = iL1_0, integ_off = iL1_0, vref_vr = vref_off_cyc ? vref_off_cyc[0] : (vref_cyc ? vref_cyc[0] : Vnom), vref_prev = vref_cyc ? vref_cyc[0] : Vnom;
  int last_off = -1;
  double E_load = 0, E_off_cond = 0, E_reg_cond = 0, E_reg_sw = 0, E_reg_in = 0, E_reg_out = 0, E_inj_in = 0, E_inj_out = 0,
         E_inj_sw = 0, E_ldo = 0, E_ad = 0, E_adq = 0, E_ldoq = 0, E_cp = 0;
  for (int r = 0; r < R; ++r) { u[1 + 2 * R + r] = 0.9; un[1 + 2 * R + r] = 0.9; }
  double vrail = 0; for (int r = 0; r < R; ++r) if (dev[r].ldo && !dev[r].ldo_series) vrail = dev[r].V_ldo;
  u[0] = iL1; u[1 + 4 * R] = vrail; un[1 + 4 * R] = vrail;
  for (int c = 0; c < ncyc; ++c) {
    double vref = vref_cyc ? vref_cyc[c] : Vnom;
    double vmin[16], vsum[16], isum[16];
    for (int r = 0; r < R; ++r) { vmin[r] = 1e9; vsum[r] = 0; isum[r] = 0; }
    int ro = role_off ? role_off[c] : 0;
    for (int s = 0; s < sub; ++s) {
      long k = (long)c * sub + s;
      double vd[16], vr_[16];
      for (int r = 0; r < R; ++r) {
        vd[r] = x[vidx[r]]; vr_[r] = ridx[r] >= 0 ? x[ridx[r]] : 0.0;
        vsf[r] += (vd[r] - vsf[r]) * dt / 1e-9; vhist[r][k & 255] = vsf[r];
      }
      // ---------------- on-die devices
      double ion_tot = 0.0;
      for (int r = 0; r < R; ++r) {
        const Dev& d = dev[r];
        int role = d.has_reg ? (role_on ? role_on[(long)c * R + r] : d.has_reg) : 0;
        bool pr = pred && pred[(long)c * R + r];
        double inj = 0.0;
        if (role != last_role[r]) {                       // bumpless transfer between roles
          if (role == 2) integ[r] = iL[r];
          if (role == 1) itgt[r] = iL[r];
          last_role[r] = role;
        }
        if (role == 1) {
          double vs = vhist[r][(k - d.dly_b) & 255];
          bool trig = (vs < vref - d.th_b) || pr;
          if (trig) itgt[r] = d.Ipk; else itgt[r] = std::max(0.0, itgt[r] - d.Ipk * dt / d.t_taper);
          if (dwell[r] > 0) --dwell[r];
          else {
            int want = sw[r];
            if (iL[r] < itgt[r] - d.band) want = 1;
            else if (iL[r] > itgt[r] + d.band || itgt[r] <= 0) want = 0;
            if (want != sw[r]) { sw[r] = want; dwell[r] = d.dwell; E_inj_sw += d.Esw; }
          }
          double vsw = sw[r] ? d.Vin : 0.0;
          if (iL[r] <= 0 && !sw[r]) iL[r] = 0.0;
          else { iL[r] += (vsw - vd[r] - iL[r] * d.Rs) / d.L * dt; if (iL[r] < 0) iL[r] = 0.0; }
          if (sw[r]) E_inj_in += d.Vin * iL[r] * dt;
          E_inj_out += vd[r] * iL[r] * dt;
          inj += iL[r];
        } else if (role == 2) {
          double vt = vref + (pr ? d.boost_on : 0.0);
          double err = vt - vd[r];
          double kp = 2 * M_PI * d.fc_on * d.Cr, ki = kp * 2 * M_PI * d.fc_on / 4;
          integ[r] += ki * err * dt;
          double icmd = kp * err + integ[r];
          double Leff = d.L / d.phases;                    // N phases in parallel: total current slope N x
          double up = (d.Vin - vd[r]) / Leff * dt, dn = vd[r] / Leff * dt;
          iL[r] += std::max(-dn, std::min(up, icmd - iL[r]));
          double duty = std::min(1.0, std::max(0.0, (vd[r] + iL[r] * d.Rs / d.phases) / d.Vin));
          E_reg_in += d.Vin * iL[r] * duty * dt;
          E_reg_out += vd[r] * iL[r] * dt;
          E_reg_cond += iL[r] * iL[r] * d.Rs / d.phases * dt;
          E_reg_sw += 2 * d.Esw * d.fsw * d.phases * dt;
          inj += iL[r];
        } else { iL[r] = 0.0; }
        ion_tot += inj;
        // LDO
        if (d.ldo) {
          double src = vr_[r];                                 // rail node: series (main path) or separate LDO rail
          double tgt = vref - d.off_ldo + (pr ? d.boost_ldo : 0.0);
          double icmd = std::max(0.0, d.gm * (tgt - vd[r]));
          double imax = std::max(0.0, (src - vd[r]) / d.R_pass);
          icmd = std::min(icmd, imax);
          ildo[r] += (icmd - ildo[r]) * std::min(1.0, dt / d.tau_ldo);
          ildo[r] = std::min(ildo[r], imax);
          E_ldo += (src - vd[r]) * ildo[r] * dt; E_ldoq += d.Iq_ldo * src * dt;
          inj += ildo[r];
          double ich = 0.0;
          if (d.cp) {                                           // N:1 pump from Vdd into the reservoir, current-limited
            ich = std::min(d.I_ch, std::max(0.0, (d.Vmul * vd[r] - src) / d.R_cp));
            double iin = d.Vmul * ich / d.eta_cp;               // drawn from Vdd
            inj -= iin; E_cp += (vd[r] * iin - src * ich) * dt;
          }
          un[1 + 3 * R + r] = ildo[r] - ich;                    // net drawn from the rail/reservoir node
        } else un[1 + 3 * R + r] = 0.0;
        un[1 + R + r] = inj;
        isum[r] += inj;
        if (d.ad_on) {
          vlp[r] += (vd[r] - vlp[r]) * dt * 2 * M_PI * d.f_lp;
          double tgt = 0.9 - d.G * (vd[r] - vlp[r]);
          vbp[r] += (tgt - vbp[r]) * dt / d.tau_a;
          vbp[r] = std::min(1.75, std::max(0.05, vbp[r]));
          un[1 + 2 * R + r] = vbp[r];
          E_ad += 1.8 * std::fabs(x[d.ad_br]) * dt; E_adq += d.Pq * dt;
        } else un[1 + 2 * R + r] = 0.9;
      }
      // ---------------- off-chip buck
      {
        double vtar = vref_off_cyc ? vref_off_cyc[c] : vref;
        double dvr = vtar - vref_vr, lim = g->slew * dt; vref_vr += std::max(-lim, std::min(lim, dvr));
        if (ro != last_off) { integ_off = iL1; last_off = ro; }
        double icmd;
        if (ro == 0) {
          double vm = 0; for (int r = 0; r < R; ++r) vm += (g->sense_rail ? vr_[r] : vd[r]); vm /= R;
          double err = vref_vr + g->rail_off - vm;
          double kp = 2 * M_PI * g->fc_off * g->Cb, ki = kp * 2 * M_PI * g->fc_off / 4;
          integ_off += ki * err * dt; icmd = kp * err + integ_off;
        } else {
          integ_off += 2 * M_PI * g->f_share * ion_tot * dt;          // absorb the on-die regulators' DC
          icmd = integ_off + g->ff_C * (vref - vref_prev) / (sub * dt);
        }
        double vo = x[vo_idx];
        double up = (g->Vin1 - vo) / g->L1 * dt, dn = vo / g->L1 * dt;
        iL1 += std::max(-dn, std::min(up, icmd - iL1));
        E_off_cond += iL1 * iL1 * g->R1 * dt;
        un[0] = iL1;
      }
      for (int r = 0; r < R; ++r) {
        double Ic = load_mA[(long)c * R + r];
        double dyn = std::max(0.0, Ic - leak_mA[r]);
        un[1 + r] = 1e-3 * (leak_mA[r] + dyn * shape[s]);
        E_load += vd[r] * un[1 + r] * dt;
      }
      for (int i = 0; i < nx; ++i) {
        double acc = 0; const double* Pi = P + (long)i * nx;
        for (int j = 0; j < nx; ++j) acc += Pi[j] * x[j];
        const double* Qi = Q + (long)i * nu;
        for (int j = 0; j < nu; ++j) acc += Qi[j] * (u[j] + un[j]);
        xn[i] = acc;
      }
      std::swap(x, xn); std::swap(u, un);
      for (int r = 0; r < R; ++r) { double v = x[vidx[r]]; if (v < vmin[r]) vmin[r] = v; vsum[r] += v; }
    }
    vref_prev = vref;
    for (int r = 0; r < R; ++r) {
      vmin_out[(long)c * R + r] = (float)vmin[r];
      if (vavg_out) vavg_out[(long)c * R + r] = (float)(vsum[r] / sub);
      if (ion_out) ion_out[(long)c * R + r] = (float)(isum[r] / sub);
    }
  }
  for (int i = 0; i < nx; ++i) x0[i] = x[i];
  stats[0] = E_load; stats[1] = E_off_cond; stats[2] = E_reg_cond; stats[3] = E_reg_sw; stats[4] = E_reg_in; stats[5] = E_reg_out;
  stats[6] = E_inj_in; stats[7] = E_inj_out; stats[8] = E_inj_sw; stats[9] = E_ldo; stats[10] = E_ad; stats[11] = E_adq; stats[12] = E_ldoq; stats[13] = E_cp;
  return 0;
}
}
