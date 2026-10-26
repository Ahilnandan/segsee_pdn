// pdn5: tiled (8x8) PDN, sparse trapezoidal solver (Eigen SparseLU), steppable with retained state for closed-loop
// co-simulation with the issue replay. Linear network: Ml x_{k+1} = Mr x_k + Bq (u_k + u_{k+1}).
// Inputs u: [0] off-chip buck inductor current into board node; per tile t: [1+t] load drawn at die node,
//           [1+T+t] on-die injection into die node (LDO bank / buck), [1+2T+t] current drawn from tile reservoir node.
// Devices: off-chip buck (averaged, voltage-setter role 0 or current-helper role 1); per tile: charge-pump bank
// (2:1 pump from Vdd into reservoir, current-limited) + regulated pass device (LDO); on-die buck (role 1 injector,
// role 2 averaged regulator). Sensors: per-tile voltage with delay/quantisation, exported every cycle.
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <vector>
#include <cmath>
#include <algorithm>
typedef Eigen::SparseMatrix<double> SpM;
struct Tile {
  int bank; double R_pass, gm, tau, off, I_ch, eta, Vmul; int ridx;      // LDO + charge-pump bank
  int buck; double Vin, L, Rs, Ipk, band, th, t_taper, fc, Cr; int phases, dwell;  // on-die buck
};
struct Sim {
  int nx, T, sub; double dt, Vnom;
  SpM Mr, Bq; Eigen::SparseLU<SpM> lu;
  Eigen::VectorXd x, u, un, rhs;
  std::vector<int> didx; int vo_idx;
  std::vector<Tile> tl;
  // off-chip buck
  double Vin1, L1, R1, fc_off, Cb, slew, f_share, iL1, integ_off, vref_vr; int last_off;
  // device state
  std::vector<double> ildo, iL, itgt, integ, vs; std::vector<int> sw, dwl, last_role;
  // energy
  double E_load, E_off, E_ldo, E_cp, E_buck, E_q; long k;
};
extern "C" {
void* sim_new(int nx, int T, int sub, double dt, double Vnom,
              int nMl, const int* Mli, const int* Mlj, const double* Mlv,
              int nMr, const int* Mri, const int* Mrj, const double* Mrv,
              int nB, const int* Bi, const int* Bj, const double* Bv, int nu,
              const int* didx, int vo_idx, const double* x0, const Tile* tiles,
              double Vin1, double L1, double R1, double fc_off, double Cb, double slew, double f_share, double iL1_0)
{
  Sim* s = new Sim(); s->nx = nx; s->T = T; s->sub = sub; s->dt = dt; s->Vnom = Vnom;
  std::vector<Eigen::Triplet<double>> t;
  SpM Ml(nx, nx); for (int i = 0; i < nMl; ++i) t.emplace_back(Mli[i], Mlj[i], Mlv[i]); Ml.setFromTriplets(t.begin(), t.end()); t.clear();
  s->Mr.resize(nx, nx); for (int i = 0; i < nMr; ++i) t.emplace_back(Mri[i], Mrj[i], Mrv[i]); s->Mr.setFromTriplets(t.begin(), t.end()); t.clear();
  s->Bq.resize(nx, nu); for (int i = 0; i < nB; ++i) t.emplace_back(Bi[i], Bj[i], Bv[i]); s->Bq.setFromTriplets(t.begin(), t.end());
  Ml.makeCompressed(); s->lu.analyzePattern(Ml); s->lu.factorize(Ml);
  s->x = Eigen::Map<const Eigen::VectorXd>(x0, nx); s->u = Eigen::VectorXd::Zero(nu); s->un = s->u; s->rhs.resize(nx);
  s->didx.assign(didx, didx + T); s->vo_idx = vo_idx; s->tl.assign(tiles, tiles + T);
  s->Vin1 = Vin1; s->L1 = L1; s->R1 = R1; s->fc_off = fc_off; s->Cb = Cb; s->slew = slew; s->f_share = f_share;
  s->iL1 = iL1_0; s->integ_off = iL1_0; s->vref_vr = Vnom; s->last_off = -1;
  s->ildo.assign(T, 0); s->iL.assign(T, 0); s->itgt.assign(T, 0); s->integ.assign(T, 0); s->vs.assign(T, Vnom);
  s->sw.assign(T, 0); s->dwl.assign(T, 0); s->last_role.assign(T, -1);
  s->E_load = s->E_off = s->E_ldo = s->E_cp = s->E_buck = s->E_q = 0; s->k = 0;
  s->u[0] = iL1_0; s->un[0] = iL1_0;
  return s;
}
void sim_free(void* p) { delete (Sim*)p; }
int sim_info(void* p) { return ((Sim*)p)->lu.info(); }
// Advance ncyc clock cycles. load_mA: ncyc*T; leak_mA: T; shape: sub; vref: ncyc (or null);
// role_off: ncyc (or null); role_on: ncyc*T (or null: tile default); arm: ncyc*T predictor flags (or null).
// Outputs per cycle per tile: vmin, vavg (float), and per tile reservoir voltage at the end of each cycle (float, or null).
int sim_run(void* p, int ncyc, const float* load_mA, const double* leak_mA, const double* shape, const double* vref_cyc,
            const unsigned char* role_off, const unsigned char* role_on, const unsigned char* arm,
            float* vmin_out, float* vavg_out, float* vres_out, double* stats)
{
  Sim* s = (Sim*)p; const int T = s->T; const double dt = s->dt;
  for (int c = 0; c < ncyc; ++c) {
    double vref = vref_cyc ? vref_cyc[c] : s->Vnom; int ro = role_off ? role_off[c] : 0;
    std::vector<double> vmin(T, 1e9), vsum(T, 0.0);
    for (int q = 0; q < s->sub; ++q, ++s->k) {
      double vm = 0, ion = 0;
      for (int t = 0; t < T; ++t) { double v = s->x[s->didx[t]]; s->vs[t] += (v - s->vs[t]) * dt / 1e-9; vm += v; }
      vm /= T;
      for (int t = 0; t < T; ++t) {
        const Tile& d = s->tl[t]; double vd = s->x[s->didx[t]], inj = 0, draw = 0;
        bool pr = arm && arm[(long)c * T + t];
        if (d.bank) {
          double vr = s->x[d.ridx];
          double tgt = vref - d.off + (pr ? d.off * 0.5 : 0.0);            // predictor pre-arm: halve the threshold
          double icmd = std::max(0.0, d.gm * (tgt - vd)), imax = std::max(0.0, (vr - vd) / d.R_pass);
          icmd = std::min(icmd, imax);
          s->ildo[t] += (icmd - s->ildo[t]) * std::min(1.0, dt / d.tau); s->ildo[t] = std::min(s->ildo[t], imax);
          double ich = std::min(d.I_ch, std::max(0.0, (d.Vmul * vd - vr) / 5.0));
          if (pr) ich = 0.0;                                               // no recharge while a burst is predicted
          double iin = d.Vmul * ich / d.eta;
          inj += s->ildo[t] - iin; draw += s->ildo[t] - ich;
          s->E_ldo += (vr - vd) * s->ildo[t] * dt; s->E_cp += (vd * iin - vr * ich) * dt; s->E_q += 0.5e-3 * (1e-9 / d.tau) * vd * dt;
        }
        int role = d.buck ? (role_on ? role_on[(long)c * T + t] : d.buck) : 0;
        if (role != s->last_role[t]) { if (role == 2) s->integ[t] = s->iL[t]; if (role == 1) s->itgt[t] = s->iL[t]; s->last_role[t] = role; }
        if (role == 1) {
          bool trig = (s->vs[t] < vref - d.th) || pr;
          if (trig) s->itgt[t] = d.Ipk; else s->itgt[t] = std::max(0.0, s->itgt[t] - d.Ipk * dt / d.t_taper);
          if (s->dwl[t] > 0) --s->dwl[t];
          else { int want = s->sw[t]; if (s->iL[t] < s->itgt[t] - d.band) want = 1; else if (s->iL[t] > s->itgt[t] + d.band || s->itgt[t] <= 0) want = 0;
                 if (want != s->sw[t]) { s->sw[t] = want; s->dwl[t] = d.dwell; s->E_buck += 5e-12 * (d.Vin / 1.8) * (d.Vin / 1.8) / std::max(d.Rs * 0.5, 0.1); } }
          double vsw = s->sw[t] ? d.Vin : 0.0;
          if (s->iL[t] <= 0 && !s->sw[t]) s->iL[t] = 0; else { s->iL[t] += (vsw - vd - s->iL[t] * d.Rs) / d.L * dt; if (s->iL[t] < 0) s->iL[t] = 0; }
          s->E_buck += s->iL[t] * s->iL[t] * d.Rs * dt; inj += s->iL[t]; ion += s->iL[t];
        } else if (role == 2) {
          double err = vref - vd, kp = 2 * M_PI * d.fc * d.Cr, ki = kp * 2 * M_PI * d.fc / 4;
          s->integ[t] += ki * err * dt; double icmd = kp * err + s->integ[t];
          double Leff = d.L / d.phases, up = (d.Vin - vd) / Leff * dt, dn = vd / Leff * dt;
          s->iL[t] += std::max(-dn, std::min(up, icmd - s->iL[t]));
          s->iL[t] = std::max(-d.Ipk, std::min(d.Ipk, s->iL[t]));                    // switch / inductor current rating
          s->E_buck += (s->iL[t] * s->iL[t] * d.Rs / d.phases + 2 * 5e-12 * 100e6 * d.phases) * dt; inj += s->iL[t]; ion += s->iL[t];
        } else s->iL[t] = 0;
        s->un[1 + T + t] = inj; s->un[1 + 2 * T + t] = draw;
      }
      // off-chip buck
      { double dv = vref - s->vref_vr, lim = s->slew * dt; s->vref_vr += std::max(-lim, std::min(lim, dv));
        if (ro != s->last_off) { s->integ_off = s->iL1; s->last_off = ro; }
        double icmd;
        if (ro == 0) { double err = s->vref_vr - vm, kp = 2 * M_PI * s->fc_off * s->Cb, ki = kp * 2 * M_PI * s->fc_off / 4;
                       s->integ_off += ki * err * dt; icmd = kp * err + s->integ_off; }
        else { s->integ_off += 2 * M_PI * s->f_share * ion * dt; icmd = s->integ_off; }
        double vo = s->x[s->vo_idx], up = (s->Vin1 - vo) / s->L1 * dt, dn = vo / s->L1 * dt;
        s->iL1 += std::max(-dn, std::min(up, icmd - s->iL1)); s->E_off += s->iL1 * s->iL1 * s->R1 * dt; s->un[0] = s->iL1; }
      for (int t = 0; t < T; ++t) {
        double Ic = load_mA[(long)c * T + t], dyn = std::max(0.0, Ic - leak_mA[t]);
        s->un[1 + t] = 1e-3 * (leak_mA[t] + dyn * shape[q]); s->E_load += s->x[s->didx[t]] * s->un[1 + t] * dt;
      }
      s->rhs = s->Mr * s->x + s->Bq * (s->u + s->un);
      s->x = s->lu.solve(s->rhs); s->u = s->un;
      for (int t = 0; t < T; ++t) { double v = s->x[s->didx[t]]; if (v < vmin[t]) vmin[t] = v; vsum[t] += v; }
    }
    for (int t = 0; t < T; ++t) {
      vmin_out[(long)c * T + t] = (float)vmin[t]; if (vavg_out) vavg_out[(long)c * T + t] = (float)(vsum[t] / s->sub);
      if (vres_out) vres_out[(long)c * T + t] = s->tl[t].bank ? (float)s->x[s->tl[t].ridx] : 0.f;
    }
  }
  stats[0] = s->E_load; stats[1] = s->E_off; stats[2] = s->E_ldo; stats[3] = s->E_cp; stats[4] = s->E_buck; stats[5] = s->E_q;
  return 0;
}
}
