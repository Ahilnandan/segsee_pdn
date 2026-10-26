"""Uop issue instances from a MediumBoom VCD (issue-queue slot request/grant), for dependency-aware replay.

Per IQ (int / mem / fp, 16 collapsing slots), sampled like the LACPo features (clock rise, values after all changes at
that timestamp). A uop is tracked by (iq, rob_idx) while it sits in the queue: r = first cycle its slot requests issue,
t = cycle it is granted. Each grant emits one instance (speculatively issued and re-issued uops emit two).
Output npz arrays (one entry per instance, sorted by t): iq, uopc, fu, rob, pdst, dtype, prs1, t1, prs2, t2, r, t,
plus per-cycle counts nvalid[q], nreq[q]."""
import argparse, os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_vcd2features import read_header, to_int

T = 'TOP.TestHarness.TestHarness.boom_tile.core.'
IQS = [('int', 'int_issue_unit'), ('mem', 'mem_issue_unit'), ('fp', 'fp_pipeline.fp_issue_unit')]
F = ['valid', 'request', 'grant', 'uop_fu_code', 'uop_pdst', 'uop_prs1', 'uop_prs2', 'uop_rob_idx', 'uop_uopc',
     'uop_dst_rtype', 'uop_lrs1_rtype', 'uop_lrs2_rtype']
NS = 16

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--vcd', required=True); ap.add_argument('--out', required=True); a = ap.parse_args()
    n2c = read_header(a.vcd)
    code = {}                                   # (q, k, field) -> vcd code
    for q, (_, sc) in enumerate(IQS):
        for k in range(NS):
            for f in F:
                nm = '%s%s.issue_slots_%d_%s' % (T, sc, k, f)
                if nm not in n2c: raise SystemExit('missing ' + nm)
                code[(q, k, f)] = n2c[nm]
    clk = n2c[T + 'lsu.clock']
    want = set(code.values()); cur = {c: 0 for c in want}
    slots = [[{f: code[(q, k, f)] for f in F} for k in range(NS)] for q in range(len(IQS))]
    pend = [dict() for _ in IQS]               # rob -> first request cycle
    inst = []; nvalid = []; nreq = []; cyc = 0; rose = False; started = False
    def edge():
        nv = [0] * len(IQS); nr = [0] * len(IQS)
        for q in range(len(IQS)):
            P = pend[q]
            for s in slots[q]:
                if not cur[s['valid']]: continue
                nv[q] += 1; rob = cur[s['uop_rob_idx']]
                if cur[s['request']]:
                    nr[q] += 1
                    if rob not in P: P[rob] = cyc
                if cur[s['grant']]:
                    r = P.pop(rob, cyc)
                    inst.append((q, cur[s['uop_uopc']], cur[s['uop_fu_code']], rob, cur[s['uop_pdst']], cur[s['uop_dst_rtype']],
                                 cur[s['uop_prs1']], cur[s['uop_lrs1_rtype']], cur[s['uop_prs2']], cur[s['uop_lrs2_rtype']], r, cyc))
        nvalid.append(nv); nreq.append(nr)
    with open(a.vcd) as fh:
        for line in fh:
            if not started:
                if '$enddefinitions' in line: started = True
                continue
            c0 = line[0] if line else ''
            if c0 == '#':
                if rose: edge(); cyc += 1; rose = False
                continue
            if c0 in 'bB':
                p = line[1:].split()
                if len(p) == 2 and p[1] in want: cur[p[1]] = to_int(p[0])
            elif c0 in '01xXzZ':
                c = line[1:].strip()
                if c == clk and c0 == '1': rose = True
                if c in want: cur[c] = 1 if c0 == '1' else 0
        if rose: edge(); cyc += 1
    A = np.array(inst, dtype=np.int64).reshape(-1, 12); A = A[np.argsort(A[:, 11], kind='stable')]
    names = ['iq', 'uopc', 'fu', 'rob', 'pdst', 'dtype', 'prs1', 't1', 'prs2', 't2', 'r', 't']
    np.savez_compressed(a.out, ncyc=cyc, nvalid=np.array(nvalid, np.int16), nreq=np.array(nreq, np.int16), **{n: A[:, i] for i, n in enumerate(names)})
    print('cycles', cyc, 'instances', len(A), 'per IQ', np.bincount(A[:, 0], minlength=3), '->', a.out)

if __name__ == '__main__':
    main()
