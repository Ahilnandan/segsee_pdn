"""Extract per-cycle hierarchy power from PTPX time-based output (value-change format: 'Time: t' then 'idx value')."""
import sys, numpy as np
src, out = sys.argv[1], sys.argv[2]
want = {}
with open(src) as fh:
    for line in fh:
        if not line.startswith('.index'): break
        p = line.split(); name = p[1][3:-1]
        if name.count('/') <= 2: want[int(p[2])] = name
cur = {i: 0.0 for i in want}; rows = []; t = None
with open(src) as fh:
    for line in fh:
        if line[0] == '.': continue
        if line.startswith('Time:'):
            if t is not None: rows.append([cur[i] for i in want])
            t = int(line.split()[1]); continue
        p = line.split()
        if len(p) == 2:
            i = int(p[0])
            if i in cur: cur[i] = float(p[1])
    rows.append([cur[i] for i in want])
np.savez_compressed(out, P=np.array(rows, dtype=np.float32), names=np.array([want[i] for i in want]))
print(len(rows), 'cycles', len(want), 'hierarchies')
