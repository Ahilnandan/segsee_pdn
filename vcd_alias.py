"""Rename Chipyard-1.14 scope 'element_reset_domain_boom_tile' -> 'boom_tile' in a VCD header (values untouched),
so Unnath's LACPo remap (tools/vcd2features.py, written for TOP...boom_tile) finds the tile. Streams the file."""
import sys
src, dst = sys.argv[1], sys.argv[2]
with open(src, 'rb') as f, open(dst, 'wb') as g:
    for line in f:
        if line.startswith(b'$enddefinitions'):
            g.write(line); break
        g.write(line.replace(b' element_reset_domain_boom_tile ', b' boom_tile '))
    while True:
        b = f.read(1 << 24)
        if not b: break
        g.write(b)
