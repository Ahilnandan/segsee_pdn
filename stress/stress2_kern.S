# Miss-driven di/dt stressmark: pointer chase (each load misses) -> heavy burst that depends on the returned value.
# The next chase load depends only on the pointer, so it issues right away; the core sits nearly idle until the data
# returns, then the burst fires: idle -> full burst edges once per miss latency.  a0 = iterations, a1 = start node,
# a2 = 64-byte cache-resident buffer for the burst loads.
.macro BURSTP
    mul   a3, t0, t2
    add   a4, t0, t3
    add   a5, t0, t4
    ld    a6, 0(a2)
    ld    a7, 8(a2)
    mul   t5, t0, t3
    sub   t6, t0, t4
    ld    t1, 16(a2)
    fmadd.d f0, f1, f2, f3
    fmadd.d f4, f5, f6, f7
    fmadd.d f8, f1, f6, f3
    xor   s2, t0, t3
.endm
.macro CHASE name, NB
.globl \name
\name:
    addi  sp, sp, -16
    sd    s2, 0(sp)
    li    t2, 5
    li    t3, 7
    li    t4, 11
    fcvt.d.l f1, t2
    fcvt.d.l f2, t3
    fcvt.d.l f3, t4
1:
    ld    a1, 0(a1)            # next node (cache miss)
    mv    t0, a1               # burst waits for the returned pointer
    .rept \NB
    BURSTP
    .endr
    addi  a0, a0, -1
    bnez  a0, 1b
    ld    s2, 0(sp)
    addi  sp, sp, 16
    ret
.endm
CHASE chase1, 1
CHASE chase3, 3
# two independent chains: bursts from either return
.globl chase2x
chase2x:
    addi  sp, sp, -16
    sd    s2, 0(sp)
    sd    s3, 8(sp)
    mv    s3, a3               # second chain pointer (BURSTP overwrites a3)
    li    t2, 5
    li    t3, 7
    li    t4, 11
2:
    ld    a1, 0(a1)
    ld    s3, 0(s3)
    mv    t0, a1
    BURSTP
    mv    t0, s3
    BURSTP
    addi  a0, a0, -1
    bnez  a0, 2b
    ld    s3, 8(sp)
    ld    s2, 0(sp)
    addi  sp, sp, 16
    ret
.globl idle2
idle2:
    li    t0, 1000003
    li    t1, 3
    li    t2, 1000003
3:
    div   t0, t0, t1
    add   t0, t0, t2
    addi  a0, a0, -1
    bnez  a0, 3b
    ret
