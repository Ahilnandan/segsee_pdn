# di/dt stressmark kernels for MediumBoom (2 ALU incl. 1 MUL, 1 MEM, 1 FPU issue ports).
# Each resonant phase: a serial chain of D dependent adds (low current) -> a burst of independent heavy ops that all
# depend on the chain result (high current) -> the next chain depends on the burst. Period ~ D + burst latency.
# a0 = iterations, a1 = pointer to a 64-byte cache-resident buffer.
.macro BURST
    mul   a2, t0, t2
    add   a3, t0, t3
    add   a4, t0, t4
    ld    a5, 0(a1)
    ld    a6, 8(a1)
    mul   a7, t0, t5
    sub   t6, t0, t4
    ld    s2, 16(a1)
    fmadd.d f0, f1, f2, f3
    fmadd.d f4, f5, f6, f7
    fmadd.d f8, f1, f6, f3
    xor   t1, t0, t3
.endm
.macro RES name, D
.globl \name
\name:
    addi  sp, sp, -16
    sd    s2, 0(sp)
    li    t0, 1
    li    t1, 3
    li    t2, 5
    li    t3, 7
    li    t4, 11
    li    t5, 13
1:
    .rept \D
    add   t0, t0, t1
    .endr
    BURST
    xor   t0, a2, a3            # next chain waits on the burst (mul result)
    add   t0, t0, a5
    addi  a0, a0, -1
    bnez  a0, 1b
    ld    s2, 0(sp)
    addi  sp, sp, 16
    ret
.endm
RES stress_d0, 0
RES stress_d1, 1
RES stress_d2, 2
RES stress_d3, 3
RES stress_d4, 4
RES stress_d6, 6
# second-droop phase: ~45 cycles of independent heavy work, then ~45 cycles of a dependent chain
.globl stress_slow
stress_slow:
    addi  sp, sp, -16
    sd    s2, 0(sp)
    li    t0, 1
    li    t1, 3
    li    t2, 5
    li    t3, 7
    li    t4, 11
    li    t5, 13
2:
    .rept 30
    BURST
    .endr
    .rept 45
    add   t0, t0, t1
    .endr
    addi  a0, a0, -1
    bnez  a0, 2b
    ld    s2, 0(sp)
    addi  sp, sp, 16
    ret
# idle gap / long low phase: dependent divides (iterative divider, ~1 uop in flight)
.globl stress_idle
stress_idle:
    li    t0, 1000003
    li    t1, 3
    li    t2, 1000003
3:
    div   t0, t0, t1
    add   t0, t0, t2
    addi  a0, a0, -1
    bnez  a0, 3b
    ret
