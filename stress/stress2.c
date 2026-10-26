// Miss-driven di/dt stressmark (v2). Pointer-chase rings at 128 B stride in free DRAM.
#include "util.h"
extern void chase1(long, void *, double *), chase3(long, void *, double *), idle2(long);
extern void chase2x(long, void *, double *, void *);
#define LINE 16                                 /* 16 longs = 128 B stride */
#define NBIG (2 * 1024 * 1024 / 128)            /* 2 MB ring: misses L1 (and L2 if present) */
#define NSMALL (64 * 1024 / 128)                /* 64 KB ring: misses L1 */
#define big ((long *)0x80400000UL)              /* rings in free DRAM: not loaded, not zeroed by the loader */
#define small ((long *)0x80300000UL)
static double buf[8] = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
static void ring(long *a, long n) {             // linear ring at 128 B stride: every load touches a new line
    for (long i = 0; i < n; i++) a[i * LINE] = (long)&a[((i + 1) % n) * LINE];
}
int main(void) {
    ring(small, NSMALL); ring(big, NBIG);
    idle2(150);
    chase1(600, small, buf); idle2(150);
    chase1(600, big, buf);   idle2(150);
    chase3(600, &big[(NBIG / 4) * LINE], buf);   idle2(150);
    chase2x(400, &big[(NBIG / 2) * LINE], buf, &big[(3 * NBIG / 4) * LINE]); idle2(150);
    return 0;
}
