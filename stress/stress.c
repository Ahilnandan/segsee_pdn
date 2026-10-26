// di/dt stressmark: phases separated by long idle gaps (dependent divides) so they can be found in the current trace.
#include "util.h"
extern void stress_d0(long, double *), stress_d1(long, double *), stress_d2(long, double *), stress_d3(long, double *);
extern void stress_d4(long, double *), stress_d6(long, double *), stress_slow(long, double *), stress_idle(long);
static double buf[8] = {1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5};
int main(void) {
    stress_idle(150);
    stress_d0(2500, buf); stress_idle(150);
    stress_d1(2200, buf); stress_idle(150);
    stress_d2(2000, buf); stress_idle(150);
    stress_d3(1800, buf); stress_idle(150);
    stress_d4(1600, buf); stress_idle(150);
    stress_d6(1400, buf); stress_idle(150);
    stress_slow(150, buf); stress_idle(150);
    for (int k = 0; k < 4; k++) { stress_slow(25, buf); stress_idle(60); }   // big steps: ~2k heavy / ~2k idle
    return 0;
}
