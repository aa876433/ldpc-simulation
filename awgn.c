#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include "awgn.h"
#include "mt_random.h"

#define TWO_PI (6.283185306)

double normal(void)
{
    static uint32_t generate = 0;
    static double z0, z1;
    generate = !generate;

    if (!generate)
        return z1;

    double u1, u2;
    do
    {
        u1 = gen_rand() * (1.0 / UINT32_MAX);
        u2 = gen_rand() * (1.0 / UINT32_MAX);
    } while (u1 <= DBL_MIN);

    double t1 = sqrt(-2.0 * log(u1));
    double t2 = TWO_PI * u2;
    z0 = t1 * cos(t2);
    z1 = t1 * sin(t2);
    return z0;
}

double get_awgn_noise(double sigma)
{
    return normal() * sigma;
}

double get_sigma_with_snr(double snr_db, double code_rate)
{
    return sqrt(1.0 / (2 * code_rate * pow(10, snr_db / 10)));
}
