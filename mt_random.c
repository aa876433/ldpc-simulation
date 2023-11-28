#include <stdint.h>
#include "mt_random.h"

#define MT_N 624
#define MT_M 397
#define MATRIX_A 0x9908b0dfUL   
#define UPPER_MASK 0x80000000UL 
#define LOWER_MASK 0x7fffffffUL 

static uint32_t mt[MT_N];
static int mti = MT_N + 1;

void init_gen_rand(uint32_t s)
{
    mt[0] = s & 0xffffffffUL;
    for (mti = 1; mti < MT_N; mti++)
    {
        mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
        mt[mti] &= 0xffffffffUL;
    }
}

uint32_t gen_rand(void)
{
    uint32_t y;
    static uint32_t mag01[2] = {0x0UL, MATRIX_A};
    
    if (mti >= MT_N)
    {
        int kk;

        if (mti == MT_N + 1)
            init_gen_rand(5489UL);

        for (kk = 0; kk < MT_N - MT_M; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk < MT_N - 1; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[MT_N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y;
}
