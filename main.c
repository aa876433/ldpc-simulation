#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include "mt_random.h"
#include "awgn.h"
#include "read_parity.h"
#include "ldpc.h"

int main()
{
    srand(time(NULL));
    init_gen_rand(rand());
    PARITY_INFO_T *parity_info = get_parity_matrix_info("h.txt");
    ldpc_init(parity_info);
    clock_t start, end;
    uint32_t early_termination = 1;

    for (uint32_t i = 42; i <= 46; i += 2)
    {
        double snr_db = (i / 10.0);
        start = clock();
        ERROR_RATE_T error_rate = ldpc_simulation(snr_db, 60, early_termination, MS_ALGORITHM);
        end = clock();
        printf("%.2f: %f, %f, %f s\n", snr_db, error_rate.bit_error_rate, error_rate.block_error_rate, ((double)(end - start) / CLOCKS_PER_SEC));
    }

    ldpc_release();
    release_parity_info(parity_info);
    return 0;
}
