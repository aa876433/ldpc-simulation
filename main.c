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
    PARITY_INFO_T *parity_info = get_parity_matrix_info("h_1944.txt");
    struct LDPC_INFO_T *p_ldpc_info = ldpc_init(parity_info);
    clock_t start, end;
    uint32_t times = 600;
    uint32_t max_iter = 25;
    uint32_t early_term = 1;
    SCHEDULE_METHOD sch = FLOOD_ALGORITHM;
    BP_METHOD bp = SC_MS_ALGORITHM;
    double norm_factor = 1;

    ldpc_set_cond_config(p_ldpc_info, max_iter, early_term);
    ldpc_set_dec_config(p_ldpc_info, sch, bp, norm_factor);

    for (uint32_t i = 100; i <= 200; i += 25)
    {
        double snr_db = (i / 100.0);
        start = clock();
        ERROR_RATE_T error_rate = ldpc_simulation(p_ldpc_info, snr_db, times);
        end = clock();
        printf("%.2f: %f, %f, %f s\n", snr_db, error_rate.bit_error_rate, error_rate.block_error_rate,
               ((double) (end - start) / CLOCKS_PER_SEC));
    }

    ldpc_release(p_ldpc_info);
    release_parity_info(parity_info);
    return 0;
}
