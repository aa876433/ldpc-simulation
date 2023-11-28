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
    uint32_t max_iter = 50;
    uint32_t early_term = 0;
    DECODE_METHOD dec_method = LAYERED_SPA_ALGORITHM;

    for (uint32_t i = 175; i <= 200; i += 25)
    {
        double snr_db = (i / 100.0);
        start = clock();
        ldpc_set_config(p_ldpc_info, snr_db, max_iter, early_term, dec_method);
        ERROR_RATE_T error_rate = ldpc_simulation(p_ldpc_info);
        end = clock();
        printf("%.2f: %f, %f, %f s\n", snr_db, error_rate.bit_error_rate, error_rate.block_error_rate, ((double)(end - start) / CLOCKS_PER_SEC));
    }

    ldpc_release(p_ldpc_info);
    release_parity_info(parity_info);
    return 0;
}
