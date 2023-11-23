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
    PARITY_INFO_T *parity_info = get_parity_matrix_info("parity_matrix.txt");
    ldpc_init(parity_info);
    clock_t start, end;

    for (uint32_t i = 14; i <= 20; i += 2)
    {
        double snr_db = (i / 10.0);
        start = clock();
        double error_rate = ldpc_simulation(snr_db, 20, SPA_ALGORITHM);
        end = clock();
        printf("%.2f: %f, %f s\n", snr_db, error_rate, ((double)(end - start) / CLOCKS_PER_SEC));
    }

    ldpc_release();
    release_parity_info(parity_info);
    return 0;
}
