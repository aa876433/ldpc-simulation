#ifndef LDPC_H
#define LDPC_H

#include <stdint.h>

typedef enum DECODE_METHOD_
{
    SPA_ALGORITHM,
    MS_ALGORITHM,
    MAX_ALGORITHM,
} DECODE_METHOD;

typedef struct ERROR_RATE_T
{
    double bit_error_rate;
    double block_error_rate;
} ERROR_RATE_T;

void ldpc_init(void *parity_info);

ERROR_RATE_T ldpc_simulation(double snr_db, uint32_t iter, uint32_t early_termination, DECODE_METHOD method);

void ldpc_release(void);

#endif //LDPC_H
