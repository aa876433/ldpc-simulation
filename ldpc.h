#ifndef LDPC_H
#define LDPC_H

#include <stdint.h>

struct LDPC_INFO_T;

typedef enum DECODE_METHOD_
{
    SPA_ALGORITHM,
    LAYERED_SPA_ALGORITHM,
    MS_ALGORITHM,
    LAYERED_MS_ALGORITHM,
    MAX_ALGORITHM,
} DECODE_METHOD;

typedef struct ERROR_RATE_T
{
    double bit_error_rate;
    double block_error_rate;
} ERROR_RATE_T;

struct LDPC_INFO_T *ldpc_init(void *info);

void ldpc_set_config(struct LDPC_INFO_T *p_ldpc_info, double snr_db, uint32_t max_iter, uint32_t early_term, DECODE_METHOD method);

ERROR_RATE_T ldpc_simulation(struct LDPC_INFO_T *);

void ldpc_release(struct LDPC_INFO_T *p_ldpc_info);

#endif //LDPC_H
