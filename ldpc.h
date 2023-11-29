#ifndef LDPC_H
#define LDPC_H

#include <stdint.h>

struct LDPC_INFO_T;

typedef enum SCHEDULE_METHOD
{
    FLOOD_ALGORITHM,
    LAYER_ALGORITHM,
    MAX_SCHEDULE_METHOD,
} SCHEDULE_METHOD;

typedef enum BP_METHOD
{
    SPA_ALGORITHM,
    MS_ALGORITHM,
    SC_MS_ALGORITHM,
    MAX_BP_METHOD,
} BP_METHOD;

typedef struct ERROR_RATE_T
{
    double bit_error_rate;
    double block_error_rate;
} ERROR_RATE_T;

struct LDPC_INFO_T *ldpc_init(void *info);

void ldpc_set_cond_config(struct LDPC_INFO_T *p_ldpc_info, uint32_t iter, uint32_t early_term);
void ldpc_set_dec_config(struct LDPC_INFO_T *p_ldpc_info, SCHEDULE_METHOD sch, BP_METHOD bp, double ms_norm_factor);
ERROR_RATE_T ldpc_simulation(struct LDPC_INFO_T *p_ldpc_info, double snr_db, uint32_t times);

void ldpc_release(struct LDPC_INFO_T *p_ldpc_info);

#endif //LDPC_H
