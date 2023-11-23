#ifndef LDPC_H
#define LDPC_H

#include <stdint.h>

typedef enum DECODE_METHOD_ {
    SPA_ALGORITHM,
    MS_ALGORITHM,
} DECODE_METHOD;

void ldpc_init(void *parity_info);

double ldpc_simulation(double snr_db, uint32_t iter, DECODE_METHOD method);

void ldpc_release(void);

#endif //LDPC_H
