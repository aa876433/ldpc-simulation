#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "ldpc.h"
#include "awgn.h"
#include "mt_random.h"
#include "read_parity.h"

typedef struct SPARSE_INFO_T_
{
    uint32_t len;
    uint32_t pos[0];
} SPARSE_INFO_T;

typedef struct LDPC_INFO_T_
{
    SPARSE_INFO_T **sparse_gen;
    SPARSE_INFO_T **sparse_row;
    uint8_t **parity_matrix;
    double **cv_matrix;
    uint32_t data_len;
    uint32_t parity_len;
    uint32_t codeword_len;
    double code_rate;
} LDPC_INFO_T;

static LDPC_INFO_T g_ldpc;

double sgn(double x)
{
    if (x < 0) return -1;
    if (x > 0) return 1;
    return 0;
}

double stable_atanh(double t)
{
    double epsilon = 1e-12;
    if (t > 1 - epsilon)
    {
        return 14;
    }
    else if (t < -1 + epsilon)
    {
        return -14;
    }
    return atanh(t);
}

void xor_to_row(const uint8_t *row_from, uint8_t *row_to, uint32_t len)
{
    for (uint32_t i = 0; i < len; i++)
    {
        row_to[i] ^= row_from[i];
    }
}

void swap_row(uint8_t *row_a, uint8_t *row_b, uint32_t len)
{
    for (uint32_t i = 0; i < len; i++)
    {
        uint8_t t = row_a[i];
        row_a[i] = row_b[i];
        row_b[i] = t;
    }
}

void simplify_parity_matrix(void)
{
    uint32_t k = g_ldpc.parity_len;
    uint32_t n = g_ldpc.codeword_len;
    uint32_t m = n - k;
    uint8_t **parity_matrix = g_ldpc.parity_matrix;
    uint32_t cur = 0;
    uint32_t find;
    for (uint32_t c = m; c < n; c++)
    {
        find = 0;
        for (uint32_t r = cur; r < k; r++)
        {
            if (parity_matrix[r][c])
            {
                if (r != cur)
                {
                    swap_row(parity_matrix[cur], parity_matrix[r], n);
                }
                find = 1;
                break;
            }
        }

        assert(find);

        for (uint32_t r = 0; r < k; r++)
        {
            if (parity_matrix[r][c] && r != cur)
            {
                xor_to_row(parity_matrix[cur], parity_matrix[r], n);
            }
        }

        cur++;
    }
}

void init_parity_matrix(PARITY_INFO_T *parity_info)
{

    uint32_t row = parity_info->row * parity_info->blk;
    uint32_t col = parity_info->col * parity_info->blk;

    uint8_t **parity_matrix = malloc(row * sizeof(uint8_t *));
    for (uint32_t i = 0; i < row; i++)
    {
        parity_matrix[i] = malloc(col * sizeof(uint8_t));
        memset(parity_matrix[i], 0, col * sizeof(uint8_t));
    }

    uint32_t block = parity_info->blk;
    uint32_t r, c;
    int shift;
    for (uint32_t i = 0; i < parity_info->row; i++)
    {
        r = i * block;
        for (uint32_t j = 0; j < parity_info->col; j++)
        {
            shift = parity_info->build_matirx[i][j];
            if (shift >= 0)
            {
                c = j * block;
                for (uint32_t k = 0; k < block; k++)
                {
                    uint32_t s = (shift + k) % block;
                    parity_matrix[r + k][c + s] = 1;
                }
            }
        }
    }

    g_ldpc.data_len = col - row;
    g_ldpc.parity_len = row;
    g_ldpc.codeword_len = col;
    g_ldpc.code_rate = (g_ldpc.data_len) / (1.0 * g_ldpc.codeword_len);
    g_ldpc.parity_matrix = parity_matrix;
}

void init_sparse_info(PARITY_INFO_T *parity_info)
{
    uint32_t *row_one = malloc(parity_info->row * sizeof(uint32_t));
    memset(row_one, 0, parity_info->row * sizeof(uint32_t));
    for (uint32_t i = 0; i < parity_info->row; i++)
    {
        for (uint32_t j = 0; j < parity_info->col; j++)
        {
            if (parity_info->build_matirx[i][j] >= 0)
            {
                row_one[i]++;
            }
        }
    }

    SPARSE_INFO_T **sparse_row = malloc(g_ldpc.parity_len * sizeof(SPARSE_INFO_T *));
    for (uint32_t i = 0; i < parity_info->row; i++)
    {
        for (uint32_t j = 0; j < parity_info->blk; j++)
        {
            uint32_t index = i * parity_info->blk + j;
            sparse_row[index] = malloc(sizeof(SPARSE_INFO_T) + row_one[i] * sizeof(uint32_t));
            sparse_row[index]->len = 0;
        }
    }

    for (uint32_t i = 0; i < g_ldpc.parity_len; i++)
    {
        for (uint32_t j = 0; j < g_ldpc.codeword_len; j++)
        {
            if (g_ldpc.parity_matrix[i][j])
            {
                sparse_row[i]->pos[sparse_row[i]->len] = j;
                sparse_row[i]->len++;
            }
        }
    }

    simplify_parity_matrix();
    SPARSE_INFO_T **sparse_gen = malloc(g_ldpc.parity_len * sizeof(SPARSE_INFO_T *));
    for (uint32_t i = 0; i < g_ldpc.parity_len; i++)
    {
        uint32_t cnt = 0;
        for (uint32_t j = 0; j < g_ldpc.data_len; j++)
        {
            if (g_ldpc.parity_matrix[i][j])
            {
                cnt++;
            }
        }

        sparse_gen[i] = malloc(sizeof(SPARSE_INFO_T) + cnt * sizeof(uint32_t));
        sparse_gen[i]->len = 0;

        for (uint32_t j = 0; j < g_ldpc.data_len; j++)
        {
            if (g_ldpc.parity_matrix[i][j])
            {
                sparse_gen[i]->pos[sparse_gen[i]->len] = j;
                sparse_gen[i]->len++;
            }
        }
    }

    free(row_one);
    for (uint32_t i = 0; i < g_ldpc.parity_len; i++)
    {
        free(g_ldpc.parity_matrix[i]);
    }
    free(g_ldpc.parity_matrix);
    g_ldpc.parity_matrix = NULL;
    g_ldpc.sparse_row = sparse_row;
    g_ldpc.sparse_gen = sparse_gen;
}

void reset_cv_matrix(void)
{
    for (uint32_t i = 0; i < g_ldpc.parity_len; i++)
    {
        for (uint32_t j = 0; j < g_ldpc.sparse_row[i]->len; j++)
        {
            g_ldpc.cv_matrix[i][j] = 0;
        }
    }
}

void init_cv_matrix(void)
{
    double **cv_matrix = malloc(g_ldpc.parity_len * sizeof(double *));
    for (uint32_t i = 0; i < g_ldpc.parity_len; i++)
    {
        cv_matrix[i] = malloc(g_ldpc.sparse_row[i]->len * sizeof(double));
    }
    g_ldpc.cv_matrix = cv_matrix;
}

void ldpc_init(void *info)
{
    if (info == NULL)
    {
        return;
    }

    PARITY_INFO_T *parity_info = (PARITY_INFO_T *) info;
    init_parity_matrix(parity_info);
    init_sparse_info(parity_info);
    init_cv_matrix();
}

void gen_random_data(uint8_t *data)
{
    for (uint32_t i = 0; i < g_ldpc.data_len; i++)
    {
        data[i] = gen_rand() & 1;
    }
}

void ldpc_encode(const uint8_t *data, uint8_t *code_word)
{
    for (uint32_t i = 0; i < g_ldpc.data_len; i++)
    {
        code_word[i] = data[i];
    }

    uint32_t xor;
    for (uint32_t i = 0; i < g_ldpc.parity_len; i++)
    {
        xor = 0;
        for (uint32_t j = 0; j < g_ldpc.sparse_gen[i]->len; j++)
        {
            xor ^= data[g_ldpc.sparse_gen[i]->pos[j]];
        }
        code_word[g_ldpc.data_len + i] = xor;
    }
}

void ldpc_check_codeword(uint8_t *code_word)
{
    uint32_t xor;
    uint32_t pos;
    SPARSE_INFO_T **sr = g_ldpc.sparse_row;
    for (uint32_t i = 0; i < g_ldpc.parity_len; i++)
    {
        xor = 0;
        for (uint32_t j = 0; j < sr[i]->len; j++)
        {
            pos = sr[i]->pos[j];
            xor ^= code_word[pos];
        }
        assert(xor == 0);
    }
}

void ldpc_received(const uint8_t *code_word, double *received, double snr_db)
{
    double sigma = get_sigma_with_snr(snr_db, g_ldpc.code_rate);
    double sigma_2 = sigma * sigma;

    for (uint32_t i = 0; i < g_ldpc.codeword_len; i++)
    {
        // bpsk modulation
        received[i] = 2 * code_word[i] - 1;

        // add noise
        double noise = get_awgn_noise(sigma);
        received[i] = received[i] + noise;

        // log likelihood ratio
        received[i] = (-2 * received[i]) / sigma_2;
    }
}

void ldpc_decode(uint8_t *decode_data, const double *received, uint32_t iter, DECODE_METHOD method)
{
    uint32_t pos;
    double product, lambda;
    SPARSE_INFO_T **sr = g_ldpc.sparse_row;
    double **cv_matrix = g_ldpc.cv_matrix;
    double *var_node = malloc(g_ldpc.codeword_len * sizeof(double));
    double *var_node_temp = malloc(g_ldpc.codeword_len * sizeof(double));
    double *check_node = malloc(g_ldpc.parity_len * sizeof(double));
    memset(var_node, 0, g_ldpc.codeword_len * sizeof(double));
    reset_cv_matrix();

    for (uint32_t it = 0; it < iter; it++)
    {
        for (uint32_t c_i = 0; c_i < g_ldpc.parity_len; c_i++)
        {
            product = 1;
            for (uint32_t v_i = 0; v_i < sr[c_i]->len; v_i++)
            {
                pos = sr[c_i]->pos[v_i];
                lambda = received[pos] + var_node[pos] - cv_matrix[c_i][v_i];
                product = product * tanh(lambda / 2);
            }

            check_node[c_i] = product;
        }

        memset(var_node_temp, 0, g_ldpc.codeword_len * sizeof(double));

        for (uint32_t c_i = 0; c_i < g_ldpc.parity_len; c_i++)
        {
            for (uint32_t v_i = 0; v_i < sr[c_i]->len; v_i++)
            {
                pos = sr[c_i]->pos[v_i];
                lambda = received[pos] + var_node[pos] - cv_matrix[c_i][v_i];
                cv_matrix[c_i][v_i] = 2 * stable_atanh(check_node[c_i] / tanh(lambda / 2));
                var_node_temp[pos] += cv_matrix[c_i][v_i];
            }
        }

        memcpy(var_node, var_node_temp, g_ldpc.codeword_len * sizeof(double));
    }

    for (uint32_t i = 0; i < g_ldpc.data_len; i++)
    {
        decode_data[i] = (received[i] + var_node[i]) > 0 ? 0 : 1;
    }

    free(var_node);
    free(var_node_temp);
    free(check_node);
}

double ldpc_simulation(double snr_db, uint32_t iter, DECODE_METHOD method)
{
    double bit_error_rate = 0;
    uint32_t times = 100;
    uint32_t error_cnt = 0;
    uint8_t *data = malloc(g_ldpc.data_len * sizeof(uint8_t));
    uint8_t *code_word = malloc(g_ldpc.codeword_len * sizeof(uint8_t));
    uint8_t *decode_data = malloc(g_ldpc.data_len * sizeof(uint8_t));
    double *received = malloc(g_ldpc.codeword_len * sizeof(double));

    for (uint32_t i = 0; i < times; i++)
    {
        gen_random_data(data);
        ldpc_encode(data, code_word);
        ldpc_received(code_word, received, snr_db);
        ldpc_decode(decode_data, received, iter, method);

        for (uint32_t j = 0; j < g_ldpc.data_len; j++)
        {
            error_cnt += data[j] ^ decode_data[j];
        }
    }

    bit_error_rate = error_cnt / (g_ldpc.data_len * times * 1.0);

    free(data);
    free(code_word);
    free(decode_data);
    free(received);
    return bit_error_rate;
}

void ldpc_release(void)
{
    for (uint32_t i = 0; i < g_ldpc.parity_len; i++)
    {
        free(g_ldpc.sparse_row[i]);
        free(g_ldpc.sparse_gen[i]);
        free(g_ldpc.cv_matrix[i]);
    }

    free(g_ldpc.sparse_row);
    free(g_ldpc.sparse_gen);
    free(g_ldpc.cv_matrix);

}