#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include "ldpc.h"
#include "awgn.h"
#include "mt_random.h"
#include "read_parity.h"

#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#define NORM_ACTOR  0.79

typedef struct LDPC_INFO_T LDPC_INFO_T;
typedef void (*DECODE_ALGORITHM)(LDPC_INFO_T *p_ldpc_info);

typedef struct SPARSE_INFO_T_
{
    uint32_t len;
    uint32_t pos[0];
} SPARSE_INFO_T;

typedef struct MIN_SUM_INFO_T_
{
    double first_min;
    double second_min;
    int sign;
} MIN_SUM_INFO_T;

struct LDPC_INFO_T
{
    uint32_t block;
    uint32_t data_len;
    uint32_t parity_len;
    uint32_t codeword_len;
    uint32_t early_term;
    uint32_t max_iter;
    DECODE_METHOD decode_method;
    double snr_db;
    double code_rate;
    SPARSE_INFO_T **sparse_gen;
    SPARSE_INFO_T **sparse_row;
    double *var_node;
    double *var_node_temp;
    void *check_node;
    double **cv_matrix;
    uint8_t *data;
    uint8_t *code_word;
    uint8_t *decode_data;
    double *received;
    uint8_t **parity_matrix;
    DECODE_ALGORITHM dec_algorithm[MAX_ALGORITHM];
};

void spa_algorithm(LDPC_INFO_T *p_ldpc_info);
void layer_spa_algorithm(LDPC_INFO_T *p_ldpc_info);
void ms_algorithm(LDPC_INFO_T *p_ldpc_info);
void layer_ms_algorithm(LDPC_INFO_T *p_ldpc_info);

inline int sgn(double x)
{
    if (x < 0) return -1;
    if (x > 0) return 1;
    return 0;
}

inline double stable_atanh(double t)
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

void simplify_parity_matrix(LDPC_INFO_T *p_ldpc_info)
{
    uint32_t k = p_ldpc_info->parity_len;
    uint32_t n = p_ldpc_info->codeword_len;
    uint32_t m = n - k;
    uint8_t **parity_matrix = p_ldpc_info->parity_matrix;
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

void init_parity_matrix(PARITY_INFO_T *parity_info, LDPC_INFO_T *p_ldpc_info)
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

    p_ldpc_info->data_len = col - row;
    p_ldpc_info->parity_len = row;
    p_ldpc_info->codeword_len = col;
    p_ldpc_info->code_rate = (p_ldpc_info->data_len) / (1.0 * p_ldpc_info->codeword_len);
    p_ldpc_info->parity_matrix = parity_matrix;
}

void init_sparse_info(PARITY_INFO_T *parity_info, LDPC_INFO_T *p_ldpc_info)
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

    SPARSE_INFO_T **sparse_row = malloc(p_ldpc_info->parity_len * sizeof(SPARSE_INFO_T *));
    for (uint32_t i = 0; i < parity_info->row; i++)
    {
        for (uint32_t j = 0; j < parity_info->blk; j++)
        {
            uint32_t index = i * parity_info->blk + j;
            sparse_row[index] = malloc(sizeof(SPARSE_INFO_T) + row_one[i] * sizeof(uint32_t));
            sparse_row[index]->len = 0;
        }
    }

    for (uint32_t i = 0; i < p_ldpc_info->parity_len; i++)
    {
        for (uint32_t j = 0; j < p_ldpc_info->codeword_len; j++)
        {
            if (p_ldpc_info->parity_matrix[i][j])
            {
                sparse_row[i]->pos[sparse_row[i]->len] = j;
                sparse_row[i]->len++;
            }
        }
    }

    simplify_parity_matrix(p_ldpc_info);
    SPARSE_INFO_T **sparse_gen = malloc(p_ldpc_info->parity_len * sizeof(SPARSE_INFO_T *));
    for (uint32_t i = 0; i < p_ldpc_info->parity_len; i++)
    {
        uint32_t cnt = 0;
        for (uint32_t j = 0; j < p_ldpc_info->data_len; j++)
        {
            if (p_ldpc_info->parity_matrix[i][j])
            {
                cnt++;
            }
        }

        sparse_gen[i] = malloc(sizeof(SPARSE_INFO_T) + cnt * sizeof(uint32_t));
        sparse_gen[i]->len = 0;

        for (uint32_t j = 0; j < p_ldpc_info->data_len; j++)
        {
            if (p_ldpc_info->parity_matrix[i][j])
            {
                sparse_gen[i]->pos[sparse_gen[i]->len] = j;
                sparse_gen[i]->len++;
            }
        }
    }

    free(row_one);
    for (uint32_t i = 0; i < p_ldpc_info->parity_len; i++)
    {
        free(p_ldpc_info->parity_matrix[i]);
    }
    free(p_ldpc_info->parity_matrix);
    p_ldpc_info->parity_matrix = NULL;
    p_ldpc_info->sparse_row = sparse_row;
    p_ldpc_info->sparse_gen = sparse_gen;
}

void reset_cv_matrix(LDPC_INFO_T *p_ldpc_info)
{
    for (uint32_t i = 0; i < p_ldpc_info->parity_len; i++)
    {
        for (uint32_t j = 0; j < p_ldpc_info->sparse_row[i]->len; j++)
        {
            p_ldpc_info->cv_matrix[i][j] = 0;
        }
    }
}

void init_cv_matrix(LDPC_INFO_T *p_ldpc_info)
{
    double **cv_matrix = malloc(p_ldpc_info->parity_len * sizeof(double *));
    for (uint32_t i = 0; i < p_ldpc_info->parity_len; i++)
    {
        cv_matrix[i] = malloc(p_ldpc_info->sparse_row[i]->len * sizeof(double));
    }
    p_ldpc_info->cv_matrix = cv_matrix;
}

void init_decode_algorithm(LDPC_INFO_T *p_ldpc_info)
{
    p_ldpc_info->dec_algorithm[SPA_ALGORITHM] = spa_algorithm;
    p_ldpc_info->dec_algorithm[LAYERED_SPA_ALGORITHM] = layer_spa_algorithm;
    p_ldpc_info->dec_algorithm[MS_ALGORITHM] = ms_algorithm;
    p_ldpc_info->dec_algorithm[LAYERED_MS_ALGORITHM] = layer_ms_algorithm;
}

LDPC_INFO_T *ldpc_init(void *info)
{
    if (info == NULL)
    {
        return NULL;
    }

    LDPC_INFO_T *p_ldpc_info = malloc(sizeof(LDPC_INFO_T));
    PARITY_INFO_T *parity_info = (PARITY_INFO_T *) info;
    init_parity_matrix(parity_info, p_ldpc_info);
    init_sparse_info(parity_info, p_ldpc_info);
    init_cv_matrix(p_ldpc_info);
    init_decode_algorithm(p_ldpc_info);
    p_ldpc_info->block = parity_info->blk;
    p_ldpc_info->var_node = malloc(p_ldpc_info->codeword_len * sizeof(double));
    p_ldpc_info->var_node_temp = malloc(p_ldpc_info->codeword_len * sizeof(double));
    uint32_t size = MAX(sizeof(double), sizeof(MIN_SUM_INFO_T));
    p_ldpc_info->check_node = malloc(p_ldpc_info->parity_len * size);
    printf("%d, %d, %d, %.3f\n", p_ldpc_info->parity_len, p_ldpc_info->data_len, p_ldpc_info->codeword_len, p_ldpc_info->code_rate);
    return p_ldpc_info;
}

void gen_random_data(LDPC_INFO_T *p_ldpc_info)
{
    for (uint32_t i = 0; i < p_ldpc_info->data_len; i++)
    {
        p_ldpc_info->data[i] = gen_rand() & 1;
    }
}

void ldpc_encode(LDPC_INFO_T *p_ldpc_info)
{
    for (uint32_t i = 0; i < p_ldpc_info->data_len; i++)
    {
        p_ldpc_info->code_word[i] = p_ldpc_info->data[i];
    }

    uint32_t xor;
    for (uint32_t i = 0; i < p_ldpc_info->parity_len; i++)
    {
        xor = 0;
        for (uint32_t j = 0; j < p_ldpc_info->sparse_gen[i]->len; j++)
        {
            xor ^= p_ldpc_info->data[p_ldpc_info->sparse_gen[i]->pos[j]];
        }
        p_ldpc_info->code_word[p_ldpc_info->data_len + i] = xor;
    }
}

int ldpc_check_codeword(LDPC_INFO_T *p_ldpc_info, const uint8_t *code_word)
{
    uint32_t xor;
    SPARSE_INFO_T **sr = p_ldpc_info->sparse_row;
    for (uint32_t i = 0; i < p_ldpc_info->parity_len; i++)
    {
        xor = 0;
        for (uint32_t j = 0; j < sr[i]->len; j++)
        {
            xor ^= code_word[sr[i]->pos[j]];
        }

        if (xor)
        {
            return 0;
        }
    }

    return 1;
}

void ldpc_received(LDPC_INFO_T *p_ldpc_info)
{
    double sigma = get_sigma_with_snr(p_ldpc_info->snr_db, p_ldpc_info->code_rate);
    double sigma_2 = sigma * sigma;

    for (uint32_t i = 0; i < p_ldpc_info->codeword_len; i++)
    {
        // bpsk modulation
        p_ldpc_info->received[i] = 2 * p_ldpc_info->code_word[i] - 1;

        // add noise
        p_ldpc_info->received[i] = p_ldpc_info->received[i] + get_awgn_noise(sigma);

        // log likelihood ratio
        p_ldpc_info->received[i] = (-2 * p_ldpc_info->received[i]) / sigma_2;
    }
}

void ldpc_decode(LDPC_INFO_T *p_ldpc_info)
{
    p_ldpc_info->dec_algorithm[p_ldpc_info->decode_method](p_ldpc_info);
}

void ldpc_set_config(struct LDPC_INFO_T *p_ldpc_info, double snr_db, uint32_t iter, uint32_t early_term, DECODE_METHOD method)
{
    p_ldpc_info->snr_db = snr_db;
    p_ldpc_info->max_iter = iter;
    p_ldpc_info->early_term = early_term;
    p_ldpc_info->decode_method = method;
}

ERROR_RATE_T ldpc_simulation(LDPC_INFO_T *p_ldpc_info)
{
    ERROR_RATE_T error_rate;
    uint32_t times = 80000;
    uint32_t error_cnt = 0;
    uint32_t block_err_cnt = 0;
    uint32_t tmp;
    p_ldpc_info->data = malloc(p_ldpc_info->data_len * sizeof(uint8_t));
    p_ldpc_info->code_word = malloc(p_ldpc_info->codeword_len * sizeof(uint8_t));
    p_ldpc_info->decode_data = malloc(p_ldpc_info->codeword_len * sizeof(uint8_t));
    p_ldpc_info->received = malloc(p_ldpc_info->codeword_len * sizeof(double));
    memset(p_ldpc_info->data, 0, p_ldpc_info->data_len * sizeof(uint8_t));
    memset(p_ldpc_info->code_word, 0, p_ldpc_info->codeword_len * sizeof(uint8_t));

    for (uint32_t i = 0; i < times; i++)
    {
        gen_random_data(p_ldpc_info);
        ldpc_encode(p_ldpc_info);
        ldpc_received(p_ldpc_info);
        ldpc_decode(p_ldpc_info);

        tmp = error_cnt;
        for (uint32_t j = 0; j < p_ldpc_info->data_len; j++)
        {
            error_cnt += p_ldpc_info->data[j] ^ p_ldpc_info->decode_data[j];
        }
        if (error_cnt > tmp)
        {
            block_err_cnt++;
        }
    }

    error_rate.bit_error_rate = error_cnt / (p_ldpc_info->data_len * times * 1.0);
    error_rate.block_error_rate = block_err_cnt / (times * 1.0);

    free(p_ldpc_info->data);
    free(p_ldpc_info->code_word);
    free(p_ldpc_info->decode_data);
    free(p_ldpc_info->received);
    return error_rate;
}

void ldpc_release(LDPC_INFO_T *p_ldpc_info)
{
    for (uint32_t i = 0; i < p_ldpc_info->parity_len; i++)
    {
        free(p_ldpc_info->sparse_row[i]);
        free(p_ldpc_info->sparse_gen[i]);
        free(p_ldpc_info->cv_matrix[i]);
    }

    free(p_ldpc_info->sparse_row);
    free(p_ldpc_info->sparse_gen);
    free(p_ldpc_info->cv_matrix);
    free(p_ldpc_info->var_node);
    free(p_ldpc_info->var_node_temp);
    free(p_ldpc_info->check_node);
}


void spa_algorithm(LDPC_INFO_T *p_ldpc_info)
{
    uint32_t pos;
    double product, lambda;
    SPARSE_INFO_T **sr = p_ldpc_info->sparse_row;
    double **cv_matrix = p_ldpc_info->cv_matrix;
    double *var_node = p_ldpc_info->var_node;
    double *var_node_temp = p_ldpc_info->var_node_temp;
    double *check_node = (double *) p_ldpc_info->check_node;
    memset(var_node, 0, p_ldpc_info->codeword_len * sizeof(double));
    reset_cv_matrix(p_ldpc_info);

    for (uint32_t it = 0; it < p_ldpc_info->max_iter; it++)
    {
        for (uint32_t c_i = 0; c_i < p_ldpc_info->parity_len; c_i++)
        {
            product = 1;
            for (uint32_t v_i = 0; v_i < sr[c_i]->len; v_i++)
            {
                pos = sr[c_i]->pos[v_i];
                lambda = p_ldpc_info->received[pos] + var_node[pos] - cv_matrix[c_i][v_i];
                product = product * tanh(lambda / 2);
            }

            check_node[c_i] = product;
        }

        memset(var_node_temp, 0, p_ldpc_info->codeword_len * sizeof(double));

        for (uint32_t c_i = 0; c_i < p_ldpc_info->parity_len; c_i++)
        {
            for (uint32_t v_i = 0; v_i < sr[c_i]->len; v_i++)
            {
                pos = sr[c_i]->pos[v_i];
                lambda = p_ldpc_info->received[pos] + var_node[pos] - cv_matrix[c_i][v_i];
                var_node[pos] -= cv_matrix[c_i][v_i];
                cv_matrix[c_i][v_i] = 2 * stable_atanh(check_node[c_i] / tanh(lambda / 2));
                var_node[pos] += cv_matrix[c_i][v_i];
            }
        }

        memcpy(var_node, var_node_temp, p_ldpc_info->codeword_len * sizeof(double));

        if (p_ldpc_info->early_term)
        {
            for (uint32_t i = 0; i < p_ldpc_info->codeword_len; i++)
            {
                p_ldpc_info->decode_data[i] = (p_ldpc_info->received[i] + var_node[i]) > 0 ? 0 : 1;
            }

            if (ldpc_check_codeword(p_ldpc_info, p_ldpc_info->decode_data))
            {
                return;
            }
        }
    }

    if (!p_ldpc_info->early_term)
    {
        for (uint32_t i = 0; i < p_ldpc_info->data_len; i++)
        {
            p_ldpc_info->decode_data[i] = (p_ldpc_info->received[i] + var_node[i]) > 0 ? 0 : 1;
        }
    }
}

void layer_spa_algorithm(LDPC_INFO_T *p_ldpc_info)
{
    uint32_t pos;
    uint32_t layer = p_ldpc_info->parity_len / p_ldpc_info->block;
    double product, lambda;
    SPARSE_INFO_T **sr = p_ldpc_info->sparse_row;
    double **cv_matrix = p_ldpc_info->cv_matrix;
    double *var_node = p_ldpc_info->var_node;
    double *check_node = (double *) p_ldpc_info->check_node;
    memset(var_node, 0, p_ldpc_info->codeword_len * sizeof(double));
    reset_cv_matrix(p_ldpc_info);

    for (uint32_t it = 0; it < p_ldpc_info->max_iter; it++)
    {
        uint32_t l_offset = 0;

        for (uint32_t l = 0; l < layer; l++)
        {
            for (uint32_t c_i = l_offset; c_i < l_offset + p_ldpc_info->block; c_i++)
            {
                product = 1;

                for (uint32_t v_i = 0; v_i < sr[c_i]->len; v_i++)
                {
                    pos = sr[c_i]->pos[v_i];
                    lambda = p_ldpc_info->received[pos] + var_node[pos] - cv_matrix[c_i][v_i];
                    product = product * tanh(lambda / 2);
                }

                check_node[c_i] = product;
            }

            for (uint32_t c_i = l_offset; c_i < l_offset + p_ldpc_info->block; c_i++)
            {
                for (uint32_t v_i = 0; v_i < sr[c_i]->len; v_i++)
                {
                    pos = sr[c_i]->pos[v_i];
                    lambda = p_ldpc_info->received[pos] + var_node[pos] - cv_matrix[c_i][v_i];
                    var_node[pos] -= cv_matrix[c_i][v_i];
                    cv_matrix[c_i][v_i] = 2 * stable_atanh(check_node[c_i] / tanh(lambda / 2));
                    var_node[pos] += cv_matrix[c_i][v_i];
                }
            }

            l_offset += p_ldpc_info->block;
        }

        if (p_ldpc_info->early_term)
        {
            for (uint32_t i = 0; i < p_ldpc_info->codeword_len; i++)
            {
                p_ldpc_info->decode_data[i] = (p_ldpc_info->received[i] + var_node[i]) > 0 ? 0 : 1;
            }

            if (ldpc_check_codeword(p_ldpc_info, p_ldpc_info->decode_data))
            {
                return;
            }
        }

    }

    if (!p_ldpc_info->early_term)
    {
        for (uint32_t i = 0; i < p_ldpc_info->data_len; i++)
        {
            p_ldpc_info->decode_data[i] = (p_ldpc_info->received[i] + var_node[i]) > 0 ? 0 : 1;
        }
    }
}

void ms_algorithm(LDPC_INFO_T *p_ldpc_info)
{
    uint32_t pos;
    int sign;
    double lambda, abs_lambda, fm, sm;
    SPARSE_INFO_T **sr = p_ldpc_info->sparse_row;
    double **cv_matrix = p_ldpc_info->cv_matrix;
    double *var_node = p_ldpc_info->var_node;
    double *var_node_temp = p_ldpc_info->var_node_temp;
    MIN_SUM_INFO_T *check_node = (MIN_SUM_INFO_T *) p_ldpc_info->check_node;

    memset(var_node, 0, p_ldpc_info->codeword_len * sizeof(double));
    reset_cv_matrix(p_ldpc_info);

    for (uint32_t it = 0; it < p_ldpc_info->max_iter; it++)
    {
        for (uint32_t c_i = 0; c_i < p_ldpc_info->parity_len; c_i++)
        {
            fm = DBL_MAX;
            sm = DBL_MAX;
            sign = 1;
            for (uint32_t v_i = 0; v_i < sr[c_i]->len; v_i++)
            {
                pos = sr[c_i]->pos[v_i];
                lambda = p_ldpc_info->received[pos] + var_node[pos] - cv_matrix[c_i][v_i];
                abs_lambda = fabs(lambda);
                if (abs_lambda < fm)
                {
                    sm = fm;
                    fm = abs_lambda;
                }
                else if (abs_lambda < sm)
                {
                    sm = abs_lambda;
                }
                sign *= sgn(lambda);
            }

            check_node[c_i] = (MIN_SUM_INFO_T) {fm, sm, sign};
        }

        memset(var_node_temp, 0, p_ldpc_info->codeword_len * sizeof(double));

        for (uint32_t c_i = 0; c_i < p_ldpc_info->parity_len; c_i++)
        {
            for (uint32_t v_i = 0; v_i < sr[c_i]->len; v_i++)
            {
                pos = sr[c_i]->pos[v_i];
                lambda = p_ldpc_info->received[pos] + var_node[pos] - cv_matrix[c_i][v_i];
                abs_lambda = fabs(lambda);
                cv_matrix[c_i][v_i] = (check_node[c_i].sign * sgn(lambda)) * NORM_ACTOR *
                        (abs_lambda == check_node[c_i].first_min ? check_node[c_i].second_min
                        : check_node[c_i].first_min);
                var_node_temp[pos] += cv_matrix[c_i][v_i];
            }
        }

        memcpy(var_node, var_node_temp, p_ldpc_info->codeword_len * sizeof(double));

        if (p_ldpc_info->early_term)
        {
            for (uint32_t i = 0; i < p_ldpc_info->codeword_len; i++)
            {
                p_ldpc_info->decode_data[i] = (p_ldpc_info->received[i] + var_node[i]) > 0 ? 0 : 1;
            }

            if (ldpc_check_codeword(p_ldpc_info, p_ldpc_info->decode_data))
            {
                return;
            }
        }
    }

    if (!p_ldpc_info->early_term)
    {
        for (uint32_t i = 0; i < p_ldpc_info->data_len; i++)
        {
            p_ldpc_info->decode_data[i] = (p_ldpc_info->received[i] + var_node[i]) > 0 ? 0 : 1;
        }
    }
}

void layer_ms_algorithm(LDPC_INFO_T *p_ldpc_info)
{
    uint32_t pos;
    int sign;
    uint32_t layer = p_ldpc_info->parity_len / p_ldpc_info->block;
    double lambda, abs_lambda, fm, sm;
    SPARSE_INFO_T **sr = p_ldpc_info->sparse_row;
    double **cv_matrix = p_ldpc_info->cv_matrix;
    double *var_node = p_ldpc_info->var_node;
    MIN_SUM_INFO_T *check_node = (MIN_SUM_INFO_T *) p_ldpc_info->check_node;

    memset(var_node, 0, p_ldpc_info->codeword_len * sizeof(double));
    reset_cv_matrix(p_ldpc_info);

    for (uint32_t it = 0; it < p_ldpc_info->max_iter; it++)
    {
        uint32_t l_offset = 0;

        for (uint32_t l = 0; l < layer; l++)
        {
            for (uint32_t c_i = l_offset; c_i < l_offset + p_ldpc_info->block; c_i++)
            {
                fm = DBL_MAX;
                sm = DBL_MAX;
                sign = 1;
                for (uint32_t v_i = 0; v_i < sr[c_i]->len; v_i++)
                {
                    pos = sr[c_i]->pos[v_i];
                    lambda = p_ldpc_info->received[pos] + var_node[pos] - cv_matrix[c_i][v_i];
                    abs_lambda = fabs(lambda);
                    if (abs_lambda < fm)
                    {
                        sm = fm;
                        fm = abs_lambda;
                    }
                    else if (abs_lambda < sm)
                    {
                        sm = abs_lambda;
                    }
                    sign *= sgn(lambda);
                }

                check_node[c_i] = (MIN_SUM_INFO_T) {fm, sm, sign};
            }

            for (uint32_t c_i = l_offset; c_i < l_offset + p_ldpc_info->block; c_i++)
            {
                for (uint32_t v_i = 0; v_i < sr[c_i]->len; v_i++)
                {
                    pos = sr[c_i]->pos[v_i];
                    lambda = p_ldpc_info->received[pos] + var_node[pos] - cv_matrix[c_i][v_i];
                    abs_lambda = fabs(lambda);
                    var_node[pos] -= cv_matrix[c_i][v_i];
                    cv_matrix[c_i][v_i] = (check_node[c_i].sign * sgn(lambda)) * NORM_ACTOR *
                            (abs_lambda == check_node[c_i].first_min ? check_node[c_i].second_min
                            : check_node[c_i].first_min);
                    var_node[pos] += cv_matrix[c_i][v_i];
                }
            }

            l_offset += p_ldpc_info->block;
        }

        if (p_ldpc_info->early_term)
        {
            for (uint32_t i = 0; i < p_ldpc_info->codeword_len; i++)
            {
                p_ldpc_info->decode_data[i] = (p_ldpc_info->received[i] + var_node[i]) > 0 ? 0 : 1;
            }

            if (ldpc_check_codeword(p_ldpc_info, p_ldpc_info->decode_data))
            {
                return;
            }
        }
    }

    if (!p_ldpc_info->early_term)
    {
        for (uint32_t i = 0; i < p_ldpc_info->data_len; i++)
        {
            p_ldpc_info->decode_data[i] = (p_ldpc_info->received[i] + var_node[i]) > 0 ? 0 : 1;
        }
    }
}