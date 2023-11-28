#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_parity.h"

PARITY_INFO_T *get_parity_matrix_info(const char *file_name)
{
    FILE *file;
    char *endptr;
    int i, j;
    int val;
    char line[1024];
    char *info_name[3] = {"block", "row", "column"};
    int mat_info[3];
    int blk, col, row;
    PARITY_INFO_T *info = NULL;

    file = fopen(file_name, "r");
    if (file == NULL)
    {
        printf("%s without open\n", file_name);
        return info;
    }

    if (fgets(line, sizeof(line), file) == NULL)
    {
        printf("%s matrix info invalid\n", file_name);
        goto EXIT;
    }

    char *token = strtok(line, " ,");
    for (i = 0; i < 3; i++)
    {
        if (token == NULL)
        {
            printf("%s invalid\n", info_name[i]);
            goto EXIT;
        }

        errno = 0;
        val = strtol(token, &endptr, 10);

        if (endptr == token || errno != 0 || val <= 0)
        {
            printf("%s invalid\n", info_name[i]);
            goto EXIT;
        }

        mat_info[i] = val;
        token = strtok(NULL, " ,");
    }

    blk = mat_info[0];
    row = mat_info[1];
    col = mat_info[2];

    int **data = (int **) malloc(row * sizeof(int *));
    for (i = 0; i < row; i++)
    {
        data[i] = (int *) malloc(col * sizeof(int));
    }

    for (i = 0; i < row; i++)
    {
        if (fgets(line, sizeof(line), file) == NULL)
        {
            printf("matrix invalid at row %d\n", i);
            goto FREE;
        }

        token = strtok(line, " ,");
        for (j = 0; j < col; j++)
        {
            if (token == NULL)
            {
                printf("matrix invalid at row %d, col %d\n", i, j);
                goto FREE;
            }

            errno = 0;
            val = strtol(token, &endptr, 10);

            if (endptr == token || errno != 0 || val < -1 || val >= blk)
            {
                printf("matrix invalid at row %d, col %d\n", i, j);
                goto FREE;
            }

            data[i][j] = val;
            token = strtok(NULL, " ,");
        }
    }
    
    info = malloc(sizeof(PARITY_INFO_T));
    info->blk = blk;
    info->row = row;
    info->col = col;
    info->build_matirx = (int **) malloc(row * sizeof(int *));
    for (i = 0; i < row; i++)
    {
        info->build_matirx[i] = (int *) malloc(col * sizeof(int));
        for (j = 0; j < col; j++)
        {
            info->build_matirx[i][j] = data[i][j];
        }
    }

    FREE:
    for (i = 0; i < row; i++)
    {
        free(data[i]);
    }
    free(data);

    EXIT:
    fclose(file);
    return info;
}

void release_parity_info(PARITY_INFO_T *parity_info)
{
    if (parity_info == NULL)
    {
        return;
    }

    for (uint32_t i = 0; i < parity_info->row; i++)
    {
        free(parity_info->build_matirx[i]);
    }

    free(parity_info->build_matirx);
    free(parity_info);
}

