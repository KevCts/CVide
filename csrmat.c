#include "csrmat.h" 
#include "coomat.h"
#include <stdio.h>

csrmat* coo_to_csr(coomat* coo) {
    csrmat* res = malloc(sizeof(csrmat));

    res->size_i = coo->size_i;
    res->size_j = coo->size_j;
    res->nnn = 0;
    res->ia = calloc(res->size_i + 1, sizeof(size_t));
    res->ja = malloc(sizeof(size_t));
    res->values = malloc(sizeof(size_t));

    for (size_t i = 0; i < res->size_i; i++) {
        for (size_t j = 0; j < res->size_j; j++) {
            if (coomat_read_value(coo, i, j) != 0) {
                res->nnn++;
                res->ja = realloc(res->ja, res->nnn*sizeof(size_t));
                res->ja[res->nnn - 1] = j;
                res->values = realloc(res->values, res->nnn*sizeof(double));
                res->values[res->nnn - 1] = coomat_read_value(coo, i, j);
            }
        }
        res->ia[i+1] = res->nnn;
    }

    free_coomat(coo);

    return res;
}

coomat* prod_csr_coo(csrmat* a, coomat* b) {
    coomat* res = init_coomat(a->size_i, b->size_j);

    for (int jb = 0; jb < b->size_j; jb++) {
        for (int i = 0; i < a->size_i; i++) {
            for(int j = a->ia[i]; j < a->ia[i+1]; j++) {
                coomat_set_value(res, i, jb, coomat_read_value(res, i, jb) + coomat_read_value(b, a->ja[j], jb) * a->values[j]);
            }
        }
    }

    free_csrmat(a);
    free_coomat(b);

    return res;
}

void free_csrmat(csrmat* matrix) {
    free(matrix->ia);
    free(matrix->ja);
    free(matrix->values);
    free(matrix);
}
