#include "csrmat.h" 
#include "coomat.h"
#include <stdio.h>

csrmat* coo_to_csr(coomat* coo) {
    csrmat* res = malloc(sizeof(csrmat));

    res->size_i = coo->size_i;
    res->size_j = coo->size_j;
    res->ia = calloc(res->size_i + 1, sizeof(size_t));
    res->ja = malloc(sizeof(size_t));
    res->values = malloc(sizeof(size_t));

    size_t nnn = 0;
    
    for (size_t i = 0; i < res->size_i; i++) {
        for (size_t j = 0; j < res->size_j; j++) {
            if (coomat_read_value(coo, i, j) != 0) {
                nnn++;
                res->ja = realloc(res->ja, nnn*sizeof(size_t));
                res->ja[nnn - 1] = j;
                res->values = realloc(res->values, nnn*sizeof(double));
                res->values[nnn - 1] = coomat_read_value(coo, i, j);
            }
        }
        res->ia[i+1] = nnn;
        printf("%ld \n", res->ia[i+ 1]);
    }

    free_coomat(coo);

    return res;
}

void free_csrmat(csrmat* matrix) {
    free(matrix->ia);
    free(matrix->ja);
    free(matrix->values);
    free(matrix);
}
