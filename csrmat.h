#ifndef CVIDE_CSRMAT_H
#define CVIDE_CSRMAT_H

#include "coomat.h"

typedef struct {
    size_t size_i;
    size_t size_j;
    double* values;
    size_t* ia;
    size_t* ja;
} csrmat;

csrmat* coo_to_csr(coomat*);

void free_csrmat(csrmat*);

#endif
