#ifndef CVIDE_COOMAT
#define CVIDE_COOMAT

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

typedef struct {
    size_t size_i, size_j;
    double* values;
} coomat;

coomat* init_coomat(size_t, size_t);

void free_coomat(coomat*);

void print_coomat(coomat*);

double coomat_read_value(coomat*, size_t, size_t);

double coomat_set_value(coomat*, size_t, size_t, double);

coomat* copy_coomat(coomat*);

coomat* sum_coomat(coomat*, coomat*);

coomat* dot_coomat(double, coomat*);

bool equal_coomat(coomat*, coomat*);

double scalar_coomat(coomat*, coomat*);

coomat* prod_coomat(coomat*, coomat*);

#endif
