#include "coomat.h"

#include <stdio.h>
#include <stdlib.h>

coomat init_coomat(size_t n, size_t m) {
    double* values = (double *) calloc((n-1) * (m-1), sizeof(double));
    coomat res = {n, m, values};
    return res;
}

void free_coomat(coomat* matrix) {
    free(matrix->values);
}

double coomat_read_value(coomat* matrix, size_t i, size_t j) {
    if (i >= matrix->size_i || j >= matrix->size_j) {
        return 0;
    } else {
        return matrix->values[i * matrix->size_i + j];
    }
}


double coomat_set_value(coomat* matrix, size_t i, size_t j, double value) {
    if (i >= matrix->size_i || j >= matrix->size_j) {
        return 0;
    } else {
        matrix->values[i * matrix->size_i + j] = value;
        return matrix->values[i * matrix->size_i + j];
    }
}

void print_coomat(coomat* matrix) {
    for (int i = 0; i < matrix->size_i; i++) {
        for (int j = 0; j < matrix->size_j; j++) {
            printf("%lf | ", coomat_read_value(matrix, i, j));
        }
        printf("\n");
    }
    free(matrix);
}

coomat copy_coomat(coomat* matrix) {
    coomat res = init_coomat(matrix->size_i, matrix->size_j);
    for (int i = 0; i < matrix->size_i; i++) {
        for (int j = 0; j < matrix->size_j; j++) {
            coomat_set_value(&res, i, j, coomat_read_value(matrix, i, j));
        }
    }
    return res;
}

coomat sum_coomat(coomat* a, coomat* b){
    coomat res = init_coomat(a->size_i, a->size_j);
    if (a->size_i == b->size_i && a->size_j == b->size_j){
        for (int i = 0; i < a->size_i; i++) {
            for (int j = 0; j < a->size_j; j++) {
                coomat_set_value(&res, i, j, coomat_read_value(a, i, j) + coomat_read_value(b, i, j));
            }
        }
    }

    free_coomat(a);
    free_coomat(b);
    return res;
}

coomat dot_coomat(double a, coomat* b){
    coomat res = init_coomat(b->size_i, b->size_j);
    for (int i = 0; i < b->size_i; i++) {
        for (int j = 0; j < b->size_j; j++) {
            coomat_set_value(&res, i, j, a * coomat_read_value(b, i, j));
        }
    }

    free_coomat(b);
    return res;
}

double scalar_coomat(coomat* a, coomat* b) {
    double res = 0;

    if (a->size_i == b->size_i && a->size_j == b->size_j){
        for (int i = 0; i < a->size_i; i++) {
            for (int j = 0; j < a->size_j; j++) {
                res += coomat_read_value(a, i, j) * coomat_read_value(b, i, j);
            }
        }
    }

    free_coomat(a);
    free_coomat(b);
    return res;
}

bool equal_coomat(coomat* a, coomat* b) {
    if (a->size_i != b->size_i && a->size_j != b->size_j){
        return false; 
    }

        for (int i = 0; i < a->size_i; i++) {
            for (int j = 0; j < a->size_j; j++) {
                if (coomat_read_value(a, i, j) != coomat_read_value(b, i, j)) {
                    return false; 
                }
            }
        }

    return true;
}

coomat prod_coomat(coomat* a, coomat* b) {
    coomat res = init_coomat(a->size_i, b->size_j);

    if (a->size_j == b->size_i) {
        for (int i = 0; i < a->size_i; i++) {
            for (int j = 0; j < b->size_j; j++) {
                double c = 0;
                for (int k = 0; k < a->size_j; k++) {
                    c += coomat_read_value(a, i, k) * coomat_read_value(b, k, j);
                }
                coomat_set_value(&res, i, j, c);
            }
        }
    }

    free_coomat(a);
    free_coomat(b);
    return res;
}
