#include "coomat.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

coomat* init_coomat(size_t n, size_t m) {
    coomat* res = malloc(sizeof(coomat));
    res->size_i = n;
    res->size_j = m;
    res->values = (double *) calloc(n * m, sizeof(double));
    return res;
}

void free_coomat(coomat* matrix) {
    if (matrix->values != NULL)
        free(matrix->values);
    free(matrix);
}

coomat_list* init_coomat_list(size_t size_i, size_t size_j) {
    coomat_list* res = malloc(sizeof(coomat_list));
    res->size_i = size_i;
    res->size_j = size_j;
    res->count = 0;
    res->elements = NULL;
    return res;
}

coomat_list* add_coomat_to_list(coomat_list* list, coomat* matrix) {
    if ((matrix->size_i == list->size_i) && (matrix->size_j == list->size_j)) {
        list->elements = realloc(list->elements, (list->count + 1) * sizeof(coomat*));
        list->elements[list->count] = matrix;
        list->count++;
    }

    return list;
}

coomat* coomat_from_list(coomat_list* list, size_t i) {
    if (i < list->count)
        return copy_coomat(list->elements[i]);
    return NULL;
}

void coomat_fun_to_list(coomat* (*f)(coomat*), coomat_list* list) {
    for (int i = 0; i < list->count; i++) {
        list->elements[i] = f(list->elements[i]);
    }
}

void double_fun_to_list(double (*f)(double), coomat_list* list) {
    for (int i = 0; i < list->count; i++) {
        list->elements[i] = fun_coomat(f, list->elements[i]);
    }
}

void free_coomat_list(coomat_list* list) {
    for (int i = 0; i < list->count; i++) {
        free_coomat(list->elements[i]);
    }

    if (list->elements != NULL)
        free(list->elements);
     
    free(list);
}

double coomat_read_value(coomat* matrix, size_t i, size_t j) {
    if (i >= matrix->size_i || j >= matrix->size_j) {
        return 0;
    } else {
        return matrix->values[i * matrix->size_j + j];
    }
}


double coomat_set_value(coomat* matrix, size_t i, size_t j, double value) {
    if (i >= matrix->size_i || j >= matrix->size_j) {
        return 0;
    } else {
        matrix->values[i * matrix->size_j + j] = value;
        return matrix->values[i * matrix->size_j + j];
    }
}

void print_coomat(coomat* matrix) {
    for (int i = 0; i < matrix->size_i; i++) {
        for (int j = 0; j < matrix->size_j; j++) {
            printf("\t%lf\t|", coomat_read_value(matrix, i, j));
        }
        printf("\n");
    }
    free_coomat(matrix);
}

coomat* copy_coomat(coomat* matrix) {
    coomat* res = init_coomat(matrix->size_i, matrix->size_j);
    for (int i = 0; i < matrix->size_i; i++) {
        for (int j = 0; j < matrix->size_j; j++) {
            coomat_set_value(res, i, j, coomat_read_value(matrix, i, j));
        }
    }
    return res;
}

coomat* fun_coomat(double (*f)(double), coomat* matrix) {
    coomat* res = init_coomat(matrix->size_i, matrix->size_j);

    for (int i = 0; i < matrix->size_i * matrix->size_j; i++) {
        res->values[i] = f(matrix->values[i]);
    }

    free_coomat(matrix);

    return res;
}

coomat* sum_coomat(coomat* a, coomat* b){
    coomat* res = init_coomat(a->size_i, a->size_j);
    if (a->size_i == b->size_i && a->size_j == b->size_j){
        for (int i = 0; i < a->size_i; i++) {
            for (int j = 0; j < a->size_j; j++) {
                coomat_set_value(res, i, j, coomat_read_value(a, i, j) + coomat_read_value(b, i, j));
            }
        }
    }

    free_coomat(a);
    free_coomat(b);
    return res;
}

coomat* dot_coomat(double a, coomat* b){
    coomat* res = init_coomat(b->size_i, b->size_j);
    for (int i = 0; i < b->size_i; i++) {
        for (int j = 0; j < b->size_j; j++) {
            coomat_set_value(res, i, j, a * coomat_read_value(b, i, j));
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

coomat* transpose_coomat(coomat* mat){
    coomat* res = init_coomat(mat->size_j, mat->size_i);

    for (int i = 0; i < mat->size_i; i++) {
        for (int j = 0; j < mat->size_j; j++) {
            coomat_set_value(res, j, i, coomat_read_value(mat, i, j));
        }
    }

    free_coomat(mat);

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

coomat* prod_coomat(coomat* a, coomat* b) {
    coomat* res = init_coomat(a->size_i, b->size_j);

    if (a->size_j == b->size_i) {
        for (int i = 0; i < a->size_i; i++) {
            for (int j = 0; j < b->size_j; j++) {
                double c = 0;
                for (int k = 0; k < a->size_j; k++) {
                    c += coomat_read_value(a, i, k) * coomat_read_value(b, k, j);
                }
                coomat_set_value(res, i, j, c);
            }
        }
    }

    free_coomat(a);
    free_coomat(b);
    return res;
}

double norm_coomat(coomat* a) {
    return sqrt(scalar_coomat(copy_coomat(a), a));
}

coomat* coomat_remove_col(coomat* mat, size_t col) {
    if(col >= mat->size_j) {
        return mat;
    }
    coomat* res = init_coomat(mat->size_i, mat->size_j - 1);
    for (size_t i = 0; i < mat->size_i; i++) {
        for (size_t j = 0; j < col; j++) {
            coomat_set_value(res, i, j, coomat_read_value(mat, i, j));
        }
        for (size_t j = col; j + 1 < mat->size_i; j++) {
            coomat_set_value(res, i, j, coomat_read_value(mat, i, j + 1));
        }
    }
    
    free_coomat(mat);
    return res;
}

coomat* coomat_remove_line(coomat* mat, size_t line){
    if(line >= mat->size_i) {
        return mat;
    }
    coomat* res = init_coomat(mat->size_i - 1, mat->size_j);
    for (size_t i = 0; i < mat->size_i; i++) {
        for (size_t j = 0; j < line; j++) {
            coomat_set_value(res, j, i, coomat_read_value(mat, j, i));
        }
        for (size_t j = line; j + 1 < mat->size_i; j++) {
            coomat_set_value(res, j, i, coomat_read_value(mat, j + 1, i));
        }
    }
    free_coomat(mat);
    return res;
}

coomat* coomat_from_array(size_t size_i, size_t size_j, double* values) {
    coomat* res = init_coomat(size_i, size_j);

    for (size_t i = 0; i < size_i * size_j; i++){
        res->values[i] = values[i]; 
    }

    return res;
}
