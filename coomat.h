#ifndef CVIDE_COOMAT_H
#define CVIDE_COOMAT_H

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

typedef struct {
    size_t size_i, size_j;
    double* values;
} coomat;

typedef struct {
    size_t size_i, size_j, count;
    coomat** elements;
} coomat_list;

coomat* init_coomat(size_t, size_t);

coomat_list* init_coomat_list(size_t, size_t);

coomat_list* add_coomat_to_list(coomat_list*, coomat*);

coomat* coomat_from_list(coomat_list*, size_t);

void coomat_fun_to_list(coomat* (*)(coomat*), coomat_list*);

void double_fun_to_list(double (*)(double), coomat_list*);

void free_coomat_list(coomat_list*);

void free_coomat(coomat*);

void print_coomat(coomat*);

double coomat_read_value(coomat*, size_t, size_t);

double coomat_set_value(coomat*, size_t, size_t, double);

coomat* copy_coomat(coomat*);

coomat* sum_coomat(coomat*, coomat*);

coomat* fun_coomat(double (*)(double), coomat*);

coomat* dot_coomat(double, coomat*);

bool equal_coomat(coomat*, coomat*);

double scalar_coomat(coomat*, coomat*);

coomat* transpose_coomat(coomat*);

coomat* prod_coomat(coomat*, coomat*);

double norm_coomat(coomat*);

#endif
