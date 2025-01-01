#include "solve.h"
#include "coomat.h"
#include "csrmat.h"

coomat_list* Euler_explicit(coomat* (*f)(coomat*), coomat* y0, double step, int step_count) {
    coomat_list* y = init_coomat_list(y0->size_i, y0->size_j);

    y = add_coomat_to_list(y, y0);

    for (int i = 0; i < step_count; i++) {
        y = add_coomat_to_list(y, sum_coomat(coomat_from_list(y, i), dot_coomat(step, f(coomat_from_list(y, i)))));
    }

    return y;
}

coomat_list* RK4_explicit(coomat* (*f)(coomat*), coomat* y0, double step, int step_count) {
    coomat_list* y = init_coomat_list(y0->size_i, y0->size_j);

    y = add_coomat_to_list(y, y0);

    for (int i = 0; i < step_count; i++) {
        coomat* k1 = f(coomat_from_list(y, i));
        coomat* k2 = f(sum_coomat(coomat_from_list(y, i), dot_coomat(step/2, copy_coomat(k1))));
        coomat* k3 = f(sum_coomat(coomat_from_list(y, i), dot_coomat(step/2, copy_coomat(k2))));
        coomat* k4 = f(sum_coomat(coomat_from_list(y, i), dot_coomat(step, copy_coomat(k3))));
        y = add_coomat_to_list(y, sum_coomat(coomat_from_list(y, i), dot_coomat(step/6, sum_coomat(sum_coomat(k1, dot_coomat(2, k2)), sum_coomat(dot_coomat(2, k3), k4)))));
    }

    return y;
}

coomat_list* csr_RK4_linear_explicit(csrmat* A, coomat* b, coomat* y0, double step, int step_count) {
    coomat_list* y = init_coomat_list(y0->size_i, y0->size_j);

    y = add_coomat_to_list(y, y0);

    for (int i = 0; i < step_count; i++) {
        coomat* k1 = sum_coomat(copy_coomat(b), prod_csr_coo(copy_csrmat(A), coomat_from_list(y, i)));
        coomat* k2 = sum_coomat(copy_coomat(b), prod_csr_coo(copy_csrmat(A), sum_coomat(coomat_from_list(y, i), dot_coomat(step/2, copy_coomat(k1)))));
        coomat* k3 = sum_coomat(copy_coomat(b), prod_csr_coo(copy_csrmat(A), sum_coomat(coomat_from_list(y, i), dot_coomat(step/2, copy_coomat(k2)))));
        coomat* k4 = sum_coomat(copy_coomat(b), prod_csr_coo(copy_csrmat(A), sum_coomat(coomat_from_list(y, i), dot_coomat(step, copy_coomat(k3)))));
        y = add_coomat_to_list(y, sum_coomat(coomat_from_list(y, i), dot_coomat(step/6, sum_coomat(sum_coomat(k1, dot_coomat(2, k2)), sum_coomat(dot_coomat(2, k3), k4)))));
    }

    free_csrmat(A);
    free_coomat(b);

    return y;
}

coomat_list* coo_RK4_linear_explicit(coomat* A, coomat* b, coomat* y0, double step, int step_count) {
    coomat_list* y = init_coomat_list(y0->size_i, y0->size_j);

    y = add_coomat_to_list(y, y0);

    for (int i = 0; i < step_count; i++) {
        coomat* k1 = sum_coomat(copy_coomat(b), prod_coomat(copy_coomat(A), coomat_from_list(y, i)));
        coomat* k2 = sum_coomat(copy_coomat(b), prod_coomat(copy_coomat(A), sum_coomat(coomat_from_list(y, i), dot_coomat(step/2, copy_coomat(k1)))));
        coomat* k3 = sum_coomat(copy_coomat(b), prod_coomat(copy_coomat(A), sum_coomat(coomat_from_list(y, i), dot_coomat(step/2, copy_coomat(k2)))));
        coomat* k4 = sum_coomat(copy_coomat(b), prod_coomat(copy_coomat(A), sum_coomat(coomat_from_list(y, i), dot_coomat(step, copy_coomat(k3)))));
        y = add_coomat_to_list(y, sum_coomat(coomat_from_list(y, i), dot_coomat(step/6, sum_coomat(sum_coomat(k1, dot_coomat(2, k2)), sum_coomat(dot_coomat(2, k3), k4)))));
    }

    free_coomat(A);
    free_coomat(b);

    return y;
}
