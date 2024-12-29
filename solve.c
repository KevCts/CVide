#include "solve.h"
#include "coomat.h"

coomat_list* Euler(coomat* (*f)(coomat*), coomat* y0, double step, int step_count) {
    coomat_list* y = init_coomat_list(y0->size_i, y0->size_j);

    y = add_coomat_to_list(y, y0);

    for (int i = 0; i < step_count; i++) {
        y = add_coomat_to_list(y, sum_coomat(coomat_from_list(y, i), dot_coomat(step, f(coomat_from_list(y, i)))));
    }

    return y;
}

coomat_list* RK4(coomat* (*f)(coomat*), coomat* y0, double step, int step_count) {
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
