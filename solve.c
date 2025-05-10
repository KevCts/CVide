#include "solve.h"

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

coomat* minres(coomat* a, coomat* b, coomat* x0, double eps, int max_it) {
    coomat_list* x = init_coomat_list(x0->size_i, x0->size_j);
    coomat_list* r = init_coomat_list(x0->size_i, x0->size_j);
    coomat_list* s = init_coomat_list(x0->size_i, x0->size_j);
    coomat_list* p = init_coomat_list(x0->size_i, x0->size_j);
    double alpha, beta;
    int k = 0;

    add_coomat_to_list(x, x0);
    add_coomat_to_list(r, sum_coomat(b, dot_coomat(-1, prod_coomat(copy_coomat(a), coomat_from_list(x, 0)))));
    add_coomat_to_list(p, coomat_from_list(r, 0));
    add_coomat_to_list(s, prod_coomat(copy_coomat(a), coomat_from_list(p, 0)));

    while (k < max_it){
        k++;
        alpha = scalar_coomat(coomat_from_list(r, k-1), coomat_from_list(s, k-1))/scalar_coomat(coomat_from_list(s, k-1), coomat_from_list(s, k-1));
        add_coomat_to_list(x, sum_coomat(coomat_from_list(x, k-1), dot_coomat(alpha, coomat_from_list(p, k-1))));
        add_coomat_to_list(r, sum_coomat(coomat_from_list(r, k-1), dot_coomat(-1 * alpha, coomat_from_list(s, k-1))));

        if (norm_coomat(coomat_from_list(r, k)) < eps)
            break;

        coomat* p_temp = coomat_from_list(s, k-1);
        coomat* s_temp = prod_coomat(copy_coomat(a), coomat_from_list(s, k-1));


        beta = scalar_coomat(copy_coomat(s_temp), coomat_from_list(s, k-1))/scalar_coomat(coomat_from_list(s, k-1), coomat_from_list(s, k-1));
        p_temp = sum_coomat(p_temp, dot_coomat(-1 * beta, coomat_from_list(p, k-1)));
        s_temp = sum_coomat(s_temp, dot_coomat(-1 * beta, coomat_from_list(s, k-1)));

        if (k > 1) {
            beta = scalar_coomat(coomat_from_list(s, k), coomat_from_list(s, k-2))/scalar_coomat(coomat_from_list(s, k-2), coomat_from_list(s, k-2));
            p_temp = sum_coomat(p_temp, dot_coomat(-1 * beta, coomat_from_list(p, k-2)));
            s_temp = sum_coomat(s_temp, dot_coomat(-1 * beta, coomat_from_list(s, k-2)));
        }

        add_coomat_to_list(p, p_temp);
        add_coomat_to_list(s, s_temp);
    }

    coomat* res = coomat_from_list(x, k);

    free_coomat_list(x);
    free_coomat_list(r);
    free_coomat_list(s);
    free_coomat_list(p);

    free_coomat(a);

    return res;
}
