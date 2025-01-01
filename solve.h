#ifndef CVIDE_SOLVE_H
#define CVIDE_SOLVE_H

#include "csrmat.h"
#include "coomat.h"

coomat_list* Euler_explicit(coomat* (*f)(coomat*), coomat* y0, double step, int step_count);
coomat_list* RK4_explicit(coomat* (*f)(coomat*), coomat* y0, double step, int step_count);
coomat_list* coo_RK4_linear_explicit(coomat* A, coomat* b, coomat* y0, double step, int step_count);
coomat_list* csr_RK4_linear_explicit(csrmat* A, coomat* b, coomat* y0, double step, int step_count);

#endif
