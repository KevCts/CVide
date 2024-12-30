#ifndef CVIDE_SOLVE_H
#define CVIDE_SOLVE_H

#include "csrmat.h"
#include "coomat.h"

coomat_list* Euler(coomat* (*f)(coomat*), coomat* y0, double step, int step_count);
coomat_list* RK4(coomat* (*f)(coomat*), coomat* y0, double step, int step_count);
coomat_list* RK4_linear(coomat* A, coomat* b, coomat* y0, double step, int step_count);

#endif
