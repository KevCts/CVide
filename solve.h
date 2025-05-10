#ifndef CVIDE_SOLVE_H
#define CVIDE_SOLVE_H

#include "csrmat.h"

coomat_list* Euler_explicit(coomat* (*f)(coomat*), coomat* y0, double step, int step_count);
coomat_list* RK4_explicit(coomat* (*f)(coomat*), coomat* y0, double step, int step_count);
coomat* minres(coomat*, coomat*, coomat*, double, int);

#endif
