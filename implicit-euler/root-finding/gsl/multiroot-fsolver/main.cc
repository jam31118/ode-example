#include <cstdio>
#include <cstdlib>

#include "gsl/gsl_multiroots.h"


struct f_params {
  double a;
  double b;
};

int rosenbrock_f (const gsl_vector *x, void *params, gsl_vector *f) 
{
  const double a = ((struct f_params *) params)->a;
  const double b = ((struct f_params *) params)->b;

  const double x0 = gsl_vector_get(x, 0);
  const double x1 = gsl_vector_get(x, 1);

  const double f0 = a * (1.0 - x0);
  const double f1 = b * (x1 - x0 * x0);

  gsl_vector_set(f, 0, f0);
  gsl_vector_set(f, 1, f1);

  return GSL_SUCCESS;
}

void print_state(size_t i_iter, gsl_multiroot_fsolver *s) {
  printf("[ i_iter = %03lu ] x = % .3f % .3f / f(x) = % .3e % .3e\n",
      i_iter, 
      gsl_vector_get(s->x,0), gsl_vector_get(s->x,1),
      gsl_vector_get(s->f,0), gsl_vector_get(s->f,1));
}


int main(int argc, char *argv[]) {
  
  //// Set up solver
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  const size_t N_dim = 2;
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc(T, N_dim);
  struct f_params p = {1.0, 10.0};
  gsl_multiroot_function f = {&rosenbrock_f, N_dim, &p};
  const double x_init[2] = {-10.0, -5.0};
  gsl_vector *x = gsl_vector_alloc(N_dim);
  gsl_vector_set(x, 0, x_init[0]);
  gsl_vector_set(x, 1, x_init[1]);
  gsl_multiroot_fsolver_set(s, &f, x);
  
  
  //// Iterate
  size_t i_iter = 0;
  print_state(i_iter, s); 


  //// Deallocate
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);

  return EXIT_SUCCESS;
}
