/* file SIR.c */
#include <R.h>
static double parms[2];

#define iota parms[0]
#define rho parms[1]

/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
  int N=2;
  odeparms(&N, parms);
}

/* Derivatives */
  void derivs (int *neq, double *t, double *y, double *ydot,
               double *yout, int *ip)
{
  ydot[0] = - iota * y[1] * y[0] / (y[0] + y[1] + y[2]);
  ydot[1] = iota * y[1] * y[0] / (y[0] + y[1] + y[2]) - rho * y[1];
  ydot[2] = rho * y[1];
  ydot[3] = iota * y[1] * y[0] / (y[0] + y[1] + y[2]);
}
/* END file SIR.c */