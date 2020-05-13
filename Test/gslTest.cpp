#include <stdio.h>
#include <iostream>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <vector>

double GaussPulse(double t, void* params){
    double *Params = (double *) params;
    double amp = Params[0], sigma = Params[1], freq = Params[2];
    return amp*std::exp(-t*t/sigma/sigma)*std::cos(freq*t);
}
int
main (void)
{
  double params[] = {1.0,1.0,1.0};
  gsl_function F;
  F.function = &GaussPulse;
  F.params = &params;

  double epsabs = 1e-8, epsrel = 1e-7;
  double result, abserr;
  size_t neval;
  gsl_integration_qng(&F, -2.0, 2.0, epsabs, epsrel, &result, &abserr, &neval);
  std::cout<<"result:"<<result<<"\n"<<"abserr:"<<abserr<<"\n"<<"neval:"<<neval<<"\n";
  return 0;
}