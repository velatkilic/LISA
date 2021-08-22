#ifndef UTILS_H
#define UTILS_h

#include "ParticleDist.h"
#include <vector>
#include <math.h>
#include <complex>

double trapz(std::vector<double> f, std::vector<double> x);
std::vector<double> logspace(double start, double end, int N);
double calc_qext(double x, std::complex<double> nref);

#endif // UTILS_H
