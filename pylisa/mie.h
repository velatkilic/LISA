#ifndef MIE_H
#define MIE_H

#include <math.h>
#include <stdio.h>
#include <float.h>

#define PI 		3.14159265

typedef struct {
  float       r, i;
}           complex;
complex Cform(float rl, float im);
float Creal(complex C);
float Cimag(complex C);
float Cabs(complex C);
complex Cadd(complex C1, complex C2);
complex Csub(complex C1, complex C2);
complex Cmulti(complex C1, complex C2);
complex Cdiv(complex C1, complex C2);
complex Cconj(complex C);

void BHMie(float *X, complex * RefRel, int *Nang, complex * S1, complex * S2, float *Qext, float *Qsca, float *Qback, float *Ganiso);

#endif // MIE_h