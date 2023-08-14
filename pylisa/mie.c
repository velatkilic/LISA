/*
 * Modified from: https://omlc.org/software/mie/miesphr.c
 * 
 * We translated the Fortran code written by Tony Durkin at UT, Austin
 * into ANSI C. Sept. 14, 1993. See main() for authors' information.
 * Information from main:
 * 
 * SPHERE MIE SCATTERING PROGRAM
 * Lihong Wang, Ph.D.; Steven L. Jacques, Ph.D.
 * Laser Biology Research Laboratory
 * Univ. of Texas / M.D. Anderson Cancer Center
 * Houston, Texas, USA.
 * Version 1.3, 09/08/1994 - 08/09/1995
 * Acknowledgement: Tony Durkin; Craig Gardner. Univ. of Texas, Austin.
 * Reference:       Craig F. Bohren; Donald R. Huffman.
 *                  Absorption and Scattering of Light by Small Particles.
 *                  John Wiley & Sons, Inc. 1983.
 * 
 * To compile, type: cc -o miesphr miesphr.c or acc -o miesphr miesphr.c
 * whichever compiler complies to ANSI C.
 */
/**********************************************************************
Program written to calculate the scattering coefficient (mus) and
scattering anisotropy (g) for wavelength region specified by the user
the program uses subroutine callbh written by Bohren et al:
	absorption and scattering of light by small particles.
**********************************************************************/

#include "mie.h"

complex
Cform(float rl, float im)
{
  complex     c;

  c.r = rl;
  c.i = im;
  return (c);
}

float
Creal(complex C)
{
  return (C.r);
}

float
Cimag(complex C)
{
  return (C.i);
}

float
Cabs(complex C)
{
  return (sqrt(C.r * C.r + C.i * C.i));
}

complex
Cadd(complex C1, complex C2)
{
  complex     c;

  c.r = C1.r + C2.r;
  c.i = C1.i + C2.i;
  return (c);
}

complex
Csub(complex C1, complex C2)
{
  complex     c;

  c.r = C1.r - C2.r;
  c.i = C1.i - C2.i;
  return (c);
}

/* (a + ib)(c + id) = ac-bd + i(bc+ad) */
complex
Cmulti(complex C1, complex C2)
{
  complex     c;

  c.r = C1.r * C2.r - C1.i * C2.i;
  c.i = C1.r * C2.i + C1.i * C2.r;
  return (c);
}

/* (a + ib)/(c + id) = (a+ib)(c-id)/(c^2+d^2) */
complex
Cdiv(complex C1, complex C2)
{
  float       temp;
  complex     c;

  temp = 1 / (C2.r * C2.r + C2.i * C2.i);
  c.r = (C1.r * C2.r + C1.i * C2.i) * temp;
  c.i = (C1.i * C2.r - C1.r * C2.i) * temp;
  return (c);
}

complex
Cconj(complex C)
{
  complex     ctemp;

  ctemp.r = C.r;
  ctemp.i = -C.i;
  return (ctemp);
}

/**************************************************************************
subroutine BHMie calculates amplitude scattering matrix
elements and efficiencies for extinction, total scattering
and backscattering for a given size parameter and
relative refractive index
*************************************************************************/
void
BHMie(float *X,
      complex * RefRel,
      int *Nang,
      complex * S1,
      complex * S2,
      float *Qext,
      float *Qsca,
      float *Qback,
      float *Ganiso)
{
  // static complex cd[3000];
  static complex cy;
  static complex can, cbn;
  static complex cxi;
  static complex cxi0, cxi1;
  static complex can1, cbn1, can2, cbn2;
  complex     ctemp1, ctemp2, ctemp3;

  static float dang, apsi, ymod;
  static double qsca1;
  static float apsi0, apsi1;
  static int  j, n;
  static float p, t;
  static float theta[200];
  static int  nstop;
  static float xstop;
  static double dn;
  static float fn;
  static int  jj;
  static float pi[200];
  static float pi0[200], pi1[200];
  static double dx;
  static int  nn;
  static float rn;
  static float ganisotmp;
  static int  rn1;
  static float chi, mu[200], tau[200];
  static double psi;
  static int  nmx;
  static float chi0, chi1;
  static double psi0, psi1;

  dx = *X;
  cy = Cform((*X) * RefRel->r, (*X) * RefRel->i);

  /*************************************************************************
  series terminated after nstop terms
  *************************************************************************/
  xstop = *X + 4 * pow(*X, 0.3333) + 2.0;
  nstop = xstop;
  ymod = Cabs(cy);
  nmx = (xstop > ymod) ? xstop : ymod + 14;
  complex cd[nmx];
  dang = PI / (2 * (*Nang - 1));

  for (j = 0; j < *Nang; j++) {
    theta[j] = j * dang;
    mu[j] = cos(theta[j]);
  }

  /************************************************************************
  logarithmic derivative cd(j) calculated by downward
  recurrence beginning with initial value 0.0+i*0.0
  at j=nmx
  ***********************************************************************/
  cd[nmx] = Cform(0.0, 0.0);
  nn = nmx;

  for (n = 0; n < nn; n++) {
    rn = nmx - n + 1;
    ctemp1 = Cform(rn, 0.0);
    ctemp1 = Cdiv(ctemp1, cy);
    ctemp2 = Cadd(cd[nmx - n], ctemp1);
    ctemp3 = Cform(1, 0);
    ctemp3 = Cdiv(ctemp3, ctemp2);
    cd[nmx - n - 1] = Csub(ctemp1, ctemp3);
  }

  for (j = 0; j < *Nang; j++) {
    pi0[j] = 0.0;
    pi1[j] = 1.0;
  }

  nn = (*Nang << 1) - 1;
  for (j = 0; j < nn; j++)
    S1[j] = S2[j] = Cform(0.0, 0.0);

  /************************************************************************
  riccati-bessel functions with real argument X
  calculated by upward recurrence
  *********************************************************************/
  psi0 = cos(dx);
  psi1 = sin(dx);
  chi0 = -sin(*X);
  chi1 = cos(*X);
  apsi0 = psi0;
  apsi1 = psi1;
  cxi0 = Cform(apsi0, -chi0);
  cxi1 = Cform(apsi1, -chi1);
  *Qsca = 0.0;
  *Ganiso = 0.0;

  can1 = Cform(0.0, 0.0);
  cbn1 = Cform(0.0, 0.0);
  can2 = Cform(0.0, 0.0);
  cbn2 = Cform(0.0, 0.0);
  qsca1 = 0.0;

  n = 1;
  do {
    dn = (double) n;
    rn = (float) n;
    fn = (2. * rn + 1.) / (rn * (rn + 1.));
    psi = (2. * dn - 1.) * psi1 / dx - psi0;
    apsi = psi;
    chi = (2. * rn - 1.) * chi1 / (*X) - chi0;
    cxi = Cform(apsi, -chi);

    ctemp1 = Cdiv(cd[n - 1], *RefRel);
    ctemp1 = Cform(ctemp1.r + rn / (*X), ctemp1.i);
    ctemp2 = Cform(ctemp1.r * apsi - apsi1, ctemp1.i * apsi);
    ctemp3 = Cmulti(ctemp1, cxi);
    ctemp3 = Csub(ctemp3, cxi1);
    can = Cdiv(ctemp2, ctemp3);

    ctemp1 = Cmulti(cd[n - 1], *RefRel);
    ctemp1 = Cform(ctemp1.r + rn / (*X), ctemp1.i);
    ctemp2 = Cform(ctemp1.r * apsi - apsi1, ctemp1.i * apsi);
    ctemp3 = Cmulti(ctemp1, cxi);
    ctemp3 = Csub(ctemp3, cxi1);
    cbn = Cdiv(ctemp2, ctemp3);

    can1 = can;
    cbn1 = cbn;
    *Qsca += (2 * rn + 1) * (Cabs(can) * Cabs(can) + Cabs(cbn) * Cabs(cbn));

    for (j = 0; j < *Nang; j++) {
      jj = (*Nang << 1) - j - 2;
      pi[j] = pi1[j];
      tau[j] = rn * mu[j] * pi[j] - (rn + 1) * pi0[j];
      p = pow(-1, n - 1);

      ctemp1 = Cform(fn * (can.r * pi[j] + cbn.r * tau[j]),
		     fn * (can.i * pi[j] + cbn.i * tau[j]));
      S1[j] = Cadd(S1[j], ctemp1);
      t = pow(-1., n);

      ctemp1 = Cform(fn * (can.r * tau[j] + cbn.r * pi[j]),
		     fn * (can.i * tau[j] + cbn.i * pi[j]));
      S2[j] = Cadd(S2[j], ctemp1);

      if (j == jj)
	continue;

      ctemp1 = Cform(fn * (can.r * pi[j] * p + cbn.r * tau[j] * t),
		     fn * (can.i * pi[j] * p + cbn.i * tau[j] * t));
      S1[jj] = Cadd(S1[jj], ctemp1);

      ctemp1 = Cform(fn * (can.r * tau[j] * t + cbn.r * pi[j] * p),
		     fn * (can.i * tau[j] * t + cbn.i * pi[j] * p));
      S2[jj] = Cadd(S2[jj], ctemp1);
    }

    psi0 = psi1;
    psi1 = psi;
    apsi1 = psi1;
    chi0 = chi1;
    chi1 = chi;
    cxi1 = Cform(apsi1, -chi1);
    rn1 = rn;
    n = n + 1;
    rn = n;

    for (j = 0; j < *Nang; j++) {
      pi1[j] = ((2. * rn - 1.) / (rn - 1.)) * mu[j] * pi[j];
      pi1[j] = pi1[j] - rn * pi0[j] / (rn - 1.);
      pi0[j] = pi[j];
    }

    dn = n;
    rn = n;
    fn = (2. * rn + 1.) / (rn * (rn + 1.));
    psi = (2. * dn - 1.) * psi1 / dx - psi0;
    apsi = psi;
    chi = (2. * rn - 1.) * chi1 / (*X) - chi0;
    cxi = Cform(apsi, -chi);

    ctemp1 = Cdiv(cd[n - 1], *RefRel);
    ctemp1 = Cform(ctemp1.r + rn / (*X), ctemp1.i);
    ctemp2 = Cform(ctemp1.r * apsi - apsi1, ctemp1.i * apsi);
    ctemp3 = Cmulti(ctemp1, cxi);
    ctemp3 = Csub(ctemp3, cxi1);
    can = Cdiv(ctemp2, ctemp3);

    ctemp1 = Cmulti(cd[n - 1], *RefRel);
    ctemp1 = Cform(ctemp1.r + rn / (*X), ctemp1.i);
    ctemp2 = Cform(ctemp1.r * apsi - apsi1, ctemp1.i * apsi);
    ctemp3 = Cmulti(ctemp1, cxi);
    ctemp3 = Csub(ctemp3, cxi1);
    cbn = Cdiv(ctemp2, ctemp3);

    can2 = can;
    cbn2 = cbn;

    ctemp1 = Cmulti(can1, Cconj(can2));
    ctemp2 = Cmulti(cbn1, Cconj(cbn2));
    ctemp3 = Cadd(ctemp1, ctemp2);
    ganisotmp = rn1 * (rn1 + 2.0) * Creal(ctemp3);

    *Ganiso += ganisotmp / (rn1 + 1.);
    ctemp1 = Cmulti(can1, Cconj(cbn1));
    *Ganiso += (2. * rn1 + 1.) * Creal(ctemp1) / (rn1 * (rn1 + 1.0));
  } while (n - 1 - nstop < 0);

  *Qsca = (2. / (*X * *X)) * *Qsca;
  *Qext = (4. / (*X * *X)) * Creal(S1[0]);
  *Qback = (4. / (*X * *X)) *
    Cabs(S1[2 * (*Nang) - 2]) * Cabs(S1[2 * (*Nang) - 2]);
  *Ganiso = (4. / (*X * *X)) * *Ganiso;
  *Ganiso = *Ganiso / *Qsca;
}