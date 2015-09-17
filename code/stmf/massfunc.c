#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "growth.h"
#include "pk.h"

#define DELTA_CRIT 1.686f
#define SQRT0707 0.840832920383116300f
#define MRS_1_SQRT2PI 0.398942280401432678f

/* Print out sigma */
float sigma_pp(float m, float z, GROWD growdat, PK *ppk)
{
  float dlin,density,scale;
  float w,s,nu,nut;
  float mdsdm,mdndm;

  if (z>0) {
    dlin = dlina(1.f/(1.f+z),growdat) / dlin0(growdat);
    w = DELTA_CRIT / dlin;
  } else if (z==0) {
    w = DELTA_CRIT;
  } else {
    fprintf(stderr,"stmassfunc error: z<0 is not allowed; z=%g\n",z);
    exit(EXIT_FAILURE);
  }

  density = (ppk->eh).omhh * 2.7755e11f; // density at z=0 in Msun/Mpc^3

  scale = powf(3.f*m/(4.f*M_PI*density),1.f/3.f);

  s = sigma2(scale,ppk);

  return sqrtf(s);

}

/* Sheth-Tormen halo mass function, expressed as M*dn/dM(M,z) */
float stmdndm(float m, float z, GROWD growdat, PK *ppk)
{
  float dlin,density,scale;
  float w,s,nu,nut;
  float mdsdm,mdndm;

  if (z>0) {
    dlin = dlina(1.f/(1.f+z),growdat) / dlin0(growdat);
    w = DELTA_CRIT / dlin;
  } else if (z==0) {
    w = DELTA_CRIT;
  } else {
    fprintf(stderr,"stmassfunc error: z<0 is not allowed; z=%g\n",z);
    exit(EXIT_FAILURE);
  }

  density = (ppk->eh).omhh * 2.7755e11f; // density at z=0 in Msun/Mpc^3

  scale = powf(3.f*m/(4.f*M_PI*density),1.f/3.f);

  s = sigma2(scale,ppk);

  nu = w/sqrtf(s);

  nut = SQRT0707 * nu; // sqrt(0.707) * nu

  mdsdm = mdsigma2dm(scale,ppk);

  mdndm = density / m / s * -mdsdm * MRS_1_SQRT2PI *
    0.3222f * (1.f + 1.f/powf(nut,0.6f)) * nut * expf(-nut*nut/2.f);

  return mdndm;
}


/* Press Schechter halo mass function, expressed as M*dn/dM(M,z) */
float psmdndm(float m, float z, GROWD growdat, PK *ppk)
{
  float dlin,density,scale;
  float w,s,nu;
  float mdsdm,mdndm;

  if (z>0) {
    dlin = dlina(1.f/(1.f+z),growdat) / dlin0(growdat);
    w = DELTA_CRIT / dlin;
  } else if (z==0) {
    w = DELTA_CRIT;
  } else {
    fprintf(stderr,"stmassfunc error: z<0 is not allowed; z=%g\n",z);
    exit(EXIT_FAILURE);
  }

  density = (ppk->eh).omhh * 2.7755e11f; // density at z=0 in Msun/Mpc^3

  scale = powf(3.f*m/(4.f*M_PI*density),1.f/3.f);

  s = sigma2(scale,ppk);

  nu = w/sqrtf(s);

  mdsdm = mdsigma2dm(scale,ppk);

  mdndm = density / m / s * -mdsdm * MRS_1_SQRT2PI *
    nu * expf(-nu*nu/2.f);

  return mdndm;
}
