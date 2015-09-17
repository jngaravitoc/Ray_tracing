#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "growth.h"
#include "pk.h"

#define DELTA_CRIT 1.686f
#define SQRT0707 0.840832920383116300f
#define MRS_1_SQRT2PI 0.398942280401432678f

/* Sheth-Tormen bias, expressed as M*dn/dM(M,z) */
float stbias(float m, float zz, GROWD growdat, PK *ppk)
{
  float dlin,density,scale;
  float w,s,nu;
  
  float bias;

  float z,corr;

  z=zz+0.0;

  corr =  dlina(1.f/(1.f+z),growdat)/ dlina(1.f/(1.f+zz),growdat);
  //fprintf(stderr," %f \n",corr);


  if (z>0) {
    dlin = dlina(1.f/(1.f+z),growdat) / dlin0(growdat);
    w = DELTA_CRIT / dlin;

    //fprintf(stderr,"dlin = %f \n",dlin);

  } else if (z==0) {
    w = DELTA_CRIT;
  } else {
    fprintf(stderr,"stmassfunc error: z<0 is not allowed; z=%g\n",z);
    exit(EXIT_FAILURE);
  }

  density = (ppk->eh).omhh * 2.7755e11f; // density at z=0 in Msun/Mpc^3

  scale = powf(3.f*m/(4.f*M_PI*density),1.f/3.f);

  s = sigma2(scale,ppk);
  
  nu = w*w/s;
  
   
  bias = 1.0;
  bias += corr * (SQRT0707*nu-1.0)/DELTA_CRIT; 
  bias += corr * 0.6/DELTA_CRIT/(1.0+powf(SQRT0707*nu,0.3f));
  

  //  bias = 1 + (nu-1.0)/DELTA_CRIT;   //Classic PS formula

  return bias;


}
