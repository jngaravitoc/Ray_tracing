#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "growth.h"

float mrsromb(float (*func)(float x, void *fdata), void *fdata, 
            float a, float b);

int growinit(float omegam, float omegal, float omegak, GROWD *pgrowdat)
{
  pgrowdat->omegam = omegam;
  pgrowdat->omegal = omegal;
  pgrowdat->omegak = omegak;

  return(EXIT_SUCCESS);
}

float dlin0(GROWD growdat)
{
  return mrsromb(inva3h3,(void *)&growdat,0,1);
}

float dlina(float a, GROWD growdat)
{
  float h;
  float integral;

  h = sqrt(1/(a*a*a)*growdat.omegam 
              + growdat.omegal 
              + 1/(a*a)*growdat.omegak);

  integral = mrsromb(inva3h3,(void *)&growdat,0,a);

  return (h * integral);
}

inline float inva3h3(float a, void *growdat)
{
  float ah;

  if (a>0) {

    ah = sqrt(1/a*((GROWD *)growdat)->omegam 
              + a*a*((GROWD *)growdat)->omegal 
              + ((GROWD *)growdat)->omegak);

    return 1/(ah*ah*ah);
  } else {
    return 0;
  }
}
