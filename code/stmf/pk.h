#ifndef PK_HINCLUDED
#define PK_HINCLUDED

#ifndef EH_HINCLUDED
#include "eishu.h"
#endif

#ifndef GR_HINCLUDED
#include "growth.h"
#endif

typedef struct pk_pass {
  int setnorm; // flag to indicate if normsig8 has been run
  float scale; // length scale to evaluate sigma2 at, in Mpc
  float anorm; // A in P(k) = A k^ns T(k)^2
  float ns; // primordial spectral index
  float sig8; // desired value of sigma_8
  EH eh; // structure containing transfer function variables
} PK;

float dsigtopdlnk(float lnk, void *ppk);
float mrsromb(float (*func)(float x, void *fdata), void *fdata, 
            float a, float b);
float sigma2(float scale, PK *ppk);
int normsig8(float sigma8, float h, float ns, PK *ppk);
int pkinit(float omegam, float omegab, float h, 
           float sigma8, float ns, PK *ppk);
float powerk(float k, PK *ppk);
float mdsigma2dm(float scale, PK *ppk);
float ddsigtopdmdlnk(float lnk, void *ppk);
float stmdndm(float m, float z, GROWD growdat, PK *ppk);
float stbias(float m, float z, GROWD growdat, PK *ppk);
/* Press Schechter halo mass function, expressed as M*dn/dM(M,z) */
float psmdndm(float m, float z, GROWD growdat, PK *ppk);
float sigma_pp(float m, float z, GROWD growdat, PK *ppk);


#endif
