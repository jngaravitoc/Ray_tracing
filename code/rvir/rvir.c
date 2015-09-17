#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "rvir.h"

// returns R_vir(M)

#define G 6.6742e-8
#define MSUN 1.998e33
#define KPC 3.0857e21

int main(int argc, char *argv[])
{
  RVIR_IN in;
  float mass,z,omegam,omegal,oneplusz3,omegamz,d,Delta_c,cosmofact,rvir,vc;
  int verbose;
  double PI=M_PI; // from math.h

  // parse input
  parseinput(argc, argv, &in);

  mass = in.mass;
  z = in.z;
  omegam = in.omegam;
  verbose = in.verbose;
  omegal = 1. - in.omegam;
  oneplusz3 = (1+z)*(1+z)*(1+z);

  omegamz = omegam*oneplusz3 / (omegam*oneplusz3 + omegal); // flat
  d = omegamz - 1;
  Delta_c = 18.*PI*PI + 82.*d - 39.*d*d; // Bryan & Norman formula

  cosmofact = pow( (omegam/omegamz) * (Delta_c/(18.*PI*PI)) , -0.33333 );

  if (verbose>=2) {
    printf("mass: %g z: %f omegam: %f omegal: %f\n",
           mass,z,omegam,omegal);
    printf("omegam_z: %f Delta_c: %g cosmofact: %g\n",
           omegamz,Delta_c,cosmofact);
  }

  rvir = 0.784 * pow(mass/1e8,0.33333) * cosmofact * (10./(1.+z));

  printf("R_vir(%g h^-1 M_sun, z=%f) = %g h^-1 kpc (comoving) \t %g \n",mass,z,rvir, rvir);

  vc = sqrt(G*(mass*MSUN)/(rvir*KPC))/1.e5; // in km/s

  //printf("V_c(%g h^-1 M_sun, z=%f) = %g km/s\n",mass,z,vc);

  return(EXIT_SUCCESS);
}
