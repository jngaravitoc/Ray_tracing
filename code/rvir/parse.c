#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "rvir.h"

int usage(int argc, char *argv[])
{
  fprintf(stderr,"Usage: %s [mass] [z] [-om Omega_m] [-v verbose]\n",*argv);
  fprintf(stderr,"Type %s --help for parameter details\n",*argv);
  exit(EXIT_FAILURE);
}

int longusage(int argc, char *argv[])
{
  fprintf(stderr,"Usage: %s [mass] [z] [-om Omega_m] [-v verbose]\n",*argv);
  fprintf(stderr,"Parameters (default value):\n");
  fprintf(stderr,"   mass    - halo mass [1e12 h^-1 Msun]\n");
  fprintf(stderr,"   z       - redshift of halo (0)\n");
  fprintf(stderr,"   om      - Omega_matter (0.3)\n");
  fprintf(stderr,"   v       - verbosity (1)\n");
  fprintf(stderr,"R_vir is output in units of h^-1 kpc\n");
  exit(EXIT_SUCCESS);
}

int parseinput(int argc, char *argv[], RVIR_IN *pinput)
{
  float mass = 1.e12; // halo mass
  float z = 0; // redshift to evaluate mass function
  float omegam = 0.3; // Omegas
  int iarg=1,iparam=1,verbose=1; // counter
  RVIR_IN in;

  // parse input
  while (iarg < argc) {

    switch (**(argv+iarg)) { // check to see if first character is a flag

    case '-': // the first character is a flag

      if (!strcmp(argv[iarg],"-om")) {
        ++iarg; if (iarg >= argc) usage(argc,argv);
        omegam = atof(argv[iarg]); ++iarg;
      }
      else if (!strcmp(argv[iarg],"-v")) {
        ++iarg; if (iarg >= argc) usage(argc,argv);
        verbose = atof(argv[iarg]); ++iarg;
      }
      else if (!strcmp(argv[iarg],"--help")) {
        longusage(argc,argv);
      }
      else usage(argc,argv);
      break;

    default: // the first character is not a flag

      switch(iparam++) {

      case 1:
        mass = atof(argv[iarg]); ++iarg;
	break;

      case 2:
        z = atof(argv[iarg]); ++iarg;
	break;

      default:
	usage(argc,argv);
	break;

      }
      break;

    } // end first character switch

  } // end while (iarg < argc)

  in.mass = mass;
  in.z = z;
  in.omegam = omegam;
  in.verbose = verbose;

  *pinput = in; // copy input data

  // echo command
  if (verbose>=1) {
    for (iarg=0;iarg<argc;iarg++)
      printf("%s ",argv[iarg]);
    printf("\n");
  }

  return(EXIT_SUCCESS);
}
