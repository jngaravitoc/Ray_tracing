#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "growth.h"
#include "pk.h"
#include <assert.h>

int usage(int argc, char *argv[])
{
  fprintf(stderr,"Usage: %s [z [lgm1 [lgm2 [nm]]]]\n",*argv);
  fprintf(stderr,"            [-om Omega_m] [-ol Omega_l] [-ob Omega_b]\n");
  fprintf(stderr,"            [-h hubble] [-sig8 sigma_8] [-ns spectral index]\n");
  fprintf(stderr,"            [-v verbose] [-o outfile]\n");
  fprintf(stderr,"Type %s --help for default values\n",*argv);
  exit(EXIT_FAILURE);
}

int longusage(int argc, char *argv[])
{
  fprintf(stderr,"Usage: %s [z [lgm1 [lgm2 [nm]]]]\n",*argv);
  fprintf(stderr,"            [-om Omega_m] [-ol Omega_l] [-ob Omega_b]\n");
  fprintf(stderr,"            [-h hubble] [-sig8 sigma_8] [-ns spectral index]\n");
  fprintf(stderr,"            [-v verbose] [-o outfile]\n");
  fprintf(stderr,"Default values for the parameters:\n");
  fprintf(stderr,"            z = 0\n");
  fprintf(stderr,"            lgm1 = 6\n");
  fprintf(stderr,"            lgm2 = 15\n");
  fprintf(stderr,"            nm = 1000\n");
  fprintf(stderr,"            om = 0.28  ol = 0.72  ob = 0.0462\n");
  fprintf(stderr,"            h = 0.7  sig8 = 0.817  ns = 0.96\n");
  fprintf(stderr,"            v = 1  outfile = stmf.txt\n");
  exit(EXIT_FAILURE);
}

/* Front end for Sheth-Tormen mass function */
/* returns a file with mass and m*dndm */
int main(int argc, char *argv[])
{
  float z = 0; // redshift to evaluate mass function
  int nm = 1000; // number of masses to sample mass function
  float lgm1 = 6, lgm2 = 15; // mass function range in log_10 mass
  float omegam = 0.28, omegal = 0.72, omegab = 0.0462; // Omegas
  float h = 0.7, sigma8 = 0.817, ns = 0.96; // hubble, sigma_8, P(k) init slope
  PK pk = {0}; // required to initialize setnorm and settf
  GROWD growdat;
  int im;
  float lgm,m,mdndm, bs;
  int iarg=1,iparam=1,verbose=1; // counter
  char outfile[128];
  FILE *fd;

  //keep track of number comulative density of halos and average bias
  double lgm_old;
  double avg_bs;
  double cum_nh;
  double delta_nh;
  double mdndm_old;

  FILE *fdLF;
  fdLF=fopen("stmf_LF.dat","w");
  assert(fdLF);
  fprintf(fdLF,"857 \n");

  sprintf(outfile,"stmf.txt");

  // parse input
  while (iarg < argc) {

    switch (**(argv+iarg)) { // check to see if first character is a flag

    case '-': // the first character is a flag

      if (!strcmp(argv[iarg],"-om")) {
        ++iarg; if (iarg >= argc) usage(argc,argv);
        omegam = atof(argv[iarg]); ++iarg;
      }
      else if (!strcmp(argv[iarg],"-ol")) {
        ++iarg; if (iarg >= argc) usage(argc,argv);
        omegal = atof(argv[iarg]); ++iarg;
      }
      else if (!strcmp(argv[iarg],"-ob")) {
        ++iarg; if (iarg >= argc) usage(argc,argv);
        omegab = atof(argv[iarg]); ++iarg;
      }
      else if (!strcmp(argv[iarg],"-h")) {
        ++iarg; if (iarg >= argc) usage(argc,argv);
        h = atof(argv[iarg]); ++iarg;
      }
      else if (!strcmp(argv[iarg],"-sig8")) {
        ++iarg; if (iarg >= argc) usage(argc,argv);
        sigma8 = atof(argv[iarg]); ++iarg;
      }
      else if (!strcmp(argv[iarg],"-ns")) {
        ++iarg; if (iarg >= argc) usage(argc,argv);
        ns = atof(argv[iarg]); ++iarg;
      }
      else if (!strcmp(argv[iarg],"-v")) {
        ++iarg; if (iarg >= argc) usage(argc,argv);
        verbose = atof(argv[iarg]); ++iarg;
      }
      else if (!strcmp(argv[iarg],"-o")) {
        ++iarg; if (iarg >= argc) usage(argc,argv);
        strcpy(outfile,argv[iarg]); ++iarg;
      }
      else if (!strcmp(argv[iarg],"--help")) {
        longusage(argc,argv);
      }
      else usage(argc,argv);
      break;

    default: // the first character is not a flag

      switch(iparam++) {

      case 1:
        z = atof(argv[iarg]); ++iarg;
	break;

      case 2:
        lgm1 = atof(argv[iarg]); ++iarg;
	break;

      case 3:
        lgm2 = atof(argv[iarg]); ++iarg;
	break;

      case 4:
        nm = atoi(argv[iarg]); ++iarg;
	break;

      default:
	usage(argc,argv);
	break;

      }
      break;

    } // end first character switch

  } // end while (iarg < argc)

  // open output file
  fd = fopen(outfile,"w");
  if (!fd) {
    fprintf(stderr,"Failed to open outfile %s\n",outfile);
    exit(EXIT_FAILURE);
  }

  // echo input parameters to file and possibly terminal
  fprintf(fd,"# z=%g  lgm1=%g  lgm2=%g  nm=%d\n",z,lgm1,lgm2,nm);
  fprintf(fd,"# omega_m=%g  omega_l=%g  omega_b=%g\n",omegam,omegal,omegab);
  fprintf(fd,"# h=%g  sigma_8=%g  ns=%g\n",h,sigma8,ns);
  fprintf(fd,"%10g  %10g  %10g  %10d\n",z,lgm1,lgm2,nm);
  fprintf(fd,"%10g  %10g  %10g\n",omegam,omegal,omegab);
  fprintf(fd,"%10g  %10g  %10g\n",h,sigma8,ns);
  if (verbose>=1) {
    printf("z=%g  lgm1=%g  lgm2=%g  nm=%d\n",z,lgm1,lgm2,nm);
    printf("omega_m=%g  omega_l=%g  omega_b=%g\n",omegam,omegal,omegab);
    printf("h=%g  sigma_8=%g  ns=%g\n",h,sigma8,ns);
    printf("verbose=%d  outfile=%s\n",verbose,outfile);
  }

  // initialize transfer function and power spectrum (required to call sigma2)
  pkinit(omegam, omegab, h, sigma8, ns, &pk);

  // set growth function data
  growinit(omegam,omegal,1.f-omegam-omegal,&growdat);

  // average quantities
  avg_bs = 0;
  cum_nh = 0;
  delta_nh = 0;
  lgm_old=lgm2;
  mdndm_old = 0.0;

  double ccc=0;

  for (im=nm-1;im>=0;im--) {

    lgm = im*(lgm2-lgm1)/(nm-1)+lgm1;

    m = powf(10.f,lgm);

    mdndm = stmdndm(m,z,growdat,&pk);
    //mdndm = psmdndm(m,z,growdat,&pk);

    bs = stbias(m,z,growdat,&pk);

    double sigma_print;
    sigma_print = sigma_pp(m,z,growdat,&pk);

    //average quantities
    ccc+=1.0;
    delta_nh= (lgm_old-lgm)*(mdndm+mdndm_old)*log(10.0)/2.0;
    cum_nh+= delta_nh;
    avg_bs+=delta_nh*bs;
    lgm_old = lgm;
    mdndm_old = mdndm;
    if(im<(nm-1)) //discard first point since zero cumulative number density is tricky
      {
	fprintf(fdLF,"%e %e %e \n",m,mdndm,cum_nh);
      }

    fprintf(fd,"%10d %10f %10e %10e %6e  %10e %10e  %10e\n",im,lgm,m,mdndm,bs,cum_nh ,avg_bs/cum_nh,sigma_print);
    //fprintf(fd,"%10d %10f %10e %10e %6e  %10e %10e\n",im,lgm,m,mdndm,bs,cum_nh ,avg_bs/cum_nh);
    if (verbose>=2) 
      printf("im=%d  lgm=%g  m*dndm(%g Msun) = %g bias=%g \n",im,lgm,m,mdndm,bs);


  }

  fclose(fd); // close output file
  fclose(fdLF);

  return(EXIT_SUCCESS);
}
