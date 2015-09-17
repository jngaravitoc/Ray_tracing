#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


// returns R_vir(M)

#define G 6.6742e-8
#define MSUN 1.998e33
#define KPC 3.0857e21

int main(int argc, char **argv)
{
  int n_points = 1;
  float mass = atof(argv[1]);
  float z = atof(argv[2]);
  //double mass;
  //double z;
  double omegal;
  double oneplusz3;
  double omegamz; 
  double d;
  double Delta_c;
  double cosmofact;
  double rvir;
  double vc;
  int i;
  double PI=M_PI; // from math.h
  FILE *in;
  char filename[100] = "rvir_in.dat";
  double m, Z;
  float omegam = 0.3;
  /*
  mass = malloc(sizeof(double)*n_points);
  z = malloc(n_points*sizeof(double));
  oneplusz3 = malloc(sizeof(double)*n_points);
  omegamz = malloc(sizeof(double)*n_points);
  d = malloc(sizeof(double)*n_points);
  Delta_c = malloc(sizeof(double)*n_points);
  cosmofact = malloc(sizeof(double)*n_points);
  rvir = malloc(n_points*sizeof(double));
  vc = malloc(sizeof(double)*n_points);
 */
  if(!argc){

  printf("You have to specify the number of halos \n");
  }
  
/*
  in = fopen(filename, "r");
  if(!in){
  printf("Problem opening file %s\n",filename);
  exit(1);
  }

  for(i=0;i<n_points;i++){
  fscanf(in, "%lf %lf \n",&m, &Z);
  mass[i]=m;
  z[i] = Z;
  }
  fclose(in);
*/

  omegal = 1. - omegam;

  //for(i=0;i<n_points;i++){
    

  oneplusz3 = (1+z)*(1+z)*(1+z);

  omegamz = omegam*oneplusz3 / (omegam*oneplusz3 + omegal); // flat
  d = omegamz - 1;
  Delta_c = 18.*PI*PI + 82.*d - 39.*d*d; // Bryan & Norman formula

  cosmofact = pow( (omegam/omegamz) * (Delta_c/(18.*PI*PI)) , -0.33333 );
  rvir = 0.784 * pow(mass/1e8,0.33333) * cosmofact * (10./(1.+z));

  vc = sqrt(G*(mass*MSUN)/(rvir*KPC))/1.e5; // in km/s
  
  printf("%g \t %g \n", rvir ,vc);

  
  //}
  //printf("V_c(%g h^-1 M_sun, z=%f) = %g km/s\n",mass,z,vc);

  return(EXIT_SUCCESS);
}
