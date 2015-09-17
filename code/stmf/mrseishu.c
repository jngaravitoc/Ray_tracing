/* mrseishu
 *
 * Eisenstein & Hu (1997) transfer function code
 *
 * Two changes I made to the distributed E&H code were to strip out
 * code I don't use, and move the global variables into a structure
 * that is explicitly passed.  This may have a deleterious impact on
 * performance, but makes the routine more obviously thread safe
 * (assuming the user takes appropriate care of the EH structure).
 *
 * Later, I optimized for speed, changing all of the math.h function
 * calls to their floating point equivalents, and making floating
 * point constants explicit to cut out type conversions.  I also added
 * a new 'global' variables mrs_1_alpha_c, which is the reciprocal of
 * alpha_c (to save that being calculated); I also created sone
 * temporary variable (mrstmp*) for optimizations.
 *
 * First call mrstfset to set the relevant cosmological parameters in
 * the EH structure.  Then call mrstf to return T(k).  Note that I
 * follow E&H units convention below: always Mpc, never h^-1 Mpc.
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "eishu.h"

int mrstfset(float omegam, float omegab, float h, EH *peh)
{
  float omega0hh, f_baryon, Tcmb;

  if (peh->settf) {
    fprintf(stderr,"TF parameters already set (or bad EH initialization)\n");
    exit(EXIT_FAILURE);
  }

  omega0hh = omegam * h * h;
  f_baryon = omegab / omegam;
  Tcmb = 2.728; // hardwired here

  TFset_parameters(omega0hh, f_baryon, Tcmb, peh);

  peh->settf = 1; // indicates that TFset_parameters has been run

  return(EXIT_SUCCESS);
}

float mrstf(float k, EH *peh)
{
  if (!peh->settf) {
    fprintf(stderr,"Must set TF parameters before evaluating T(k)\n");
    exit(EXIT_FAILURE);
  }

  return TFfit_onek(k, NULL, NULL, *peh);
}

void TFset_parameters(float omega0hh, float f_baryon, float Tcmb, EH *peh)
     /* Set all the scalars quantities for Eisenstein & Hu 1997
        fitting formula */
     /* Input: omega0hh -- The density of CDM and baryons, in units of
        critical dens, multiplied by the square of the Hubble
        constant, in units of 100 km/s/Mpc */
     /*       f_baryon -- The baryon fraction of the total matter density */
     /*       Tcmb -- The temperature of the CMB in Kelvin.  Tcmb<=0 forces use
              of the COBE value of  2.728 K. */
     /* Output: Nothing, but set many global variables used in TFfit_onek(). 
        You can access them yourself, if you want. */
     /* Note: Units are always Mpc, never h^-1 Mpc. */
{
  // EH global vars 
  float omhh,obhh,theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality;
  float sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,beta_node,k_peak;
  float sound_horizon_fit,alpha_gamma,mrs_1_alpha_c;
  // local vars
  float z_drag_b1, z_drag_b2;
  float alpha_c_a1, alpha_c_a2, beta_c_b1, beta_c_b2, alpha_b_G, y;

  if (f_baryon<=0.0 || omega0hh<=0.0) {
    fprintf(stderr, "TFset_parameters(): Illegal input.\n");
    exit(1);
  }
  omhh = omega0hh;
  obhh = omhh*f_baryon;
  if (Tcmb<=0.0) Tcmb=2.728;	/* COBE FIRAS */
  theta_cmb = Tcmb/2.7;

  z_equality = 2.50e4*omhh/POW4(theta_cmb);  /* Really 1+z */
  k_equality = 0.0746*omhh/SQR(theta_cmb);

  z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
  z_drag_b2 = 0.238*pow(omhh,0.223);
  z_drag = 1291*pow(omhh,0.251)/(1+0.659*pow(omhh,0.828))*
    (1+z_drag_b1*pow(obhh,z_drag_b2));
    
  R_drag = 31.5*obhh/POW4(theta_cmb)*(1000/(1+z_drag));
  R_equality = 31.5*obhh/POW4(theta_cmb)*(1000/z_equality);

  sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*
    log((sqrt(1+R_drag)+sqrt(R_drag+R_equality))/(1+sqrt(R_equality)));

  k_silk = 1.6*pow(obhh,0.52)*pow(omhh,0.73)*(1+pow(10.4*omhh,-0.95));

  alpha_c_a1 = pow(46.9*omhh,0.670)*(1+pow(32.1*omhh,-0.532));
  alpha_c_a2 = pow(12.0*omhh,0.424)*(1+pow(45.0*omhh,-0.582));
  alpha_c = pow(alpha_c_a1,-f_baryon)*
    pow(alpha_c_a2,-CUBE(f_baryon));
  mrs_1_alpha_c = 1.f/alpha_c; // introduced for optimization

  beta_c_b1 = 0.944/(1+pow(458*omhh,-0.708));
  beta_c_b2 = pow(0.395*omhh, -0.0266);
  beta_c = 1.0/(1+beta_c_b1*(pow(1-f_baryon, beta_c_b2)-1));

  y = z_equality/(1+z_drag);
  alpha_b_G = y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
  alpha_b = 2.07*k_equality*sound_horizon*pow(1+R_drag,-0.75)*alpha_b_G;

  beta_node = 8.41*pow(omhh, 0.435);
  beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*sqrt(pow(17.2*omhh,2.0)+1);

  k_peak = 2.5*3.14159*(1+0.217*omhh)/sound_horizon;
  sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1+10.0*pow(obhh,0.75));

  alpha_gamma = 1-0.328*log(431.0*omhh)*f_baryon + 0.38*log(22.3*omhh)*
    SQR(f_baryon);

  // set values in structure
  (*peh).omhh = omhh;
  (*peh).obhh = obhh;
  (*peh).theta_cmb = theta_cmb;
  (*peh).z_equality = z_equality;
  (*peh).k_equality = k_equality;
  (*peh).z_drag = z_drag;
  (*peh).R_drag = R_drag;
  (*peh).R_equality = R_equality;
  (*peh).sound_horizon = sound_horizon;
  (*peh).k_silk = k_silk;
  (*peh).alpha_c = alpha_c;
  (*peh).mrs_1_alpha_c = mrs_1_alpha_c;
  (*peh).beta_c = beta_c;
  (*peh).alpha_b = alpha_b;
  (*peh).beta_b = beta_b;
  (*peh).beta_node = beta_node;
  (*peh).k_peak = k_peak;
  (*peh).sound_horizon_fit = sound_horizon_fit;
  (*peh).alpha_gamma = alpha_gamma;
    
  return;
}

float TFfit_onek(float k, float *tf_baryon, float *tf_cdm, EH eh)
     /* Input: k -- Wavenumber at which to calculate transfer
        function, in Mpc^-1.
        *tf_baryon, *tf_cdm -- Input value not used; replaced on
        output if the input was not NULL. */
     /* Output: Returns the value of the full transfer function
        fitting formula.  This is the form given in Section 3 of
        Eisenstein & Hu (1997).
        *tf_baryon -- The baryonic contribution to the full fit.
        *tf_cdm -- The CDM contribution to the full fit. */
     /* Notes: Units are Mpc, not h^-1 Mpc. */
{
  // EH global vars 
  float omhh,obhh,theta_cmb,z_equality,k_equality,z_drag,R_drag,R_equality;
  float sound_horizon,k_silk,alpha_c,beta_c,alpha_b,beta_b,beta_node,k_peak;
  float sound_horizon_fit,alpha_gamma,mrs_1_alpha_c;
  // local vars
  float T_c_ln_beta, T_c_ln_nobeta, T_c_C_alpha, T_c_C_noalpha;
  float q, xx, xx_tilde;
  float T_c_f, T_c, s_tilde, T_b_T0, T_b, f_baryon, T_full;
  float mrstmp1,mrstmp2,mrstmp3,mrstmp4,mrstmp5;

  // set values from structure
  omhh = eh.omhh;
  obhh = eh.obhh;
  theta_cmb = eh.theta_cmb;
  z_equality = eh.z_equality;
  k_equality = eh.k_equality;
  z_drag = eh.z_drag;
  R_drag = eh.R_drag;
  R_equality = eh.R_equality;
  sound_horizon = eh.sound_horizon;
  k_silk = eh.k_silk;
  alpha_c = eh.alpha_c;
  mrs_1_alpha_c = eh.mrs_1_alpha_c;
  beta_c = eh.beta_c;
  alpha_b = eh.alpha_b;
  beta_b = eh.beta_b;
  beta_node = eh.beta_node;
  k_peak = eh.k_peak;
  sound_horizon_fit = eh.sound_horizon_fit;
  alpha_gamma = eh.alpha_gamma;

  k = fabs(k);	/* Just define negative k as positive */
  if (k==0.0) {
    if (tf_baryon!=NULL) *tf_baryon = 1.0;
    if (tf_cdm!=NULL) *tf_cdm = 1.0;
    return 1.0;
  }

  q = k/13.41f/k_equality;
  xx = k*sound_horizon;

  T_c_ln_beta = logf(2.718282f+1.8f*beta_c*q);
  T_c_ln_nobeta = logf(2.718282f+1.8f*q);
  mrstmp1 = 386.f/(1.f+69.9f*powf(q,1.08f)); // tmp variable for speed
  //  T_c_C_alpha = 14.2f/alpha_c + mrstmp1;
  T_c_C_alpha = 14.2f*mrs_1_alpha_c + mrstmp1; // for optimization
  T_c_C_noalpha = 14.2f + mrstmp1;

  //  T_c_f = 1.f/(1.f+POW4(xx/5.4f));
  mrstmp2 = xx/5.4f; // opt
  T_c_f = 1.f/(1.f+mrstmp2*mrstmp2*mrstmp2*mrstmp2); // opt
  T_c = T_c_f*T_c_ln_beta/(T_c_ln_beta+T_c_C_noalpha*SQR(q)) +
    (1-T_c_f)*T_c_ln_beta/(T_c_ln_beta+T_c_C_alpha*SQR(q));

  //  s_tilde = sound_horizon*powf(1.f+CUBE(beta_node/xx),-1.f/3.f);
  mrstmp3 = beta_node/xx; // opt
  s_tilde = sound_horizon*powf(1.f+mrstmp3*mrstmp3*mrstmp3,-1.f/3.f); // opt
  xx_tilde = k*s_tilde;

  T_b_T0 = T_c_ln_nobeta/(T_c_ln_nobeta+T_c_C_noalpha*SQR(q));
  //  T_b = sinf(xx_tilde)/(xx_tilde)*(T_b_T0/(1+SQR(xx/5.2f)) +
  //     alpha_b/(1.f+CUBE(beta_b/xx))*expf(-powf(k/k_silk,1.4f)));
  mrstmp4 = xx/5.2f; // opt
  mrstmp5 = beta_b/xx; // opt
  T_b = sinf(xx_tilde)/(xx_tilde)*(T_b_T0/(1+mrstmp4*mrstmp4) + 
        alpha_b/(1.f+mrstmp5*mrstmp5*mrstmp5)*expf(-powf(k/k_silk,1.4f)));//opt

  f_baryon = obhh/omhh;
  T_full = f_baryon*T_b + (1-f_baryon)*T_c;

  /* Now to store these transfer functions */
  if (tf_baryon!=NULL) *tf_baryon = T_b;
  if (tf_cdm!=NULL) *tf_cdm = T_c;
  return T_full;
}

