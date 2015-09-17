#ifndef EH_HINCLUDED
#define EH_HINCLUDED

typedef struct eh_globals {
  int settf;
  float omhh,          /* Omega_matter*h^2 */
    obhh,              /* Omega_baryon*h^2 */
    theta_cmb,         /* Tcmb in units of 2.7 K */
    z_equality,        /* Redshift of matter-radiation equality, really 1+z */
    k_equality,        /* Scale of equality, in Mpc^-1 */
    z_drag,            /* Redshift of drag epoch */
    R_drag,            /* Photon-baryon ratio at drag epoch */
    R_equality,        /* Photon-baryon ratio at equality epoch */
    sound_horizon,     /* Sound horizon at drag epoch, in Mpc */
    k_silk,            /* Silk damping scale, in Mpc^-1 */
    alpha_c,           /* CDM suppression */
    mrs_1_alpha_c,     // recipricol of alpha_c (for optimization)
    beta_c,            /* CDM log shift */
    alpha_b,           /* Baryon suppression */
    beta_b,            /* Baryon envelope shift */
    beta_node,         /* Sound horizon shift */
    k_peak,            /* Fit to wavenumber of first peak, in Mpc^-1 */
    sound_horizon_fit, /* Fit to sound horizon, in Mpc */
    alpha_gamma;       /* Gamma suppression in approximate TF */
} EH;

#define SQR(a) ((a)*(a))
#define CUBE(a) ((a)*(a)*(a))
#define POW4(a) ((a)*(a)*(a)*(a))

int mrstfset(float omegam, float omegab, float h, EH *peh);
float mrstf(float k, EH *peh);
void TFset_parameters(float omega0hh, float f_baryon, float Tcmb, EH *peh);
float TFfit_onek(float k, float *tf_baryon, float *tf_cdm, EH eh);

#endif
