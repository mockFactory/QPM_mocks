#include <stdlib.h>
#include "header.h"


/* These are the globale variables. When needed, each is set to the
 * default value. Those without defualt values are either calculated
 * or set in the batfile.
 */


double GAMMA=0.2,         /* Shape parameter of EBW power spectrum */
  HUBBLE=0.7,             /* Hubble in 100 km/s/Mpc (not realy used in code). */
  SIGMA_8=0.95,           /* Normalization of power spectrum */
  SIGMA_8Z0=0.8,          /* Normalization of power spectrum */
  RHO_CRIT=2.78e11,       /* Critial mass density  in h^2 M_sun/Mpc^3 */
  SPECTRAL_INDX=1.0,      /* n_s -> P(k) = k^n_s */
  OMEGA_M=0.1,            /* Matter density */
  OMEGA_TEMP=0.3,         /* For M/L minimization */
  OMEGA_B=0.0,            /* Baryon density */
  DELTA_CRIT=1.686,       /* Critical overdensity for linear collapse */
  MSTAR,                  /* Mass scale at which sigm(M) = DELTA_CRIT */
  GALAXY_DENSITY,         /* Number density of galaxies (Mpc/h)^-3 */
  GALAXY_DENSITY2,        /* Number density of SECOND SET of galaxies (Mpc/h)^-3 (for x-corr)*/
  GRAVITY,                /* Newton's constant in internal units. */
  BOX_SIZE,               /* Size of box, if comparing to simulations. */
  RESOLUTION,             /* Simulations: BOX_SIZE/np^(1/3) */
  R_MIN_2HALO=-1.0,       /* Minimum scale of two-halo pairs, set by halo exclusion */
  MASS_PER_PARTICLE,      /* Simulations: mass in M_sol/h */
  GALAXY_BIAS,            /* Large scale bias of galaxies w.r.t linear matter distribution */
  DELTA_HALO=200,         /* Overdensity which defines the edge of a halo. */
  VBIAS=1.0,              /* Velocity bias of satellite galaxies. */
  VBIAS_SLOPE=0.00,       /* Slope of velocity bias relation with halo mass. */
  VBIAS_MASS_THRESHOLD=0, /* Threshold mass above which there is no velocity bias. */
  VBIAS_C=0.00,           /* Velocity bias of central galaxies. */
  CVIR_FAC=1.0,           /* Ratio between galaxy and dark matter concentrations */
  BIAS_A=0.707,           /* parameter for halo bias function (Tinker et all 2005 App. A) */
  BIAS_B=0.35,            /* parameter for halo bias function (adapted from Sheth Mo Tormen) */
  BIAS_C=0.8,             /* parameter for halo bias function */
  JENKINS_A=0.315,        /* Jenkins mass function--> normalization */
  JENKINS_B=0.61,         /* Jenkins mass function--> constant in exponential */
  JENKINS_C=3.8,          /* Jenkins mass function--> exponent in exponential */
  DNDM_PARAMS[10],        /* Holds parameters for 5-param dndM fit */
  MAXVEL=4000,            /* maximum velocity in the P(vz) lookup table (set later) */
  SIGV=500,               /* Free parameter for Kaiser model */
  BETA,                   /* Redshift-space distortion parameter. */
  XI_MAX_RADIUS=0,        /* Maximum radial value of 2-halo term */
  LOCAL_DENSITY=0,        /* for Hubble Bubble calculations--> change the mass function */
  MASS_THRESHOLD=0,       /* Mass at which you turn on change in HOD */
  REDSHIFT=0,             /* redshift */
  SATELLITE_FRACTION=0,   /* fraction of galaxies that are satellites in the HOD */
  DENSITY_THRESHOLD=0,    /* Density at which you turn on change in HOD. */
  HALO_EXCLUSION_RADIUS=2;/* 2== spherical exclusion. 1== no-overlap */

int ITRANS=4,             /* Type of transfer function to be used */
  RESET_KAISER=0,         /* Flag to reset tabulated Kaiser distortion quantities */
  RESET_COSMOLOGY=0,      /* Flag to recalculate tabulated quantities (e.g. b(M), dn/dM...) */
  SOFT_CENTRAL_CUTOFF=0,  /* Whether N_cen is a step function or not */
  NUM_POW2MASS_BINS,      /* Number of mass bins for velocity PDF calculation */
  RESET_PVZ=0,            /* Flag to recalculate velocity PDF table */
  COVAR=1,                /* Flag to use covariance matrix for wp chi^2 */
  PCA=0,                  /* Flag to use principle component analysis for chi^2 calc */
  COVARZ=0,               /* Flag to use covariance matrix for xi(s,p) chi^2*/
  EXCLUSION=2,            /* Type of two-halo exclusion ("spherical" is default) */
  FIX_PARAM=1,            /* Which parameter to leave free for fixing ng (1=M_min)*/
  DEPROJECTED=0,          /* Fit 3-space data rather than wp(rp) */
  OUTPUT=1,               /* Flag to output diagnostic information */
  POWELL=1,               /* minimization technique 1=POWELL, 0=AMOEBA */
  MCMC=0,                 /* Flag for doing Markov Chain analysis */
  LINEAR_PSP=0,           /* Flag for using linear P(k) for two-halo term. (Instead of Smith)*/
  KAISER=0,               /* Flag for using Kaiser linear+exponential model */
  ERROR_FLAG=0,           /* Global notification of problem with qromo/zbrent/etc */
  BEST_FIT=0,             /* Use best-fit mass function, bias from Tinker et al's analysis */
  RESET_FLAG_2H=0,        /* Flag to recalculate 2-halo term */
  RESET_FLAG_1H=0,        /* Flag to recalculate 1-halo term */
  IVFLAG=0,
  WP_ONLY=0,              /* =1 use only wp, =0 use n(N) and M/L (+others in ml_min.c) */
  GAO_EFFECT=0,           /* =1 means that HOD changes in high density, =0 means lo den change.*/
  N_HOD_PARAMS=9,         /* number of possible parameters for HOD. for use w/ MCMC */
  XCORR=0,                /* flag for cross-correlation */  
  RESTART = 0,            /* For restarting an MCMC chain (ml only) */
  USE_ERRORS = 0,         /* Flag for using systematic errors on massfunc, bias, etc in MCMC */
  DENSITY_DEPENDENCE = 0; /* Flag for den-dep HOD. */

/* Stuff for wtheta
 */
double MSHIFT,
  MALL,
  NGAL_DRG = 6.5e-4;

int *used_halo;

int ARGC;
char **ARGV;

char RESTART_FILE[100];   /* previous output file from which to restart chain */

long IDUM_MCMC=-555;      /* Random seed to start Markov Chain. */

/* These are for the parallel implementation, an ID for each
 * processor and the total number of processors.
 */
int ThisTask=0,
  NTask=1;

struct hod_parameters HOD,HOD2,HODt;
struct file_parameters Files;
struct perform_tasks Task;
struct workspace Work;
struct COLOR_DATA wp_color;
struct m2n_workspace M2N;
struct WPL wpl;

int LENSING_OUTPUT_FLAG = 0;   /* flag for outputting all data+model fits to a file during MCMC */
int DONT_FIT_SMF        = 0;   // don't use the smf in the fit
int DONT_FIT_LENSING    = 0;   // don't use the lensing in the fit
int DONT_FIT_CLUSTERING = 0;
int BLUE_FLAG           = 1;
int JPL_FLAG            = 0;
int RESET_FLAG_DS       = 1;
int RESET_FLAG_XGM1     = 1;
int FISHER              = 0;
int SMF_CORRECTION      = 0;
int RESET_HALO_BIAS     = 1; 
int ASCII_OUTPUT        = 0;
int NO_FOF_HALOS        = 0;
float SUBFRAC = -1;
double MASS_MIN = -1;
