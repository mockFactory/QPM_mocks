#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/* This is a calculation of the differential halo mass function of
 * Jenkins et al 2001 (MNRAS 321, 372). 
 *
 */
void initialize_mass_function(float *a1, float *a2, float *a3, float *a4);

float halo_mass_function(float mass)
{
  float sig,logm,a,slo,shi,rm,rlo,rhi,mlo,mhi,dsdM,n,nuprime,nufnu,p,A;
  static int flag=0,SO180=0,SO324=0,WARREN=0,ST=0,JENKINS=0;
  static float pnorm, prev_delta, prev_cosmo;
  float btemp = -1;

  static float
    a1 = 0.325277,
    a2 = 0.492785,
    a3 = 0.310289,
    a4 = 1.317104,
    a5 = 2.425681;
    
  /* Jenkins et al. SO 180 best-fit
   */
  if(SO180)
    {
      JENKINS_A = 0.301;
      JENKINS_B = 0.64;
      JENKINS_C = 3.82;
    }      
  if(SO324)
    {
      JENKINS_A = 0.316;
      JENKINS_B = 0.67;
      JENKINS_C = 3.82;
    }      
  if(JENKINS) //their .2 fof function
    {
      JENKINS_A = 0.315;
      JENKINS_B = 0.61;
      JENKINS_C = 3.8;
    }      


  /* First normalize the power spectrum
   */
  pnorm=SIGMA_8/sigmac(8.0);
  rm=pow(3.0*mass/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  sig=pnorm*sigmac(rm);
  logm=log10(mass);
  
  mlo=0.99*mass;
  mhi=1.01*mass;
  rlo=pow(3.0*mlo/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  rhi=pow(3.0*mhi/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);

  slo=pnorm*sigmac(rlo);
  shi=pnorm*sigmac(rhi);
  dsdM=(shi-slo)/(mhi-mlo);

  if(SO324)goto JENKINS_FUNCTION;
  if(SO180)goto JENKINS_FUNCTION;
  if(WARREN)goto WARREN_FUNCTION;
  if(ST)goto ST_FUNCTION;
  if(JENKINS)goto JENKINS_FUNCTION;

  /* Tinker et al. (in prep) for SO 200
   */
  if(DELTA_HALO != prev_delta || prev_cosmo != RESET_COSMOLOGY)
    {
      initialize_mass_function(&a1,&a2,&a3,&a4);
      prev_delta = DELTA_HALO;
      prev_cosmo = RESET_COSMOLOGY;

      // if we're using systematic errors in an MCMC, adjust parameter a1 (amplitude)
      if(USE_ERRORS)
	a1 *= M2N.mf_amp;
      fprintf(stderr,"MF PARAMS for DELTA=%f %f %f %f %f\n",DELTA_HALO,a1,a2,a3,a4);
    }
  
  n = -a1*(pow(sig/a3,-a2)+1)*exp(-a4/sig/sig)*OMEGA_M*RHO_CRIT/mass/sig*dsdM;
  return(n);

  /* Jenkins et al. FOF .2 best-fit (unless SO180==1)
   */
 JENKINS_FUNCTION:
  a=-JENKINS_A*OMEGA_M*RHO_CRIT/mass/sig;
  n=a*dsdM*exp(-pow(fabs(JENKINS_B-log(sig)),JENKINS_C));
  return(n);

  /* Warren et al. (calibrated only on concordance cosmology, FOF.2)
   */
 WARREN_FUNCTION:
  n = -0.7234*(pow(sig,-1.625)+0.2538)*exp(-1.198/sig/sig)*OMEGA_M*RHO_CRIT/mass/sig*dsdM;
  return(n);



  /* Need to find the derivative dlog(sig)/dlog(M)
   */
  mlo=0.99*logm;
  mhi=1.01*logm;
  rlo=pow(3.0*pow(10.0,mlo)/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  rhi=pow(3.0*pow(10.0,mhi)/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  slo=log10(pnorm*sigmac(rlo));
  shi=log10(pnorm*sigmac(rhi));
  dsdM=(shi-slo)/(mhi-mlo);

 ST_FUNCTION:

  /* This is a bunch of Sheth-Tormen stuff.
   * NB! because I'm skipping the above derivative (dlogs/dlogM), i'm using the lower
   */
  nuprime=0.841*DELTA_CRIT/sig;
  nufnu=0.644*(1+1.0/pow(nuprime,0.6))*(sqrt(nuprime*nuprime/2/PI))*exp(-nuprime*nuprime/2);
  //n=RHO_CRIT*OMEGA_M/mass*mass*nufnu*fabs(dsdM);
  n=RHO_CRIT*OMEGA_M/mass*nufnu*fabs(dsdM)/sig;
  return(n);




}


/* It may be a bit costly to run the above function every time you need
 * dn/dM, so here we put the values into an array and then interpolate. 
 *
 * The currentrange of masses calculated is 10^9 to 10^16.7. The tabulation is
 * done in log(M), so the spline interpolation will perform a power-law fit
 * to masses outside this range.
 */
float dndM_interp(float m)
{
  static int flag=0,prev_cosmo=0, n;
  static float *x,*y,*y2;
  int i;
  float dm,max=16.7,min=8,a,m1,m2,dm1;

  if(!flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      n = 200;
      if(!ThisTask && OUTPUT)
	fprintf(stdout,"RESET: resetting mass function for %f %f\n",OMEGA_M,SIGMA_8);
      fflush(stdout);

      if(!flag)
	{
	  x=vector(1,n);
	  y=vector(1,n);
	  y2=vector(1,n);
	}
      flag=1;
      dm=(float)(max-min)/n;
      for(i=1;i<=n;++i)
	{
	  x[i]=pow(10.0,min+i*dm);
	  y[i]=log(halo_mass_function(x[i]));
	  //printf("MF%d %e %e\n",RESET_COSMOLOGY,x[i],exp(y[i]));fflush(stdout);
	  if(isnan(y[i])) { n = i-1; break; }
	  if(isinf(y[i])) { n = i-1; break; }
	  x[i]=log(x[i]);
	  continue;
	}
      spline(x,y,n,2.0E+30,2.0E+30,y2);
      prev_cosmo=RESET_COSMOLOGY;
      //fprintf(stderr,"MMAX %e\n",exp(x[n]));
    }
  m=log(m);
  if(m>x[n])return 0;
  splint(x,y,y2,n,m,&a);
  return(exp(a));

}



void initialize_mass_function(float *a1, float *a2, float *a3, float *a4)
{
  int n = 9, i;
  float *x, *y, *z, at, ztemp;

  x = vector(1,n);
  y = vector(1,n);
  z = vector(1,n);

  // initialize the overdensities
  for(i=1;i<=9;i+=2)
    x[i] = log(200*pow(2.0,(i-1.0)/2.0));
  for(i=2;i<=9;i+=2)
    x[i] = log(300*pow(2.0,(i-2.0)/2.0));

  //first parameter
  y[1] = 1.858659e-01 ;
  y[2] = 1.995973e-01 ;
  y[3] = 2.115659e-01 ;
  y[4] = 2.184113e-01 ;
  y[5] = 2.480968e-01 ;
  y[6] = 2.546053e-01 ;
  y[7] = 2.600000e-01 ;
  y[8] = 2.600000e-01 ;
  y[9] = 2.600000e-01 ;

  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a1);
  if(DELTA_HALO>=1600) *a1 = 0.26;

  //second parameter
  y[1] = 1.466904e+00 ;
  y[2] = 1.521782e+00 ;
  y[3] = 1.559186e+00 ;
  y[4] = 1.614585e+00 ;
  y[5] = 1.869936e+00 ;
  y[6] = 2.128056e+00 ;
  y[7] = 2.301275e+00 ;
  y[8] = 2.529241e+00 ;
  y[9] = 2.661983e+00 ;

  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a2);

  //third parameter
  y[1] = 2.571104e+00 ;
  y[2] = 2.254217e+00 ;
  y[3] = 2.048674e+00 ;
  y[4] = 1.869559e+00 ;
  y[5] = 1.588649e+00 ;
  y[6] = 1.507134e+00 ;
  y[7] = 1.464374e+00 ;
  y[8] = 1.436827e+00 ;
  y[9] = 1.405210e+00 ;

  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a3);


  //fourth parameter
  y[1] = 1.193958e+00;
  y[2] = 1.270316e+00;
  y[3] = 1.335191e+00;
  y[4] = 1.446266e+00;
  y[5] = 1.581345e+00;
  y[6] = 1.795050e+00;
  y[7] = 1.965613e+00;
  y[8] = 2.237466e+00;
  y[9] = 2.439729e+00;

  spline(x,y,n,1.0E+30,1.0E+30,z);
  splint(x,y,z,n,log(DELTA_HALO),a4);

  // now adjust for redshift
  if(!(REDSHIFT>0))return;

  ztemp = REDSHIFT;
  if(REDSHIFT>3) ztemp = 3.0;
  *a1 *= pow(1+ztemp,-0.14);
  *a2 *= pow(1+ztemp,-0.14);
  at = -pow(0.75/log10(DELTA_HALO/75),1.2);
  at = pow(10.0,at);
  *a3 *= pow(1+ztemp,-at);

}

