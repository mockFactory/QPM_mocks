#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

float collapse_redshift(float z);
float cvir_pnorm_g1,
  cvir_sigma_g1;

/* This calculates and tabulates the halo concentrations
 * as a function of halo mass. Uses the "Bullock model", 
 * described in a little more detail below.
 */

float halo_concentration(float m)
{
  static int flag=1,n,prev_cosmo=0;
  static float *x,*y,*y2;
  int i;
  float x1,x2,cfac;
  float a,dm,x3,x4;
  FILE *fp;
  char fname[1000];

  if(flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      MSTAR = mstar();
      if(OUTPUT)
	fprintf(stdout,"Calc cvir with DELTA_HALO= %f\n",DELTA_HALO);
      prev_cosmo=RESET_COSMOLOGY;
      n=50;
      if(flag)
	{
	  x=vector(1,n);
	  y=vector(1,n);
	  y2=vector(1,n);
	}
      cvir_pnorm_g1=SIGMA_8/sigmac(8.0);

      flag=0;
      
      dm=(log(HOD.M_max)-log(1.0E8))/(n-1);
      for(i=1;i<=n;++i)
	{
	  x[i]=exp((i-1)*dm+log(1.0E8));
	  // ALEXIE changed to new model here
	  //y[i]=cvir_model(x[i]);
	  //printf("first : %f \n",y[i]);
	  y[i]=munoz_cuartas_cvir_model(x[i]);
	  //printf("second : %f \n \n",y[i]);
	  x[i]=log(x[i]);
	  y[i]=log(y[i]);
	}
      spline(x,y,n,1.0E+30,1.0E+30,y2);
    }
  m=log(m);
  splint(x,y,y2,n,m,&a);
  return(exp(a));
}



float munoz_cuartas_cvir_model(float mass)
{
  float alpha, beta, gamma, b_cvir, a_cvir, w_cvir, m_cvir,cvir;
  
  w_cvir = 0.029;
  m_cvir = 0.097;
  alpha  = -110.001;
  beta   = 2469.720;
  gamma  = 16.885;

  b_cvir = (alpha/(REDSHIFT+gamma))+(beta/pow(REDSHIFT+gamma,2));
  a_cvir = w_cvir*REDSHIFT - m_cvir;

  cvir = pow(10.0, a_cvir*log10(mass)+b_cvir);
  //fprintf("cvir and mass %e %e \n",cvir, mass);

  return(cvir);
}


