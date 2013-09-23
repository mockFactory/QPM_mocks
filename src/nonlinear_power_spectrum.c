#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/* This is my implementation of the Halo Model fitting function
 * described in Appendix C of Smith, Peacock, et al. 2003 MNRAS 341, 1311
 *
 * Comparison to Smith et al. program halofit.f output is successful.
 */

/* Internal functions
 */
float func_nl(float r);
void tabulate_pk_nl(float *kk, float *pknl, int nk);

/* This returns the linear power 
 *   Delta = 4pi*k^3 P(k)/(2pi)^3
 */
float linear_power_spectrum(float xk)
{
  static float *kk,*pknl,*y2,pnorm=-1,ahi,bhi;
  static int flag=1,nk=1000,prev_cosmology=0;
  float a,psp,x1[4],y1[4];
  int i;

  if(pnorm<0 || prev_cosmology!=RESET_COSMOLOGY)
    {
      pnorm=SIGMA_8/sigmac(8.0);  
      pnorm*=pnorm;
      prev_cosmology=RESET_COSMOLOGY;
    }
  if(ITRANS>0)
    psp=pow(xk,SPECTRAL_INDX)*pow(transfnc(xk),2.);
  else
    psp=pow(xk,SPECTRAL_INDX);
  psp=psp*pnorm*xk*xk*xk/(2*PI*PI);
  return(psp);
}

float nonlinear_power_spectrum(float xk)
{
  static float *kk,*pknl,*y2,pnorm=-1,ahi,bhi;
  static int flag=1,nk=1000,prev_cosmology=0;
  float a,psp,x1[4],y1[4],xklog;
  int i;
  
  if(flag || RESET_COSMOLOGY!=prev_cosmology)
    {
      prev_cosmology=RESET_COSMOLOGY;
      pnorm=SIGMA_8/sigmac(8.0);   
      flag=0;
      kk=vector(1,nk);
      pknl=vector(1,nk);
      y2=vector(1,nk);
      tabulate_pk_nl(kk,pknl,nk);
      spline(kk,pknl,nk,1.0E+30,1.0E+30,y2);

      /* This takes the last four points in the power spectrum at high k
       * and fits a power law to them for extrapolation.
       */
      for(i=0;i<4;++i)
	{
	  x1[i]=(kk[nk-3+i]);
	  y1[i]=(pknl[nk-3+i]);
	}
      least_squares(x1,y1,4,&ahi,&bhi);

    }

  xklog=log(xk);

  /* If xk is less than the smallest k value tabulates, return linear power.
   */
  if(xklog<kk[1])
    return(linear_power_spectrum(xk));

  /* If xk larger than highest k, return extrapolation.
   */
  if(xklog>kk[nk])
    return((exp((ahi+bhi*xklog))));

  splint(kk,pknl,y2,nk,xklog,&a);
  return(exp(a));

}
void tabulate_pk_nl(float *kk, float *pknl, int nk)
{
  float r_nl,pnorm,DeltaQ,DeltaLin,DeltaH,psp,neff,s1,s2,n1,n2,ncurve,
    lnk_lo,lnk_hi,dlnk,xk1,y,xk,fy;
  float an,bn,cn,f1b,f2b,f3b,alpha,beta,gamma,nu,mu;
  int i,j,n=10000;

  /* First normalize the linear power spectrum.
   */
  pnorm=SIGMA_8/sigmac(8.0);

  /* Calculate the non-linear scale.
   */
  r_nl = exp(zbrent(func_nl,log(0.01),log(10.0),1.0E-4));
  if(OUTPUT)
    fprintf(stdout,"R_NL= %f\n",r_nl);

  /* Calculate the effective spectral index at the non-linear scale.
   */
  s1=pnorm*sigmac(-r_nl*0.999);
  s2=pnorm*sigmac(-r_nl*1.001);
  neff=-(3+2*(s2-s1)/0.002);
  if(OUTPUT)
    fprintf(stderr,"neff= %f\n",neff);

  /* Spectral curvature.
   */
  lnk_hi=10.0;
  lnk_lo=-10.0;
  dlnk=(lnk_hi-lnk_lo)/n;
  s1=0;
  for(i=1;i<=n;++i)
    {
      xk1=exp(lnk_lo+(i-0.5)*dlnk);
      y=xk1*r_nl;
      DeltaLin=linear_power_spectrum(xk1);
      s1+=DeltaLin*y*y*(1-y*y)*exp(-y*y)*dlnk;
    }
  ncurve=(3+neff)*(3+neff)+4*s1;
  if(OUTPUT)
    fprintf(stderr,"ncurve= %f\n",ncurve);

  /* Coefficients of the model.
   */
  an=pow(10.0,1.4861 + 1.8369*neff + 1.6762*neff*neff + 0.7940*neff*neff*neff + 
	 0.1670*pow(neff,4.0) - 0.6202*ncurve);
  bn=pow(10.0,0.9463 + 0.9466*neff + 0.3084*neff*neff - 0.9400*ncurve);

  cn= pow(10.0,-0.2807 + 0.6669*neff + 0.3214*neff*neff - 0.0793*ncurve);

  gamma =  0.8649 + 0.2989*neff + 0.1631*ncurve;
  alpha =  1.3884 + 0.3700*neff - 0.1452*neff*neff;
  beta  =  0.8291 + 0.9854*neff + 0.3401*neff*neff;
  mu    = pow(10.0,-3.5442 + 0.1908*neff);
  nu    = pow(10.0, 0.9589 + 1.2857*neff);

  /* Testing the parameter dependence.
   * TESTING TESTING
   */
  /*
  alpha *= 0.3;
  beta *= 0.3;
  */

  /* Omega-dependent functions (FLAT LAMBDA COSMOLOGY)
   */
  f1b = pow(OMEGA_M,-0.0307);
  f2b = pow(OMEGA_M,-0.0585);
  f3b = pow(OMEGA_M,+0.0743);


  /* Tabulate the power spectrum
   */
  lnk_lo=log(0.001);
  lnk_hi=log(1000.0);
  dlnk=(lnk_hi-lnk_lo)/(nk-1);

  for(i=1;i<=nk;++i)
    {
      xk=exp(lnk_lo+(i-1)*dlnk);
      y=xk*r_nl;
      fy=y/4.0 + y*y/8.0;
      
      /* TEST */
      /* fy*=1.2;  */

      DeltaLin=linear_power_spectrum(xk);
      
      DeltaQ = DeltaLin*pow(1+DeltaLin,beta)/(1+alpha*DeltaLin)*exp(-fy);

      DeltaH = an*pow(y,3*f1b)/(1+bn*pow(y,f2b)+pow(cn*f3b*y,3-gamma));

      DeltaH*= 1.0/(1+mu/y+nu/(y*y));

      kk[i]=log(xk);
      pknl[i]=log(DeltaQ + DeltaH);
    }
}

/* This is a function to find the non-linear scale. Similar to M_star, 
 * but here we're using a Guassian smoothing window rather than top hat.
 * (And we're finding the scale at which the variance = 1, not 1.686).
 */
float func_nl(float r)
{
  static int prev_cosmology=0;
  float sig;
  static float pnorm=-1;

  if(pnorm<0 || RESET_COSMOLOGY!=prev_cosmology)
    pnorm=SIGMA_8/sigmac(8.0);
  prev_cosmology=RESET_COSMOLOGY;

  r=exp(r);
  sig=pnorm*sigmac(-r);
  return sig-1.0;
}
