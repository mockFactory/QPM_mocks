#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

/* This file holds several useful functions for the HOD, such as the number of
 * central and satellite halos, their second moment.
 *
 * This also intializes the HOD parameters.
 */

/* internal routines.
 */
float fixfunc1(float m);
float fixfunc2(float m);
float func_mlow(float m);
float func_mhi(float m);
float one_halo_from_file(float m);

/*****************
 * DIFFERENT HOD PARAMETERIZATIONS: (1) Central galaxies
 *
 * This integer HOD.pdfc controls the parameterization of the central occupation.
 * (I realize 'pdf' is not the right acronym, but it's a holdover from old code.)
 *
 * For soft central cutoffs, the PDF is the nearest integer distribution, AKA
 * Bernoulli distribution.
 *
 * 0 = No galaxies, just halos. (always returns 1)
 * 1 = Hard cutoff at M_min. Unity above M_min. (DEFAULT)
 * 2 = Soft cutoff of the form 0.5*(1+erf((log10(m) - log10(HOD.M_min))/HOD.sigma_logM)
 * 3 = Soft cutoff of the form N_cen = exp(-HOD.M_min/m)
 * 4 = Hard cutoff, but N_cen -> 0 for M>M_min (i.e., blue galaxies) 
 *      MaxCen*exp(-(lgM-lgM_min)**2/2/(sigma_logM)**2)
 *      NB! -> The value of N_cen(M_min) can be < 1. (It will== MaxCen)
 * 5 = Hard Cutoff, Step function, but N_cen != 1, but some number < 1.
 * 6 = Same as 4, but symmetric Gaussian around M_min (This is for mag bins.)
 * 7 = A sqaure function. N_cen=1 for M_min <= m <= M_cen_max (for mag bin data)
 * 8 = A sqaure function (like case 7) but with Gaussian cutoffs on either edge
 *      instead of sharp cutoffs.
 *
 * 9 = Magnitude bin model, but the upper mass cutoff is defined by the lower
 *      mass cutoff of the next-highest bin. (see function in ml_minimization.c)
 *
 * M_low is a parameter created to keep the integrals from having to go to M=0 for
 * soft central cutoffs. N_cen(mlow) = 1.0E-3;
 * If the value of sigma_logM gets above 1 or so, then M_low can be rediculously 
 * small, so I have placed a lower limit on M_low of 1.0E+7 Msol/h
 */
float N_cen(float m)
{
  float x,f=1;

  switch(HOD.pdfc) {
  case -1:
    return 0;
  case 0:
    return 1;
  case 1: 
    if(m<HOD.M_min)return(0);
    return(f*1.0);
    break;
  case 2:
    //return((1-0.3*0.5)*HOD.MaxCen*0.5*(1+erf((log10(m) - log10(HOD.M_min))/HOD.sigma_logM)));
    return(f*HOD.MaxCen*0.5*(1+erf((log10(m) - log10(HOD.M_min))/HOD.sigma_logM)));
    break;
  case 3:
    return(f*exp(-HOD.M_min/m));
    break;
  case 4:
    if(m<HOD.M_min)return(0);
    x = (log10(m) - log10(HOD.M_min))/HOD.sigma_logM;
    return(f*HOD.MaxCen*exp(-x*x/2));
  case 5:
    if(m<HOD.M_min)return(0);
    return(f*HOD.MaxCen);
  case 6:
    x = (log10(m) - log10(HOD.M_min))/HOD.sigma_logM;
    return(f*HOD.MaxCen*exp(-x*x/2));    
  case 7:
    if(m>=HOD.M_min && m<=HOD.M_cen_max)return(f);
    return(0);
  case 8:
    if(m>=HOD.M_min && m<=HOD.M_cen_max)return(f);
    if(m<HOD.M_low)return(0);
    if(m<HOD.M_min)
      x = (log10(m) - log10(HOD.M_min))/HOD.sigma_logM;
    else
      x = (log10(m) - log10(HOD.M_cen_max))/HOD.sigma_logM;
    return(f*exp(-x*x/2));    
  case 10: // for DRGs
    if(m<HOD.M_min)return 0;
    return(f*exp(-HOD.M_min*HOD.mass_shift/(m-HOD.M_min)));
    //return(f*exp(-pow(HOD.M_min*HOD.mass_shift/(m-HOD.M_min),HOD.shift_alpha)));
    break;
  case 11: // for 1-DRGs
    if(m<HOD.M_min)return 0;
    return 1 - (f*exp(-HOD.M_min*HOD.mass_shift/(m-HOD.M_min)));
    break;

  case 20: // for DRGs
    if(m<HOD.M_min*HOD.mshift2)return 0;
    return(f*exp(-HOD.M_min*HOD.mass_shift*HOD.mshift2/(m-HOD.M_min*HOD.mshift2)));
    //return(f*exp(-pow(HOD.M_min*HOD.mass_shift/(m-HOD.M_min),HOD.shift_alpha)));
    break;
  case 21: // for 1-DRGs
    if(m<HOD.M_min)return 0;
    if(m<HOD.M_min*HOD.mshift2)return 1;
    return 1 - (f*exp(-HOD.M_min*HOD.mass_shift*HOD.mshift2/(m-HOD.M_min*HOD.mshift2)));
    break;

  default:
    endrun("Illegal value of HOD.pdfc.");
  }
  return 0;
}

/*****************
 * DIFFERENT HOD PARAMETERIZATIONS: (1) Satellite galaxies
 *
 * This integer HOD.pdfs controls the parameterization of the satellite occupation.
 *
 * 0 = halos only, no galaxies (always returns zero)
 * 1 = power law: N_sat = pow(m/HOD.M1,HOD.alpha) [cut off at M_low]
 * 2 = power law with soft cutoff: N_sat = pow((m-HOD.M_cut)/HOD.M1,alpha) 
 * 3 = power law with exp cutoff: N_sat = pow(m/HOD.M1,alpha)*exp(-M_cut/(m-M_min))
 * 4 = broken power law with exp cutoff (alpha changes at some high mass value)
 * 5 = broken power law (m-mcut)/m1 parameterization
 * 6 = power with alternate exp cutoff: N_sat = pow(m/M1)*exp(-M_min/m) -> Alexie changed this
 * 7 = from Zheng's new paper lognormal cutoff power law
 * 8 = power with Conroy exp cutoff: N_sat = pow(m/M1)*exp(-M_cut/m)
 * 
 * The PDF is always assumed to be Poisson for satellite galaxies.
 */
float N_sat(float m)
{
  float m1,f=1,n,nc;

  switch(HOD.pdfs) {
  case 0:
    return 0;
  case 1: 
    if(m<HOD.M_min)return(0);
    return(f*pow(m/HOD.M1,HOD.alpha));
    break;
  case 2:
    if(m<HOD.M_low || m<HOD.M_cut)return(0);
    return(f*pow((m-HOD.M_cut)/HOD.M1,HOD.alpha));
    break;
  case 3:
    if(m<HOD.M_min)return(0);
    return(f*exp(-HOD.M_cut/(m-HOD.M_min))*pow(m/HOD.M1,HOD.alpha));
    break;
  case 4:
    if(m<HOD.M_min)return(0);
    if(m<HOD.M_sat_break)
      return(f*exp(-HOD.M_cut/(m-HOD.M_min))*pow(m/HOD.M1,HOD.alpha));
    else
      {
	m1 = exp(log(HOD.M_sat_break) - HOD.alpha/HOD.alpha1*log(HOD.M_sat_break/HOD.M1));
	/*
	m1 = exp(log(HOD.M_sat_break) - 1/HOD.alpha1*
		 (HOD.alpha*log(HOD.M_sat_break/HOD.M1) - HOD.M_cut/(HOD.M_sat_break-HOD.M_min)));
	*/
	return(f*exp(-HOD.M_cut/(m-HOD.M_min))*pow(m/m1,HOD.alpha1));
      }
    break;
  case 5:
    if(m<HOD.M_low || m<HOD.M_cut)return(0);
    if(m<HOD.M_sat_break)
      return(f*pow((m-HOD.M_cut)/HOD.M1,HOD.alpha));
    else
      {
	m1 = HOD.M_sat_break*pow((HOD.M_sat_break-HOD.M_cut)/HOD.M1,-HOD.alpha/HOD.alpha1);
	return(f*pow(m/m1,HOD.alpha1));
      }
    break;
  case 6:
    return(pow(m/HOD.M1,HOD.alpha)*exp(-HOD.M_cut/m)); // Alexie changed from Mmin to Mcut here
  case 66:
    return(pow(m/HOD.M1,HOD.alpha)*exp(-(HOD.M_cut+HOD.M_min)/m)); // JEREMY CHANGED FROM BELOW!!
    //return(pow(m/HOD.M1,HOD.alpha)*exp(-HOD.M_cut/m)*Ncen_CLF(m)); // Alexie changed from Mmin to Mcut here
  case 7:
    if(m<HOD.M_low || m<HOD.M_cut)return(0);
    return(0.5*(1+erf((log10(m) - log10(HOD.M_min))/HOD.sigma_logM))
           *pow((m-HOD.M_cut)/HOD.M1,HOD.alpha));
    break;
  case 8:
    n = (pow(m/HOD.M1,HOD.alpha)*exp(-HOD.M_cut/m));
    if(m<HOD.M_min)
      {
	nc = N_cen(m);
	  if(n>nc)return(nc); 
      }
    return(n);
  case 9:
    if(m<HOD.M_sat_break)
      {
	n = (pow(m/HOD.M1,HOD.alpha)*exp(-HOD.M_cut/m));
	if(m<HOD.M_min)
	  {
	    nc = N_cen(m);
	    if(n>nc)return(nc); 
	  }
	return n;
      }
    m1 = (pow(HOD.M_sat_break/HOD.M1,HOD.alpha)*exp(-HOD.M_cut/HOD.M_sat_break));
    n = m1*pow(m/HOD.M_sat_break,HOD.alpha1);
    return n;
    break;

  case 10:
    if(m<HOD.M_low)return(0);
    return(0.5*(1+erf((log10(m) - log10(HOD.M_min))/HOD.sigma_logM))
           *pow((m)/HOD.M1,HOD.alpha));
    break;
  case 11:
    if(m<HOD.M_low)return(0);
    return(f*exp(-HOD.M_cut/m)*pow(m/HOD.M1,HOD.alpha)*
	   0.5*(1+erf((log10(m) - log10(HOD.M_min))/HOD.sigma_logM)));
    break;

  case 101:
    if(m<HOD.M_min)return 0;
    n = pow(m/HOD.M1,HOD.alpha);
    if(n<1)return 0;
    return n-1;

  case 102:
    if(m<HOD.M_min)return 0;
    return pow(m/HOD.M1,HOD.alpha);

  default:
    endrun("Illegal value of HOD.pdfs.");
  }
  return 0;
}

/* If the sample is split up by color, this function
 * returns the blue fraction at a given mass for satellite galaxies.
 * This function is parameterized as a log-normal.
 * See equation (11) in Zehavi et al Apj 2005, 360, 1
 *
 * if HOD.color == 1, returns blue fraction
 * if HOD.color == 2, returns red fraction
 */
float satellite_blue_fraction(float m)
{
  float x;
  x = log10(m) - log10(HOD.M_low0);
  x=(HOD.fblue0_sat*exp(-x*x/(2*HOD.sigma_fblue_sat*HOD.sigma_fblue_sat)));
  if(x>1)x=1;
  if(HOD.color==2)
    return(1-x);
  return(x);
}

/* If the sample is split up by color, this function
 * returns the blue fraction at a given mass for central galaxies.
 * This function is parameterized as a log-normal.
 * See equation (11) in Zehavi et al Apj 2005, 360, 1
 *
 * if HOD.color == 1, returns blue fraction
 * if HOD.color == 2, returns red fraction
 */
float central_blue_fraction(float m)
{
  float x;
  x = log10(m) - log10(HOD.M_low0);
  x=(HOD.fblue0_cen*exp(-x*x/(2*HOD.sigma_fblue_cen*HOD.sigma_fblue_cen)));
  if(x>1)x=1;
  if(HOD.color==2)
    return(1-x);
  return(x);
}

/* If what is needed is the total number of galaxies in a halo.
 */
float N_avg(float m)
{
  return(N_cen(m)+N_sat(m));
}


/* This is the <(N_sat-1)(N_sat)> moment, which is
 * for the number of pairs of satellite galaxies.
 * For a Poisson distribution, <N(N-1)> = N^2
 */
float moment_ss(float m)
{
  float n,x;
  n=N_sat(m);
  return(n*n);

  // sub-poisson model from millennium run
  // doesn't seem to effect -20.5 stats much
  x = 0.6 + 0.4*exp(-0.1/pow(n,1.5));
  return(n*n*x*x);

  
  // HAMANA et al
  if(n>1) return n*n;
  if(n>0.25) return n*n*log10(4*n)/log10(4);
  return 0;

}

/* This is a function to set the HOD parameters until 
 * I get some input code or batch file set up.
 * 
 */
void set_HOD_params()
{
  int i,j=1;
  float m,error=1.0,tol=1.0E-4,prev,s1,mlo;

  /* If the mass function at M_max is undefined, reset M_max
   */
  //LOCAL_DENSITY = -0.2;
  if(LOCAL_DENSITY!=0)
    {
      HOD.M_max = 1.0E16;
      while(dndM_interp(HOD.M_max)<=0) 
	{
	  HOD.M_max*=0.9;
	}
      fprintf(stderr,"NEW M_max= %e\n",HOD.M_max);
    }

  if(HOD.pdfc == 2 || HOD.pdfc == 3 || HOD.pdfc == 6 || HOD.pdfc == 8 || HOD.pdfc == 9)
    SOFT_CENTRAL_CUTOFF=1;

  /* Error trap both the galaxy density and M_min both left unspecified.
   */
  if(HOD.M_min<=0 && GALAXY_DENSITY<=0 && HOD.free[0]==0)
    endrun("ERROR: Must specify either M_min or GALAXY_DENSITY");

  /* If the user has specified M_min and M1, calculate the galaxy density.
   */
  if(HOD.M_min>0 && HOD.M1>0)
    {
      HOD.M_low = -1;
      if(SOFT_CENTRAL_CUTOFF)
	HOD.M_low0 = HOD.M_low = exp(zbrent(func_mlow,log(HOD.M_min*1.0E-6),log(HOD.M_min*1.1),1.0E-5));
      else
	HOD.M_low0 = HOD.M_low = HOD.M_min;
      GALAXY_DENSITY=qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
      if(OUTPUT) {
	fprintf(stdout,"M_low= %e\n",HOD.M_low);
	fprintf(stdout,"ng= %e\n",GALAXY_DENSITY); }
    }      

  /* If M_min<=0 then use the specified galaxy density to calculate M_min
   */
  if(HOD.M_min<=0)
    {
      if(HOD.pdfc==7 && HOD.pdfc==8)
	HOD.M_min=pow(10.0,zbrent(fixfunc1,8.0,log10(HOD.M_cen_max*0.99),1.0E-5));
      else
	HOD.M_min=pow(10.0,zbrent(fixfunc1,7.0,14.8,1.0E-5));

      HOD.M_low = -1;
      if(SOFT_CENTRAL_CUTOFF)
	HOD.M_low = exp(zbrent(func_mlow,log(HOD.M_min)-5*HOD.sigma_logM*2.3,
			       log(HOD.M_min),1.0E-5));
      else
	HOD.M_low = HOD.M_min;
      if(HOD.M_low<1.0E7)HOD.M_low=1.0E+7;
      if(OUTPUT) {
	fprintf(stdout,"M_min %e [ng= %e]\n",HOD.M_min,GALAXY_DENSITY);
	fprintf(stdout,"M_low= %e\n",HOD.M_low); }
    }

  /* If M1<=0 then use the specified galaxy density to calculate M1
   */
  if(HOD.M1<=0)
    {
      HOD.M_low = -1;
      if(SOFT_CENTRAL_CUTOFF)
	HOD.M_low = exp(zbrent(func_mlow,log(HOD.M_min)-5*HOD.sigma_logM*2.3,
			       log(HOD.M_min*1.1),1.0E-5));
      else
	HOD.M_low = HOD.M_min;
      if(HOD.M_low<1.0E7)HOD.M_low=1.0E+7;
      HOD.M1=pow(10.0,zbrent(fixfunc2,log10(HOD.M_low),15.8,1.0E-5));
      if(OUTPUT) {
	fprintf(stdout,"M1 %e [ng= %e]\n",HOD.M1,GALAXY_DENSITY);
	fprintf(stdout,"M_min = %e M_low= %e\n",HOD.M_min,HOD.M_low); }
    }

  /* Set the number of halo mass bins we've got.
   */
  NUM_POW2MASS_BINS=log(HOD.M_max/HOD.M_low)/LN_2+1;

  if(HOD.pdfc==6)
    HOD.M_hi = set_high_central_mass();
  HOD.M_low0 = set_low_mass();
  return;

}

/* If soft central cutoff, then put a lower limit on mass integrals
 * that begins at the mass where N_cen=0.001.
 */
float func_mlow(float m)
{
  /* Have a check in case the passed mass is equal to M_min, but the value of
   * N_cen is < 0.001 (which might be the case for very small sigma_logM)
   */
  //fprintf(stdout,"MLOTOP %e %e %e\n",exp(m),HOD.M_min,N_cen(exp(m)));
  if(fabs(exp(m)-HOD.M_min)<0.001*HOD.M_min)
    if(N_cen(exp(m))<0.001*N_cen(HOD.M_min))return(0);
  //fprintf(stdout,"MLO %e %e %e %e %e\n",exp(m),N_cen(exp(m)),HOD.M_min,N_cen(HOD.M_min),HOD.M_low); 
  //return(N_cen(exp(m))-0.001); //changing to 0.001 of the value at HOD.M_Min
  return(N_cen(exp(m))/N_cen(HOD.M_min)*1000 - 1);
}

/* this is the brute force, where i just go from M_min and step down.
 */
float mlow_brute_force()
{
  int i=0;
  float mass;
  mass = HOD.M_min;
  while(N_cen(mass)>1.0E-3) { mass*=0.9; }
  return mass;
}

/* It is straightforward to calculate what M_low
 * should be given the other parameters on N_cen.
 */
float set_low_mass()
{
  int i;
  float m,test_n;

  if(ERROR_FLAG)
    {
      fprintf(stdout,"uncaught ERROR at top of set low mass\n");
      ERROR_FLAG = 0;
    }

  if(!SOFT_CENTRAL_CUTOFF)return(HOD.M_min);
  switch(HOD.pdfc){
  case 8:
  case 6:
    m = log10(HOD.M_min) - sqrt(-2*HOD.sigma_logM*HOD.sigma_logM*log(0.001));
    m = pow(10.0,m);
    return(m);
  case 9:
    m = exp(zbrent(func_mlow,log(HOD.M_min)-5*HOD.sigma_logM*2.3,
		   log(HOD.M_min),1.0E-5));
    return(m);
  case 100:
  case 101:
    if(wpl.mlow_flag)
      return HOD.M_max*0.98;
    /*
    test_n=N_cen(HOD.M_min*1.0E-6);
    if(test_n > 0.001)
      return HOD.M_min*1.0E-7; // this is not ideal but should avoid getting some zbrent error here (notice returning 10^-7 here)
    */
    m = exp(zbrent(func_mlow,log(HOD.M_min*1.0E-6),
		   log(HOD.M_min),1.0E-5)); // changed to 1.0E-6. Could be a problem here is BOTH the data points are above 10^-3
    printf("LOW %e %e\n",HOD.M_min,m);
    if(m>HOD.M_max)
      {
	if(ERROR_FLAG) ERROR_FLAG = 0;
	return HOD.M_max*1.0E-3;
      }    
    if(ERROR_FLAG && HOD.M_min > 1.0E+14)
      {
	ERROR_FLAG = 0;
	return HOD.M_min*1.0E-6;
      }
    if(ERROR_FLAG)
      {
	ERROR_FLAG = 0;
	return HOD.M_min*1.0E-3;
      }
    if(m>HOD.M_min)
      {
	return HOD.M_min*1.0E-3;
      }
    
    if(ERROR_FLAG)
      {
	fprintf(stdout," ERROR in set_low_mass: %e %e %e %d \n",HOD.M_min, HOD.sigma_logM, GALAXY_DENSITY,wpl.mlow_flag);
	fprintf(stdout,"n wpl.mlow_flag = %d \n",wpl.mlow_flag);
	test_n=N_cen(HOD.M_max);
	fprintf(stdout,"ncen at HOD max = %e \n",test_n);
	fprintf(stdout,"ncen at HOD.M_min*1.0E-6 = %e \n",N_cen(HOD.M_min*1.0E-6));
	fprintf(stdout,"ncen at HOD.M_min and Mmin = %e %e \n",N_cen(HOD.M_min),HOD.M_min);
	for(i=1000;i<=1600;++i)
	 printf("HODCEN %f %e\n",i/100.0,N_cen(pow(10.0,i/100.0)));
	exit(0);
      }

    return(m);
  default:
    m = exp(zbrent(func_mlow,log(HOD.M_min*1.0E-6),
		   log(HOD.M_min),1.0E-5));
    return(m);
  }
  return(0);
}

/* If modeling magnitude bin samples, then there will be a high mass
 * scale for central galaxies as well as a low mass scale. This finds
 * the mass at which N_cen(m)=0.001, where m>M_min.
 */
float set_high_central_mass()
{
  float m,n;

  if(HOD.pdfc==7)
    return(HOD.M_cen_max);
  if(!(HOD.pdfc==6 || HOD.pdfc==8 || HOD.pdfc==9))
    return(HOD.M_max);

  m = HOD.M_min;
  n = N_cen(m);

  while(n>0.001)
    {
      m*=2;
      n = N_cen(m);
      if(m>HOD.M_max)return(HOD.M_max);
    }
  m = exp(zbrent(func_mhi,log(m/2),log(m),1.0E-5));
  return(m);
}

float func_mhi(float m)
{
  m=exp(m);
  return(N_cen(m)-0.001);
}

/* This is a copy of the above function that can be called from any routine.
 * (But this one integrates over dlogm).
 */
float func_halo_density(float m)
{
  float n1;

  m=exp(m);
  n1=dndM_interp(m);
  return(n1*m);
}


/* This is a copy of the above function that can be called from any routine.
 * (But this one integrates over dlogm
 */
float func_galaxy_density(float m)
{
  float n1,n2,m0;

  m=exp(m);
  n1=dndM_interp(m);
  n2=N_avg(m);
  //fprintf(stdout,"%e %e %e\n",m,n1,n2);
  
  return(n1*n2*m);
}

float func_mean_halo_mass(float m)
{
  float n1,n2,m0;

  m=exp(m);
  n1=dndM_interp(m);
  n2=N_avg(m);
  return(n1*n2*m*m);
}

/* This is the equation for zbrent to solve. What value
 * of M_min gives the correct galaxy density?
 * For HODs that have sharp central cutoffs, use M_min as the lower
 * limit of the integration. For soft cutoffs, first find M_low.
 */
float fixfunc1(float m)
{
  float n,mlo;

  HOD.M_min=m=pow(10.0,m);
  HOD.M_low=0;
  if(SOFT_CENTRAL_CUTOFF)
    mlo = (zbrent(func_mlow,log(HOD.M_min)-5*HOD.sigma_logM*2.3,log(HOD.M_min),1.0E-5));
  else
    mlo = log(HOD.M_min);
  if(exp(mlo)<1.0E7)mlo=log(1.0E+7);
  
  n=qromo(func_galaxy_density,mlo,log(HOD.M_max),midpnt);
  //    fprintf(stderr,"MMIN %e %e %e %e %e %e\n",m,exp(mlo),n,GALAXY_DENSITY,N_sat(2*m),N_cen(m)); 
  return(GALAXY_DENSITY-n);
}

/* This function is sent to zbrent to determine what M1 is based
 * on the number density. Both M_min and M_low have already been specified.
 */
float fixfunc2(float m)
{
  float n;
  HOD.M1=m=pow(10.0,m);  
  n=qromo(func_galaxy_density,log(HOD.M_low),log(HOD.M_max),midpnt);
  return(GALAXY_DENSITY-n);
}


/* This function is to be passed to qromo to integrate the number density
 * of satellite galaxies.
 */
float func_satfrac(float m)
{
  m=exp(m);
  return(N_sat(m)*dndM_interp(m)*m);
}

/* This function is to be passed to qromo to integrate the number density
 * of satellite galaxies.
 */
float func_satellite_density(float m)
{
  m=exp(m);
  return(N_sat(m)*dndM_interp(m)*m);
}

/* This function is to be passed to qromo to integrate the number density
 * of satellite galaxies.
 */
float func_central_density(float m)
{
  m=exp(m);
  return(N_cen(m)*dndM_interp(m)*m);
}

