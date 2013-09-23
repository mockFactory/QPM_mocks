#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "header.h"

#define NBINLOOKUP 10000

/* This function reads in a halo file and creates a mock galaxy distribution 
 * by populating the dark matter halos with galaxies with the currently specified
 * HOD functions.
 *
 * The format of the halo file needs to be in ASCII:
 * col 1 - halo ID (not used)
 * col 2 - number of particles in halo
 * cols 3-5 - x,y,z position of halo, in Mpc/h
 * col 6 - space velocity of halo (not used)
 * cols 7-9 - vx,vy,vz velocity of halo, in km/s
 *
 * The values of RESOLUTION, BOX_SIZE, OMEGA_M from the batch file are used
 * to convert particle number into mass, and to box-wrap galaxies.
 * If you want another file format, by all means edit.
 *
 * Output: mock galaxies are printed to [root].mock
 * in format: x,y,z [Mpc/h] vz,vy,vz [km/s]
 *
 * Output: HOD calculated from populated halos (for cross-checking with true HOD)
 * [root].binned_HOD in format: bin id, log10(mass), <N_gal>, N_halo
 * 
 * Satellite positions are chosen from a random sampling of the NFW profile
 * for the mass of the halo being populated. If CVIR_FAC is specified, then that value
 * will be used to adjust the NFW profile. Satellite velocities are chosen
 * from a Gaussin distribution with width = to virial dispersion of halo (plus
 * halo center-of-mass).
 *
 * NB-- If you are using the code to calculate M_min from a desired space density,
 * be sure that the linear matter power spectrum is the same as that used in the
 * simulation, or else the space densities won't match. [Mass function is used for this
 * purpose].
 *
 */
double NFW_central_velocity(double mass, double v[], double mag);
void calc_nbody_two_halo(float **gal, int *id, int ngal);
void calc_nbody_one_halo(float **gal, int *id, int ngal);

double *g21_rad, *g21_xi;
int g21_nrad;

/* External functions.
 */
void nbrsfind2(float smin,float smax,float rmax,int nmesh,float xpos,float ypos,float zpos,
               int *nbrmax,int *indx,float *rsqr,float *x,float *y,float *z,
               int *meshparts,int ***meshstart,int ip);
void meshlink2(int np1,int *nmesh,float smin,float smax,float rmax,float *x1,float *y1,float *z1,
	       int **meshparts,int ****meshstart,int meshfac);
void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);


void populate_simulation_hod()
{
  FILE *fp,*fpa[9],*fp2,*fpb[9],*fpc[9],*fps[9],*fpt,*fpblue,*fpcen;
  int i,j,k,n,imass,n1,j_start=0,i1,galcnt[1000],halocnt[1000],imag;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,r,fac,sigv;
  char aa[1000];
  float x1,xh[3],vh[3],vgf[3];
  long IDUM3 = -445;

  float **galarr;
  int *galid,id1=0,id2=0,j1;
  float dx,dy,dz,dr,drh,rv1,rv2,rmin,rmax;
  float **haloarr;
  int ngal,nsati[9],ALL_FILES=0,TRACK_GALAXIES=0,WARREN_MASS_CORRECTION=0,haloid;

  float *xt,*yt,*zt,*vxt,*vyt,*vzt;

  //static int nhalo, *used_halo, iuse_flag = 1;
  int ihalo=0;

  int SO_FILE = 1,
    hodflag = 0,
    JEANS_DISPERSION = 0;

  // fof stuff
  float *xxh, *yxh, *zxh, *vxh, *vyh, *vzh, *positions, *hmass;
  double *rvir, min_fofmass, p, rr, errfac,  ntot, massfac, volume, vfac, hubblez, bias1;
  static double galaxy_density_ref = -1;
  int skipflag, *numgrp, igrp, ii, nhalo, nfof, npbh;


  srand48(555);
  fprintf(stderr,"halopop> starting the population of the halo catalog...\n");

  // the fof halo file
  /*
  fprintf(stderr,"Opening [%s]\n",Files.FOFHaloFile);
  fp=fopen(Files.FOFHaloFile,"r");
  // get number of fofhalos
  fread(&nfof,sizeof(int),1,fp);
  */
  nfof = 0; // not reading in the FOF halos any more

  // the fake halo catalog
  fprintf(stderr,"Opening [%s]\n",Files.HaloFile);
  fp2=fopen(Files.HaloFile,"r");
  // get that number of halos
  fread(&npbh,sizeof(int),1,fp2);

  // get total number of halos and allocate memory
  nhalo = nfof + npbh;
  hmass = vector(1,nhalo);
  xxh = vector(1,nhalo);
  yxh = vector(1,nhalo);
  zxh = vector(1,nhalo);
  vxh = vector(1,nhalo);
  vyh = vector(1,nhalo);
  vzh = vector(1,nhalo);
  //positions = vector(0,3*nfof-1);

  // read in FOF halos
  /*
  fprintf(stderr,"Reading in the FOF halos...\n");
  vfac = BOX_SIZE*sqrt(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M))*100/(1+REDSHIFT);
  MASS_PER_PARTICLE = OMEGA_M*RHO_CRIT*pow(RESOLUTION,3.0);
  fprintf(stderr,"FOF Mass per particle: %e\n",MASS_PER_PARTICLE);
  fread(&hmass[1],sizeof(float),nfof,fp);
  for(i=1;i<=nfof;++i)
    hmass[i] = hmass[i]*MASS_PER_PARTICLE;

  fread(positions,sizeof(float),3*nfof,fp);
  j = -1;
  for(i=1;i<=nfof;++i)
    {
      xxh[i] =positions[++j]*BOX_SIZE;
      yxh[i] =positions[++j]*BOX_SIZE;
      zxh[i] =positions[++j]*BOX_SIZE;
    }
  fread(positions,sizeof(float),3*nfof,fp);
  j = -1;
  for(i=1;i<=nfof;++i)
    {
      vxh[i] =positions[++j]*vfac;
      vyh[i] =positions[++j]*vfac;
      vzh[i] =positions[++j]*vfac;
    }
  free_vector(positions,0,3*nfof-1);
  */

  if(ASCII_OUTPUT==2)
    {
      fp = fopen("testhalo.dat","w");
      for(i=1;i<=nfof;++i)
	fprintf(fp,"%d %e %e %e %e %e %e %e %e %e\n",
		i,hmass[i],0.0,0.0,xxh[i],yxh[i],zxh[i],vxh[i],vyh[i],vzh[i]);
      fclose(fp);
    }

  // read in the PB-halos
  fprintf(stderr,"Reading in the PB-halos...\n");
  fread(&hmass[1+nfof],sizeof(float),npbh,fp2);
  fread(&xxh[1+nfof],sizeof(float),npbh,fp2);
  fread(&yxh[1+nfof],sizeof(float),npbh,fp2);
  fread(&zxh[1+nfof],sizeof(float),npbh,fp2);
  fread(&vxh[1+nfof],sizeof(float),npbh,fp2);
  fread(&vyh[1+nfof],sizeof(float),npbh,fp2);
  fread(&vzh[1+nfof],sizeof(float),npbh,fp2);
  fclose(fp2);


  /*
  fp = fopen("testhalo.dat","w");
  for(i=1;i<=nhalo;++i)
    fprintf(fp,"%d %e %e %e %e %e %e %e %e %e\n",i,hmass[i],0.0,0.0,xxh[i],yxh[i],zxh[i],vxh[i],vyh[i],vzh[i]);
  fclose(fp);
  exit(0);
  */
  volume = BOX_SIZE*BOX_SIZE*BOX_SIZE;

  if(HOD.M_min>0)hodflag = 1;
  if(galaxy_density_ref<0)galaxy_density_ref = GALAXY_DENSITY;
  massfac = 1.5;

  ntot = 0;
  bias1 = 0;
  set_HOD_params();  
  HOD.M_low = 1.0E+10;
  for(i=1;i<=nhalo;++i) {
    ntot += N_cen(hmass[i]) + N_sat(hmass[i]);
    bias1 += (N_cen(hmass[i]) + N_sat(hmass[i]))*bias_interp(hmass[i],-1);
  }
  bias1 = bias1/ntot;
  GALAXY_DENSITY = ntot/volume;
  errfac = (GALAXY_DENSITY/galaxy_density_ref)-1.0;
  errfac = fabs(errfac);
  printf("%e %e %e %e\n",GALAXY_DENSITY,errfac,HOD.M_min,bias1);

  //if the full HOD is specified, skip the calibration step and go straight to populating halos
  if(hodflag)goto STARTPOP;

  // now constrain the number density
  //errfac = 1;
  while(errfac>0.01)
    {
      if(GALAXY_DENSITY>galaxy_density_ref)
	{
	  HOD.M_min *= massfac;
	  //HOD.M1 *= massfac;
	  //HOD.M_cut *= massfac;
	}	  
      if(GALAXY_DENSITY<=galaxy_density_ref)
	{
	  HOD.M_min /= massfac;
	  //HOD.M1 /= massfac;
	  //HOD.M_cut /= massfac;
	}	
      ntot = 0;
      bias1 = 0;
      for(i=1;i<=nhalo;++i)
	ntot += N_cen(hmass[i]) + N_sat(hmass[i]);
      for(i=1;i<=nhalo;++i)
	bias1 += (N_cen(hmass[i]) + N_sat(hmass[i]))*bias_interp(hmass[i],-1);
      GALAXY_DENSITY = ntot/volume;
      bias1 = bias1/ntot;
      errfac = (GALAXY_DENSITY/galaxy_density_ref)-1.0;
      errfac = fabs(errfac);
      printf("%e %e %e %e\n",GALAXY_DENSITY,errfac,HOD.M_min,bias1);
      massfac = (massfac-1)/2 + 1;
    }

 STARTPOP:
  printf("HOD PARAMS: %e %e %e %e %e\n",HOD.M_min, HOD.M1, HOD.M_cut, HOD.alpha, HOD.sigma_logM);

  // create mock output file
  sprintf(aa,"%s.mock",Task.root_filename);      
  fp2 = fopen(aa,"w");

  for(ii=1;ii<=nhalo;++ii)
    {
      haloid = ii;
      mass = hmass[ii];
      xh[0] = xxh[ii];
      xh[1] = yxh[ii];
      xh[2] = zxh[ii];
      vh[0] = vxh[ii];
      vh[1] = vyh[ii];
      vh[2] = vzh[ii];

      ncen=N_cen(mass);
      if(drand48()>ncen)
	goto SATELLITES;

      for(i=0;i<3;++i)
	{
	  if(xh[i]<0)xh[i]+=BOX_SIZE;
	  if(xh[i]>BOX_SIZE)xh[i]-=BOX_SIZE;
	  vg[i] = 0;
	}
      if(VBIAS_C>0)
	NFW_central_velocity(mass,vg,mag);
      fprintf(fp2,"%e %e %e %e %e %e %e 0\n",xh[0],xh[1],xh[2],vh[0]+vg[0],vh[1]+vg[1],vh[2]+vg[2],mass);

    SATELLITES:
      nsat = N_sat(mass);
      if(nsat>250)
	n1 = gasdev(&IDUM3)*sqrt(nsat) + nsat;
      else
	n1 = poisson_deviate(nsat);      
      
      for(i=1;i<=n1;++i)
	{
	  r = NFW_position(mass,xg);
	  NFW_velocity(mass,vg,mag);
	  for(k=0;k<3;++k)
	    {
	      xg[k]+=xh[k];
	      if(xg[k]<0)xg[k]+=BOX_SIZE;
	      if(xg[k]>BOX_SIZE)xg[k]-=BOX_SIZE;
	      vg[k]+=vh[k];
	    }	
	  fprintf(fp2,"%e %e %e %e %e %e %e 1\n",xg[0],xg[1],xg[2],vg[0],vg[1],vg[2],mass);
	}
      fflush(fp2);
    }
  fclose(fp2);

  return ;
  
}
/* Generate a random integer based on a Poisson distribution 
 * with mean given as input.
 */
int poisson_deviate(double nave)
{
  static int flag=0;
  double p,pp;
  int n;

  p=0;
  pp=1;

  while(p<pp)
    {
      if(nave<1)
	n=(int)(drand48()*20);
      else
	n=(int)(drand48()*30*nave);
      p=poisson_prob(n,nave);
      pp=drand48();
    }
  return(n);
}

/* Poisson probability of n given n_average
 */
double poisson_prob(int n, double nave)
{
  int i;
  double fac=1;

  if(n>0)
    for(i=1;i<=n;++i)
      fac*=nave/i;

  return((float)(fac*exp(-nave)));
}

/* Randomy generates a position away from the origin with 
 * a probability given by the NFW profile for a halo of the input
 * mass (and including the CVIR_FAC)
 */
double NFW_position(double mass, double x[])
{
  double r,pr,max_p,costheta,sintheta,phi1,signs,rvir,rs,cvir;
  
  cvir=halo_concentration(mass)*CVIR_FAC;
  rvir=pow(3*mass/(4*DELTA_HALO*PI*RHO_CRIT*OMEGA_M),1.0/3.0);
  rs=rvir/cvir;
  max_p=NFW_density(rs,rs,1.0)*rs*rs*4.0*PI;

  for(;;) {
    r=drand48()*rvir;
    pr=NFW_density(r,rs,1.0)*r*r*4.0*PI/max_p;
    
    if(drand48()<=pr)
      {
	costheta=2.*(drand48()-.5);
	sintheta=sqrt(1.-costheta*costheta);
	signs=2.*(drand48()-.5);
	costheta=signs*costheta/fabs(signs);
	phi1=2.0*PI*drand48();
	
	x[0]=r*sintheta*cos(phi1);
	x[1]=r*sintheta*sin(phi1);
	x[2]=r*costheta;
	return r;
      }
  }
}

/* This is the NFW density profile
 */
double NFW_density(double r, double rs, double ps)
{
  return(ps*rs/(r*(1+r/rs)*(1+r/rs)));
}

/* This sets the velocity to be isotropic Gaussian.
 */
double NFW_velocity(double mass, double v[], double mag)
{
  static long IDUM2=-455;
  static double fac = -1;
  double sigv,vbias=1;
  int i;

  if(fac<0)
      fac=sqrt(4.499E-48)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  sigv=fac*pow(mass,1.0/3.0)/sqrt(2.0);
  for(i=0;i<3;++i)
    v[i]=gasdev(&IDUM2)*sigv*VBIAS;
  return(0);
}

/* This sets the velocity to be isotropic Gaussian.
 */
double NFW_central_velocity(double mass, double v[], double mag)
{
  static long IDUM2=-455;
  static double fac = -1;
  double sigv,vbias=1;
  int i;

  if(fac<0)
      fac=sqrt(4.499E-48)*pow(4*DELTA_HALO*PI*OMEGA_M*RHO_CRIT/3,1.0/6.0)*3.09E19;
  sigv=fac*pow(mass,1.0/3.0)/sqrt(2.0);
  for(i=0;i<3;++i)
    v[i]=gasdev(&IDUM2)*sigv*VBIAS_C;
  return(0);
}
