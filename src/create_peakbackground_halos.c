#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "header.h"

#define NBINLOOKUP 10000
#define cnint(x) ((x-floor(x)) < 0.5 ? floor(x) : ceil(x))




// internal functions for the 2LPT stuff-- density estimation
double density_field_halo_mass(double delta);

double *g21_rad, *g21_xi;
int g21_nrad;
double CMASS_GALAXY_DENSITY,MASS_MIN,MASS_MAX;

/* External functions.
 */
void nbrsfind2(float smin,float smax,float rmax,int nmesh,float xpos,float ypos,float zpos,
               int *nbrmax,int *indx,float *rsqr,float *x,float *y,float *z,
               int *meshparts,int ***meshstart,int ip);
void meshlink2(int np1,int *nmesh,float smin,float smax,float rmax,float *x1,float *y1,float *z1,
	       int **meshparts,int ****meshstart,int meshfac);
void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);
int sample_particles(int np, float *x, float *y, float *z, float *vx, float *vy, float *vz, float *delta, float mass_limit);

void create_peakbackground_halos()
{
  FILE *fp,*fpa[9],*fp2,*fp4,*fpb[9],*fpc[9],*fps[9],*fpt,*fpblue, *fpsub, *fp_part,*fp3,*fp1;
  int i,j,k,n,imass,n1,j_start=0,i1,galcnt[1000],halocnt[1000], nfof;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,r,fac,sigv,mbar;
  char aa[1000];
  float x1,xh[3],vh[3],vgf[3],vv;
  long IDUM3 = -445;

  float **galarr;
  int *galid,id1=0,id2=0,j1;
  float dx,dy,dz,dr,drh,rv1,rv2,rmin,rmax, fred;
  float **haloarr, msample[11], mmax;
  int ngal,nsati[9],ALL_FILES=0,TRACK_GALAXIES=0,WARREN_MASS_CORRECTION=0,haloid, iflag;

  float *xt,*yt,*zt,*vxt,*vyt,*vzt;

  //static int nhalo, *used_halo, iuse_flag = 1;
  int ihalo=0, nhalo, npart, idum;
  double min_halo_mass, max_halo_mass, nexp, delta, t0, t1;

  int SO_FILE = 0,
    JEANS_DISPERSION = 0;
  
  // new for mcmc_lensing;
  double m1, m2, mfslope, rad, mstar, rproj, xs[20], vhalo, imag;
  long IDUM=-555;
  float xx[20], xd[3], xv[3], *subi, *submass, *temp, fraction, dbar;
  int ii;

  // Gadget stuff
  float *xp, *yp, *zp, *vxp, *vyp, *vzp, *mvect, mlim;
  int np, npcheck, *iii;
  char filename[1000], fname[1000];

  // nbrsfind
  float *rsqr, rcube, *density, delta_box;
  int *meshparts, ***meshstart,nmesh,meshfac,nbrmax,*indx;
  float *rsqr2, rmax2;
  int *meshparts2, ***meshstart2,nmesh2,meshfac2,*indx2;
  
  // fof stuff
  float *xxh, *yxh, *zxh, *vxh, *vyh, *vzh;
  float *fofmass, *rvir, min_fofmass, p, rr, *iotemp, *dvect, *mass2;
  int skipflag, skipflag2, *numgrp, igrp, *numgrp2, *imerge, *id, icnt=0, nhalo_tmp;

  //TPM input
  struct tpmhdr1 {
    int npart;
    int nsph;
    int nstar;
    float aa;
    float softlen; 
  } tpmhdr;
  float *positions, vfac;

  ASCII_OUTPUT = 0;
  
  MASS_MIN = log(3.0E11);
  MASS_MAX = log(3.0E14);
  min_fofmass = 32*OMEGA_M*RHO_CRIT*RESOLUTION*RESOLUTION*RESOLUTION;
  if(NO_FOF_HALOS) 
    { 
      MASS_MAX = log(HOD.M_max); 
      min_fofmass = HOD.M_max;
    }

  fp = fopen(Files.pmfile,"r");
  // get header
  fread(&idum,sizeof(int),1,fp);
  fread(&idum,sizeof(int),1,fp);
  fread(&tpmhdr,1,idum,fp);
  printf("particles in TPM file: %d\n",tpmhdr.npart);
  np = tpmhdr.npart;
  printf("done allocation. %d\n",np);
  //np = np*3;
  printf("done allocation. %d\n",np);
  //positions = vector(1,(long)(3*np));
  //positions = calloc((long)(np*3),(long)(sizeof(float)));
  // positions = malloc((long)(np*3*sizeof(float)));
  positions = malloc(((long)np)*3*(long)sizeof(float));
  printf("done allocation.\n");

  fflush(stdout);
  xp = vector(1,np);
  yp = vector(1,np);
  zp = vector(1,np);
  vxp = vector(1,np);
  vyp = vector(1,np);
  vzp = vector(1,np);
  density = vector(1,np);
  fread(positions,(long)sizeof(float),3*(long)np,fp);
  //fread(positions,(long)sizeof(float),(long)(3*np),fp);
  //fread(positions,(long)(sizeof(float)),(long)(3*np),fp);
  
  j = -1;
  for(i=1;i<=np;++i)
    {
      xp[i] =positions[++j]*BOX_SIZE;
      yp[i] =positions[++j]*BOX_SIZE;
      zp[i] =positions[++j]*BOX_SIZE;
    }
  
  fread(positions,(long)sizeof(float),3*(long)np,fp);
  //fread(positions,(long)sizeof(float),(long)(3*np),fp);
  j = -1;
  vfac = sqrt(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M))*100/(1+REDSHIFT);
  if(SIGV>0)
    fprintf(stdout,"Adding SIGV=%f random motions to all particles\n",SIGV);
  for(i=1;i<=np;++i)
    {
      vxp[i] =positions[++j]*BOX_SIZE*vfac + gasdev(&IDUM)*SIGV;
      vyp[i] =positions[++j]*BOX_SIZE*vfac + gasdev(&IDUM)*SIGV;
      vzp[i] =positions[++j]*BOX_SIZE*vfac + gasdev(&IDUM)*SIGV;
    }


  fread(&density[1],sizeof(float),np,fp);
  fclose(fp);
  //free_vector(positions,1,3*np);
  free(positions);
  for(i=1;i<=np;++i)
    density[i] = density[i] - 1;

  // now subsample the particles based on their density
  //np = sample_particles(np, xp,yp,zp,vxp,vyp,vzp,density,exp(MASS_MIN));

  // how many halos do we need?
  nexp = qromo(func_halo_density,MASS_MIN,MASS_MAX,midpnt)*BOX_SIZE*BOX_SIZE*BOX_SIZE;
  fraction = nexp/np;
  fprintf(stderr,"Fraction of particles to use: %e\n",fraction);

  // how many halos are we going to output
  x1 = qromo(func_halo_density,log(1.0E12),MASS_MAX,midpnt)*BOX_SIZE*BOX_SIZE*BOX_SIZE;
  nhalo = nexp*1.05; //extra 5% for the FOF halos-- just in case
  fprintf(stderr,"Number of halos (approximately) to output (logM>12): %e (%f)\n",x1,x1/np);
  
  // allocate memory for output
  xxh = vector(1,nhalo);
  yxh = vector(1,nhalo);
  zxh = vector(1,nhalo);
  vxh = vector(1,nhalo);
  vyh = vector(1,nhalo);
  vzh = vector(1,nhalo);
  rvir = vector(1,nhalo);
  mvect = vector(1,nhalo);
  dvect = vector(1,nhalo);
    

  // open the output file (incase we want ASCII output)
  sprintf(fname,"%s.PBHaloFile",Task.root_filename);
  fp3 = fopen(fname,"w");

  ii= 0 ;
  for(i=1;i<=np;++i)
    {
      if(drand48()>fraction)continue; 

      // monte carlo a halo mass
      mass = density_field_halo_mass(delta);

      // sample 10 times for this particle
      /*
      iflag = 0;
      for(j=1;j<=10;++j)
	{
	  msample[j] = density_field_halo_mass(delta);
	  if(msample[j]>1.0E12)iflag++;
	  if(msample[j]>mmax)mmax = msample[j];
	}
      if(iflag)
	{
	  j = (int)(drand48()*9.999) + 1;
	  while(msample[j]<1.0E12)
	    j = (int)(drand48()*9.999) + 1;
	  mmax = msample[j];
	}
      mass = mmax;
      // done with sampling
      */

      // take only halos where the number within the fof file is less than the Tinker MF
      if(mass<1.0E+12)continue;
      p = 1 - exp(-pow(4.5E+13/mass,0.75))*1.2;
      if(NO_FOF_HALOS)p = 1.0;
      if(mass<min_fofmass)p = 1.0;
      if(drand48()>p)continue;

      ii++;
      mvect[ii] = mass;
      xxh[ii] = xp[i];
      yxh[ii] = yp[i];
      zxh[ii] = zp[i];
      vxh[ii] = vxp[i];
      vyh[ii] = vyp[i];
      vzh[ii] = vzp[i];
      dvect[ii] = delta+1;
      rvir[ii] = pow(4*mvect[ii]/(3*OMEGA_M*RHO_CRIT*DELTA_HALO*PI),THIRD);

    }
  nhalo = ii;

  /* Free up the memory from the particles
   */
  fprintf(stderr,"PBhalo> free memory from particles\n");
  free_vector(xp,1,np);
  free_vector(yp,1,np);
  free_vector(zp,1,np);
  free_vector(vxp,1,np);
  free_vector(vyp,1,np);
  free_vector(vzp,1,np);
  free_vector(density,1,np);


  /* Now we need to remove the overlapping halos.
   */
  indx=malloc(nhalo*sizeof(int));
  rsqr=malloc(nhalo*sizeof(float));
  nmesh=0;
  rmax = 5;
  rcube = BOX_SIZE;
  fprintf(stderr,"PBhalos> creating mesh for halo cat\n");
  meshlink2(nhalo,&nmesh,0.0,rcube,rmax,&xxh[1],&yxh[1],&zxh[1],&meshparts,&meshstart,meshfac);
  
  /* make temp vector for sorting the halos
   */
  id = ivector(1,nhalo);
  imerge = ivector(1,nhalo);
  mass2 = vector(1,nhalo);
  for(i=1;i<=nhalo;++i)mass2[i] = -mvect[i];
  for(i=1;i<=nhalo;++i) { imerge[i] = 0; id[i] = i; }

  /* Sorting the halos
   */
  fprintf(stderr,"PBhalos> Sorting the halos\n");
  sort2(nhalo,mass2,id);
  fprintf(stderr,"PBhalos> Done sorting the halos\n");

  for(ii=1;ii<=nhalo;++ii)
    {
      i = id[ii];
      nbrmax=nhalo;
      rmax = rvir[i];
      if(imerge[i])continue;
      //printf("%d %e %e %e\n",i,xxh[i],yxh[i],zxh[i]);fflush(stdout);
      nbrsfind2(0.0,rcube,rmax,nmesh,xxh[i],yxh[i],zxh[i],&nbrmax,indx,rsqr,&xxh[1],&yxh[1],&zxh[1],
		meshparts,meshstart,i);
      for(j1=0;j1<nbrmax;++j1)
	indx[j1]++;
      for(j1=0;j1<nbrmax;++j1)
	{
	  /* Let's check our indices
	   */
	  if(indx[j1]==i) {
	    continue; }
	  icnt++;
	  j = indx[j1];
	  if(mvect[i]>mvect[j])imerge[j]=1;
	  else imerge[i]=1;
	  //printf("%d %d %d %f %f %f %f %f %f\n",icnt,i,j,xxh[i],yxh[i],zxh[i],xxh[j],yxh[j],zxh[j]);
	  fflush(stdout);
	}
    }
  fprintf(stderr,"PBHalos> Removing %d overlapping halos\n",icnt);

  j = 0;
  for(i=1;i<=nhalo;++i)
    if(!imerge[i])++j;
  nhalo_tmp = j;
  for(j=0,i=1;i<=nhalo;++i)
    if(!imerge[i])rvir[++j] = mvect[i];
  for(j=0,i=1;i<=nhalo_tmp;++i)
    mvect[i] = rvir[i];
  for(j=0,i=1;i<=nhalo;++i)
    if(!imerge[i])rvir[++j] = xxh[i];
  for(j=0,i=1;i<=nhalo_tmp;++i)
    xxh[i] = rvir[i];
  for(j=0,i=1;i<=nhalo;++i)
    if(!imerge[i])rvir[++j] = yxh[i];
  for(j=0,i=1;i<=nhalo_tmp;++i)
    yxh[i] = rvir[i];
  for(j=0,i=1;i<=nhalo;++i)
    if(!imerge[i])rvir[++j] = zxh[i];
  for(j=0,i=1;i<=nhalo_tmp;++i)
    zxh[i] = rvir[i];
  for(j=0,i=1;i<=nhalo;++i)
    if(!imerge[i])rvir[++j] = vxh[i];
  for(j=0,i=1;i<=nhalo_tmp;++i)
    vxh[i] = rvir[i];
  for(j=0,i=1;i<=nhalo;++i)
    if(!imerge[i])rvir[++j] = vyh[i];
  for(j=0,i=1;i<=nhalo_tmp;++i)
    vyh[i] = rvir[i];
  for(j=0,i=1;i<=nhalo;++i)
    if(!imerge[i])rvir[++j] = vzh[i];
  for(j=0,i=1;i<=nhalo_tmp;++i)
    vzh[i] = rvir[i];
  for(j=0,i=1;i<=nhalo;++i)
    if(!imerge[i])rvir[++j] = dvect[i];
  for(j=0,i=1;i<=nhalo_tmp;++i)
    dvect[i] = rvir[i];
  nhalo = nhalo_tmp;

  if(NO_FOF_HALOS) goto SKIP_FOF_HALOS;

  // the fof halo file
  fprintf(stderr,"Opening [%s]\n",Files.FOFHaloFile);
  fp=fopen(Files.FOFHaloFile,"r");
  // get number of fofhalos
  fread(&nfof,sizeof(int),1,fp);
  positions = vector(0,3*nfof-1);

  fprintf(stderr,"FOF halos: %d (+%d=%d)\n",nfof, nhalo, nhalo+nfof);
  // read in FOF halos
  fprintf(stderr,"Reading in the FOF halos...\n");
  vfac = BOX_SIZE*sqrt(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M))*100/(1+REDSHIFT);
  MASS_PER_PARTICLE = OMEGA_M*RHO_CRIT*pow(RESOLUTION,3.0);
  fprintf(stderr,"FOF Mass per particle: %e\n",MASS_PER_PARTICLE);
  fread(positions,sizeof(float),nfof,fp);
  j = -1;
  for(i=nhalo;i<=nhalo+nfof;++i)
    mvect[i] = positions[++j]*MASS_PER_PARTICLE;

  fread(positions,sizeof(float),3*nfof,fp);
  j = -1;
  for(i=nhalo+1;i<=nhalo+nfof;++i)
    {
      xxh[i] =positions[++j]*BOX_SIZE;
      yxh[i] =positions[++j]*BOX_SIZE;
      zxh[i] =positions[++j]*BOX_SIZE;
    }
  fread(positions,sizeof(float),3*nfof,fp);
  j = -1;
  for(i=nhalo;i<=nhalo+nfof;++i)
    {
      vxh[i] =positions[++j]*vfac;
      vyh[i] =positions[++j]*vfac;
      vzh[i] =positions[++j]*vfac;
      dvect[i] = -1; // set the density of the FOF halos
    }
  free_vector(positions,0,3*nfof-1);

  nhalo += nfof; 

 SKIP_FOF_HALOS:
  if(ASCII_OUTPUT)
    {
      for(i=1;i<=nhalo;++i)
	{
	  fprintf(fp3,"%d %e %e %e %e %e %e %e %e %e\n",ii,mvect[i],-1.0,dvect[i],
		  xxh[i],yxh[i],zxh[i],vxh[i],vyh[i],vzh[i]);
	  fflush(fp3);
	}
    }


  if(!ASCII_OUTPUT)
    {
      fwrite(&nhalo,sizeof(int),1,fp3);
      fwrite(mvect,sizeof(float),nhalo,fp3);
      fwrite(xxh,sizeof(float),nhalo,fp3);
      fwrite(yxh,sizeof(float),nhalo,fp3);
      fwrite(zxh,sizeof(float),nhalo,fp3);
      fwrite(vxh,sizeof(float),nhalo,fp3);
      fwrite(vyh,sizeof(float),nhalo,fp3);
      fwrite(vzh,sizeof(float),nhalo,fp3);
      fwrite(dvect,sizeof(float),nhalo,fp3);
    }
  fclose(fp3);
  fprintf(stderr,"PBhalos> done with creating PBHaloFile.\n");
  
}

/* 
 * n(M|delta) = n(M)*(1+delta*b(M))
 */
double density_field_halo_mass(double delta)
{
  double m,p,pmax,f,logm;
  static int iter=0;
  pmax = dndM_interp(exp(MASS_MIN))*exp(MASS_MIN)*
    (1+delta*bias_interp(exp(MASS_MIN),-1));
  p = -1;
  while (drand48()>p) {
    m = exp(drand48()*(MASS_MAX-MASS_MIN) + MASS_MIN);
    logm = log10(m);
    f = 1;
    if(logm>13)
      f = 1.5*pow((logm-13)/2,2.0) + 1;
    p = dndM_interp(m)*m*(1+delta*bias_interp(m,-1)*f)/pmax;
  }
  return m;

}



