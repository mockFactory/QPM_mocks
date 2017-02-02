#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <stdint.h>

#include "header.h"

#define NBINLOOKUP 10000
#define cnint(x) ((x-floor(x)) < 0.5 ? floor(x) : ceil(x))
#define SIG 0.1
#define SIG2INV (1/(SIG*SIG))

#define xPrhoMhalo(x,y) exp(-0.5*(x-y)*(x-y)*SIG2INV)


/* Internal functions.
 */
float PrhoMhalo(const float logrho, const float logMhalo) ;


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
int ftread(char *ptr,unsigned size, unsigned nitems,FILE *stream);

void create_lognormal_halos()
{
  FILE *fp,*fpa[9],*fp2,*fp4,*fpb[9],*fpc[9],*fps[9],*fpt,*fpblue, *fpsub, *fp_part,*fp3,*fp1;
  int j_start=0,i1,galcnt[1000],halocnt[1000], nfof;
  double mass,xg[3],vg[3],nsat,nc[10],ncen,mlo,mag,err1,err2,r,fac,sigv,mbar;
  char aa[1000];
  float x1,xh[3],vh[3],vgf[3],vv,x1_;
  long IDUM3 = -445;

  long i,j,k,n,n1,imass,ngal,j1,jstep,jstep_max = 6;

  float **galarr;
  int *galid,id1=0,id2=0;
  float dx,dy,dz,dr,drh,rv1,rv2,rmin,rmax, fred;
  float **haloarr;
  int nsati[9],ALL_FILES=0,TRACK_GALAXIES=0,WARREN_MASS_CORRECTION=0,haloid;

  float *xt,*yt,*zt,*vxt,*vyt,*vzt;
  float lm, logMhalo, mu, sig;

  //static int nhalo, *used_halo, iuse_flag = 1;
  long ihalo=0, nhalo, npart, idum, nmax, nadd;
  double min_halo_mass, max_halo_mass, nexp, delta, t0, t1;

  int SO_FILE = 0,
    JEANS_DISPERSION = 0;
  
  // new for mcmc_lensing;
  double m1, m2, mfslope, rad, mstar, rproj, xs[20], vhalo, imag;
  long IDUM=-555;
  float xx[20], xd[3], xv[3], *subi, *submass, *temp, fraction, dbar;
  long ii, ntotread=0, np1, ncounter=0;

  // Gadget stuff
  float *xp, *yp, *zp, *vxp, *vyp, *vzp, *mvect, mlim;
  int npcheck, *iii;
  long np;
  char filename[1000], fname[1000];

  // nbrsfind
  float *rsqr, rcube, *density, delta_box;
  int *meshparts, ***meshstart,nmesh,meshfac,nbrmax,*indx;
  float *rsqr2, rmax2;
  int *meshparts2, ***meshstart2,nmesh2,meshfac2,*indx2;
  double xran;
  
  // fof stuff
  float *xxh, *yxh, *zxh, *vxh, *vyh, *vzh;
  float *fofmass, *rvir, min_fofmass, p, rr, *iotemp, *dvect, *mass2;
  int skipflag, skipflag2, *numgrp, igrp, *numgrp2, *imerge, *id, icnt=0, nhalo_tmp;

  // cpm stuff
  float *fdat, znow;
  int *idat;
  

  //TPM input
  struct tpmhdr1 {
    int npart;
    int nsph;
    int nstar;
    float aa;
    float softlen; 
  } tpmhdr;
  float *positions, vfac;

  float logprho, dlnMbin, prob, *pnorm, logM, lnRhoMin, lnRhoMax, *rhohist, logrho, *mux, *mux1,m0;
  int NmassBin=256, jbin, NrhoBin=2048, ibin;
  double avg, ntar, ntot;

  long ichunk, nchunk, nchunk_reads, nread;

  // openmp stuff
  int irank=0,nrank=1;
  float *logMx, **pnormx;
  int *iix;
  unsigned int iseed = -555;
  struct drand48_data drand_buf;


  lnRhoMin = log(1.0E-3);
  lnRhoMax = log(100);
  rhohist = vector(0,NrhoBin-1);
  for(i=0;i<NrhoBin;++i)rhohist[i] = 0;

  ASCII_OUTPUT = 0;

  
  if(MASS_MIN<=0)
    MASS_MIN = log(1.0E12);
  else
    MASS_MIN = log(MASS_MIN);
  MASS_MAX = log(HOD.M_max); 
  printf("LNHalos> mass min, mass max: %e %e\n",exp(MASS_MIN), exp(MASS_MAX));

  NmassBin = (int)((MASS_MAX-MASS_MIN)/0.05+1);
  dlnMbin = (MASS_MAX-MASS_MIN)/NmassBin; 
  pnorm = vector(0,NmassBin-1);


  // if we're reading in from the CPM code, check here.
  if(Files.NumFiles>=0) goto TPM_INPUT;

  sprintf(fname,"%s",Files.pmfile);
  fp = fopen(fname,"r");
  fprintf(stderr,"%s\n",fname);

  idat=(int *)calloc(5,sizeof(int));
  fdat=(float *)calloc(9,sizeof(float));

  ftread(idat,sizeof(int),5,fp);
  ftread(fdat,sizeof(float),9,fp);
  np=idat[1];
   
  printf("LNHalos: total number of partciles: %ld\n",np);

  xp = vector(1,np);
  yp = vector(1,np);
  zp = vector(1,np);
  vxp = vector(1,np);
  vyp = vector(1,np);
  vzp = vector(1,np);
  density = vector(1,np);
  printf("done allocation.\n");
  fflush(stdout);

  ftread(&znow,sizeof(float),1,fp) ;
  ftread(&xp[1],sizeof(float),np,fp) ;
  ftread(&yp[1],sizeof(float),np,fp) ;
  ftread(&zp[1],sizeof(float),np,fp) ;
  ftread(&vxp[1],sizeof(float),np,fp) ;
  ftread(&vyp[1],sizeof(float),np,fp) ;
  ftread(&vzp[1],sizeof(float),np,fp) ;
  fclose(fp);

  // add in the artificial dispersion to the particle velocities
  for(i=1;i<=np;++i)
    {
      vxp[i] += gasdev(&IDUM)*SIGV;
      vyp[i] += gasdev(&IDUM)*SIGV;
      vzp[i] += gasdev(&IDUM)*SIGV;
    }

  // now get the density file
  sprintf(fname,"%s.den",Files.pmfile);
  fp = fopen(fname,"r");
  fprintf(stderr,"%s\n",fname);
  ftread(&density[1],sizeof(float),np,fp) ;
  fclose(fp);
  fprintf(stderr,"here\n");

  goto SKIP_TPM;

 TPM_INPUT:
  // get the total number of particles in all files
  np = 0;
  for(i=0;i<Files.NumFiles;++i)
    {
      if(Files.NumFiles>1)
	sprintf(fname,"%s.%02d",Files.pmfile,i);
      else
	sprintf(fname,"%s",Files.pmfile);
      fp = fopen(fname,"r");
      fprintf(stderr,"%s\n",fname);
      // get header
      fread(&idum,sizeof(int),1,fp);
      fread(&idum,sizeof(int),1,fp);
      fread(&tpmhdr,1,idum,fp);
      np += tpmhdr.npart;
      printf("particles in TPM file (%d): %d %ld %d %d\n",i,tpmhdr.npart,np,idum,sizeof(tpmhdr));
      printf("aexp: %f\n",tpmhdr.aa);
      fflush(stdout);
      fclose(fp);
    }
  printf("LNHalos: total number of partciles: %ld\n",np);

  xp = vector(1,np);
  yp = vector(1,np);
  zp = vector(1,np);
  vxp = vector(1,np);
  vyp = vector(1,np);
  vzp = vector(1,np);
  density = vector(1,np);
  printf("done allocation.\n");
  fflush(stdout);

  // read in the data in 10 discrete chunks

  ncounter = 0;
  for(i1=0;i1<Files.NumFiles;++i1)
    {
      if(Files.NumFiles>1)
	sprintf(fname,"%s.%02d",Files.pmfile,i1);
      else
	sprintf(fname,"%s",Files.pmfile);
      fp = fopen(fname,"r");
      fprintf(stderr,"%s\n",fname);
      // get header
      fread(&idum,sizeof(int),1,fp);
      fread(&idum,sizeof(int),1,fp);
      //fread(&tpmhdr,1,idum,fp);
      fread(&tpmhdr,sizeof(tpmhdr),1,fp);


      np1 = tpmhdr.npart;
      nchunk = tpmhdr.npart/10;
      positions = malloc(((long)nchunk)*3*(long)sizeof(float));
      nchunk_reads = tpmhdr.npart/nchunk + 1;
      if(tpmhdr.npart%nchunk==0)nchunk_reads--;
      printf("done allocation of temp array. nchunk= %d %d %d %ld\n",
	     nchunk,nchunk_reads,tpmhdr.npart,np);
      fflush(stdout);
      
      ntotread = 0;
      ii = ncounter + 1;
      for(ichunk=1;ichunk<=nchunk_reads;++ichunk)
	{
	  nread = nchunk;
	  if(ntotread + nread>np1) nread = np1-ntotread;
	  fprintf(stderr,"%d %ld %ld %ld %d %ld\n",ichunk,ii,nread,ntotread,tpmhdr.npart,np);
	  fread(positions,(long)sizeof(float),3*(long)nread,fp);
	  
	  j = -1;
	  for(i=1;i<=nread;++i)
	    {
	      xp[ii] =positions[++j]*BOX_SIZE;
	      yp[ii] =positions[++j]*BOX_SIZE;
	      zp[ii] =positions[++j]*BOX_SIZE;
	      //if(drand48()<0.0001)printf("TEST %f %f %f %d %d\n",xp[ii],yp[ii],zp[ii],ii,np);
	      ii++;
	    }
	  ntotread += nread;
	}
      fprintf(stderr,"Read %ld positions (%ld total out of %ld)\n",ntotread,ii-1,np);
      

      // now do the same thing for the velocities
      ii = ncounter + 1;
      ntotread = 0;
      vfac = sqrt(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M))*100/(1+REDSHIFT);
      if(SIGV>0)
	fprintf(stdout,"Adding SIGV=%f random motions to all particles\n",SIGV);
      for(ichunk=1;ichunk<=nchunk_reads;++ichunk)
	{
	  fprintf(stderr,"velread: chunk %d/%d\n",ichunk,nchunk_reads);
	  nread = nchunk;
	  if(ntotread + nread>np1) nread = np1-ntotread;
	  fread(positions,(long)sizeof(float),3*(long)nread,fp);
	  
	  j = -1;
	  for(i=1;i<=nread;++i)
	    {
	      vxp[ii] =positions[++j]*BOX_SIZE*vfac + gasdev(&IDUM)*SIGV;
	      vyp[ii] =positions[++j]*BOX_SIZE*vfac + gasdev(&IDUM)*SIGV;
	      vzp[ii] =positions[++j]*BOX_SIZE*vfac + gasdev(&IDUM)*SIGV;
	      ii++;
	    }
	  ntotread += nread;
	}
      fprintf(stderr,"Read %d velocities\n",ntotread);
      free(positions);

      positions = malloc(((long)np1)*(long)sizeof(float));
      fread(positions,sizeof(float),np1,fp);
      for(ii=ncounter+1, i=0 ;i<np1;++ii, ++i)
	density[ii] = positions[i];
      free(positions);

      fclose(fp);
      i = ncounter + 1;
      //fprintf(stdout,"%f %f %f %f %f %f %f %ld\n",xp[i],yp[i],zp[i],vxp[i],vyp[i],vzp[i],density[i],i); 
      ncounter += np1;
    }

 SKIP_TPM:

  if(SUBFRAC>0)
    {
      if(ARGC>3)
	fp2 = fopen(ARGV[3],"w");
      else 
	fp2 = fopen("pm_ascii.dat","w");
      n = 0;
      fprintf(stderr,"Opening file for random sample. Subfrac= %f\n",SUBFRAC);

      if(SUBFRAC>1)
	{
	  for(i=1;i<=np;++i)
	    { 
	      if(xp[i]<100 && yp[i]<100 && zp[i]<10) {
		n++;
		fprintf(fp2,"%f %f %f %f %f %f %f %ld\n",xp[i],yp[i],zp[i],vxp[i],vyp[i],vzp[i],density[i],n); }
	    }
	  exit(0);
	}

      for(i=1;i<=np;++i)
	if(drand48()<SUBFRAC) { n++;
	  fprintf(fp2,"%f %f %f %f %f %f %f %d\n",xp[i],yp[i],zp[i],vxp[i],vyp[i],vzp[i],density[i],n); }
      fprintf(stdout,"selected %d particles\n",n);
      exit(0);
    }

  //REDSHIFT = 0; // do we keep this? TODO
  //RESET_COSMOLOGY++;

  // how many halos do we need?
  nexp = qromo(func_halo_density,MASS_MIN,MASS_MAX,midpnt)*BOX_SIZE*BOX_SIZE*BOX_SIZE;
  fraction = nexp/np;
  fprintf(stderr,"Fraction of particles to use: %e\n",fraction);

  // how many halos are we going to output
  nhalo = nexp*1.01; 
  
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
    
  // loop for histogramming the rhos over all sampled particles
  ii = 0;
  for (j=1; j<=np; j++) {
    if (density[j]>=0) {
      logprho = log(density[j]);
      if (logprho>lnRhoMin && logprho<lnRhoMax) {
	ibin=(int)( NrhoBin*(logprho-lnRhoMin)/(lnRhoMax-lnRhoMin) );
	rhohist[ibin] += 1;
	ii++;
      }
    }
  }
  fprintf(stderr,"Done getting histogram of densities. %d %d\n",np,ii);

  // notes: when i put this inside the parallel region, it never terminates. when i put this outside the parallel region, it doesn't work. WTF. idea: try putting it inside but only executing for itask==0?
    mux = vector(0,NmassBin-1);
    for(i=0;i<NmassBin;++i)
      {
	// my function, z=1.5
	lm = exp(MASS_MIN+dlnMbin*(i+0.5));
	if(REDSHIFT>1)
	  {
	    m0 = pow(10,13.25);
	    //mux[i] = 1 + 0.13*log10(lm/m0) + pow(lm/m0,0.35)/(1+pow(lm/m0,-0.7)) - 0.5;
	    mux[i] = 1 + 0.115*log10(lm/m0) + pow(lm/m0,0.35)/(1+pow(lm/m0,-0.7)) - 0.5;
	  }
	else
	  {
	    m0 = pow(10,13.36);
	    mux[i] = 1 + 0.08*log10(lm/m0) + pow(lm/m0,0.35)/(1+pow(lm/m0,-0.7)) - 0.5;
	  }
	printf("MUX %d %e %e %e\n",i,lm,log10(lm),mux[i]);
      }

  // get the normalization for each bin

  // First compute the "average" number of particles selected per mass bin,
  // then the normalization is the ratio to what we want.
  for (jbin=0; jbin<NmassBin; jbin++) {
    logM = MASS_MIN+dlnMbin*(jbin+0.5);
    //ntar = exp(nofm->val(logM))*(MASS_MAX-MASS_MIN)/NmassBin * cc->vol; //density?
    ntar = qromo(func_halo_density, logM-0.5*dlnMbin, logM+0.5*dlnMbin, midpnt)*
      BOX_SIZE*BOX_SIZE*BOX_SIZE;
    ntot += ntar;
    avg  = 1e-30;
    for ( ibin=0; ibin<NrhoBin; ibin++) {
      logrho= lnRhoMin+(ibin+0.5)*(lnRhoMax-lnRhoMin)/NrhoBin;
      avg += xPrhoMhalo(logrho,mux[jbin])*rhohist[ibin];
    }
    pnorm[jbin]=ntar/avg;
    printf("NORM: %d %e %e %e %e %e %e\n",jbin,exp(logM),pnorm[jbin],
	   xPrhoMhalo(log(0.5),mux[jbin]),xPrhoMhalo(log(1),mux[jbin]),
    	   xPrhoMhalo(log(2),mux[jbin]),xPrhoMhalo(log(4),mux[jbin]));
    if (pnorm[jbin]>1.0) {
      fprintf(stderr,"Error in normalization for bin %d %e\n",jbin,pnorm[jbin]);
      fprintf(stderr,"%e %e %e\n",ntar, avg, exp(logM));
      pnorm[jbin]=1;
      //exit(1);
    }
  }
 //exit(0);
  fprintf(stderr,"Done getting histogram normalization. %e %e\n",nexp,ntot);

  /*----- Martin's code for converting the sampled particles into
   * masses via their densities.
   */
  fprintf(stderr,"%d\n",nrank);
  //logMx = vector(0,nrank-1); // Tinker is confused why this is nrank and not NmassBin
  logMx = vector(0,NmassBin-1); // JLT
  iix = ivector(0,nrank-1);
  for(i=0;i<nrank;++i)
    iix[i] = -1;

  // initialize the masses
  for(i=1;i<=nhalo;++i)
    mvect[i] = -1;

  system("date");
  t0 = second();
  //np = 10000000;
  //OUTPUT = 0;
  muh(NmassBin);

  for(i=0;i<NmassBin;++i)
    {
      lm = MASS_MIN+dlnMbin*(i+0.5);
      logMx[i] = lm;
      //mux[i] = 0.4863-0.0007793*lm+0.07161*lm*lm; // martin's function
      // my function, z=1.5
      //lm = exp(logMx[i]);
      //m0 = pow(10,13.25);
      //mux[i] = 1 + 0.13*log10(lm/m0) + pow(lm/m0,0.35)/(1+pow(lm/m0,-0.7)) - 0.5;
      //printf("%d %e %e %e\n",i,lm,logMx[i]/log(10),mux[i]);
    }
  fflush(stdout);
  muh(0);

  nmax = 0;
  sig = 1/0.1;
#pragma omp parallel private(ii, logM, logprho, jbin, prob, irank, nrank, j, lm, mu,x1, iseed, sig, xran, drand_buf,i,mux, mux1) \
  shared(np, xp, yp, zp, vxp, vyp, vzp, mvect, dvect, rvir, density, xxh, yxh, zxh, vxh, vyh, vzh,nhalo,pnorm,pnormx, logMx)
  {
    muh(1);
    nrank = omp_get_num_threads();
    muh(2);
    fprintf(stderr, "rank: %d %d\n",irank=omp_get_thread_num(),nrank);
    iseed = iseed - irank;
    srand48_r (iseed, &drand_buf);

    mux1 = vector(0,NmassBin-1);
    for(i=0;i<NmassBin;++i)
      {
	lm = exp(MASS_MIN+dlnMbin*(i+0.5));
	if(REDSHIFT>1)
	  {
	    m0 = pow(10,13.25);
	    //mux1[i] = 1 + 0.13*log10(lm/m0) + pow(lm/m0,0.35)/(1+pow(lm/m0,-0.7)) - 0.5;
	    mux1[i] = 1 + 0.115*log10(lm/m0) + pow(lm/m0,0.35)/(1+pow(lm/m0,-0.7)) - 0.5;
	  }
	else
	  {
	    m0 = pow(10,13.36);
	    mux1[i] = 1 + 0.08*log10(lm/m0) + pow(lm/m0,0.35)/(1+pow(lm/m0,-0.7)) - 0.5;
	  }
	//printf("BOO1 %d %d %e %e\n",irank,i,lm,mux1[i]);
      }
    #pragma omp barrier

    ii = irank+1;

    for(j=ii;j<=np;j+=nrank) {
      if( (j%1000000==0)&&(OUTPUT) ){ 
	fprintf(stderr,"%ld %e %e %d %f\n",j,nexp,j*1./np,jbin,ii*1./(j)); }
    if(density[j]<0){ printf("ZERO %ld %f\n",j,density[j]); fflush(stdout); }
    if (density[j]>0) {
      logprho = log(density[j]);

      for ( jbin=NmassBin-1; jbin>=0; jbin--) {
	//logM = MASS_MIN+dlnMbin*(jbin+0.5);
	//prob = PrhoMhalo(logprho,logM) * pnorm[jbin];

	//mu=mux1[jbin];// 0.4863-0.0007793*lm+0.07161*lm*lm; //martin's version
	//sig = 0.1;
	//x1=(logprho-mu)/sig;
	//prob=exp(-0.5*x1*x1)*pnorm[jbin];

	prob = xPrhoMhalo(mux1[jbin],logprho)*pnorm[jbin];
	/*
	if(!irank) fprintf(stdout,"%d %e %e %e %e %e %e %e %e %e\n",
			   jbin,mux1[jbin],mu,logprho,logMx[jbin]/log(10),prob,pnorm[jbin],x1,exp(-0.5*x1*x1),xPrhoMhalo(mu,logprho));
	if(!irank && !jbin)exit(0);
	*/
	if (prob>1.0) {
	  fprintf(stderr,"Out-of-bounds in create. %e %e %e %e\n",prob,mux1[jbin],logprho,pnorm[jbin]);
	  exit(1);
	}
	drand48_r(&drand_buf, &xran);
	if (xran<prob) {       // Select this as a halo.
	  xxh[ii] = xp[j];
	  yxh[ii] = yp[j];
	  zxh[ii] = zp[j];
	  
	  vxh[ii] = vxp[j];
	  vyh[ii] = vyp[j];
	  vzh[ii] = vzp[j];
	  
	  drand48_r (&drand_buf, &xran);
	  mvect[ii] = exp(logMx[jbin]+(xran-0.5)*dlnMbin);
	  dvect[ii] = density[j];
	  rvir[ii] = pow(4*mvect[ii]/(3*OMEGA_M*RHO_CRIT*DELTA_HALO*PI),THIRD);
	  drand48_r(&drand_buf, &xran);
	  
	  if(xran<-0.0001)
	    printf("BOO %d %d %d %e %e %e %d\n",irank, j, jbin, logprho, logM, prob, ii);
	  ii+=nrank;

	  //	  fprintf(stdout,"HALO %d %d %e %e %e %e %e %e %e %e %e\n",j,ii,mvect[ii],rvir[ii],dvect[ii],
	  //	  xxh[ii],yxh[ii],zxh[ii],vxh[ii],vyh[ii],vzh[ii]);
	  //fflush(stdout);
	  break;
	}
      }
      //if(j==3)      exit(0);
    }
    } 

    ii -= nrank;
    fprintf(stdout,"nhalo(irank=%d)= %d\n",irank,ii); 
    fflush(stdout);
    nhalo = nmax;
    
#pragma omp critical 
  {
    if(ii<nhalo) nhalo = ii;
    if(ii>nmax) nmax = ii;
  }
  #pragma omp barrier
  }

  t1 = second();
  printf("TIME %.3f\n",timediff(t0,t1));
  system("date");

  // remove the array entries with no halos assigned, get the total number of halos
  muh(nhalo);
  muh(nmax);
  nadd = nhalo;
  for(j=nhalo+1;j<=nmax;++j)
    {
      if(mvect[j]>0) 
	{
	  nadd++;
	  mvect[nadd] = mvect[j];
	  xxh[nadd] = xp[j];
	  yxh[nadd] = yp[j];
	  zxh[nadd] = zp[j];
	  
	  vxh[nadd] = vxp[j];
	  vyh[nadd] = vyp[j];
	  vzh[nadd] = vzp[j];
	  
	  dvect[nadd] = density[j];
	  rvir[nadd] = rvir[j];
	}
    }

  nhalo = nadd;

  /*-- end martin's code
   */
  fprintf(stderr,"Number of halos to be outputted (pre-overlaps): %d\n",nhalo);

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
  nmesh=0;
  rmax = 3;
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

#pragma omp parallel private(ii, i, nbrmax, rmax, j1, icnt, j, indx, rsqr, irank, nrank)
  {
    // allocate these internally
    indx=malloc(nhalo*sizeof(int));
    rsqr=malloc(nhalo*sizeof(float));

   irank=omp_get_thread_num() + 1;
   nrank = omp_get_num_threads();

   for(ii=irank;ii<=nhalo;ii+=nrank)
     {
       i = id[ii];
       nbrmax=nhalo;
       rmax = rvir[i];
       if(imerge[i])continue;
       //printf("%d %e %e %e\n",i,xxh[i],yxh[i],zxh[i]);fflush(stdout);
       nbrsfind2(0.0,rcube,rmax,nmesh,xxh[i],yxh[i],zxh[i],&nbrmax,indx,rsqr,
		 &xxh[1],&yxh[1],&zxh[1],meshparts,meshstart,i);
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
	   //fflush(stdout);
	 }
     }
  }
   #pragma omp barrier

  for(icnt=0,i=1;i<=nhalo;++i)
    if(imerge[i])icnt++;
  fprintf(stderr,"LNHalos> Removing %d overlapping halos\n",icnt);

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
  muh(nhalo);

  // open the output file (incase we want ASCII output)
  sprintf(fname,"%s.HaloFile",Task.root_filename);
  fp3 = fopen(fname,"w");



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
  fprintf(stderr,"PBhalos> done with creating LNHaloFile.\n");
  exit(0);
}



// Returns P(rho|Mhalo), which is a Gaussian/lognormal.
float PrhoMhalo(const float logrho, const float logMhalo) 
{
  float lm=logMhalo-28.0;     // Conveniently centered.
  float mu=0.4863-0.0007793*lm+0.07161*lm*lm; //martin's version
  //float lm=logMhalo-27.0;     // Conveniently centered. (JLT)
  //float mu=0.3+0.07161*lm*lm/1.7+0.002*lm*lm*lm; //my modification
  float sig=0.1; //martin
  //float sig=0.7;

  // try2 (meant to go with sig=0.1) results: 22% selected 
  //float lm=logMhalo-25.0;     // Conveniently centered. (JLT)
  //float mu=-0.2+0.07161*lm*lm/2.7+0.0015*lm*lm*lm;

  // try first modified, but with variabl sigma
  //sig = 0.1;
  //if(logMhalo<29.9) sig = -(lm-1)*0.4 + 0.86;

  float xx=(logrho-mu)/sig;
  float pp=exp(-0.5*xx*xx);
  return(pp);
}
