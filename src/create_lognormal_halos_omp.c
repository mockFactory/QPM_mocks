#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "header.h"

#define NBINLOOKUP 10000
#define cnint(x) ((x-floor(x)) < 0.5 ? floor(x) : ceil(x))


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

void create_lognormal_halos()
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
  float **haloarr;
  int ngal,nsati[9],ALL_FILES=0,TRACK_GALAXIES=0,WARREN_MASS_CORRECTION=0,haloid;

  float *xt,*yt,*zt,*vxt,*vyt,*vzt;
  float lm, logMhalo, mu, sig;

  //static int nhalo, *used_halo, iuse_flag = 1;
  int ihalo=0, nhalo, npart, idum, nmax, nadd;
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
  double xran;
  
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

  float logprho, dlnMbin, prob, *pnorm, logM, lnRhoMin, lnRhoMax, *rhohist, logrho, *mux, m0;
  int NmassBin=256, jbin, NrhoBin=2048, ibin;
  double avg, ntar, ntot;

  int ichunk, nchunk, nchunk_reads, nread, ntotread=0;

  // openmp stuff
  int irank=0,nrank=1;
  float *logMx, **pnormx;
  int *iix;
  unsigned int iseed = -555;
  struct drand48_data drand_buf;

/*
  t0 = second();
  np = 100000000;
  xp = vector(0,np-1);
#pragma omp parallel private(nrank, irank, i, x1, IDUM, iseed)
  {
    nrank = omp_get_num_threads();
    fprintf(stdout, "rank: %d %d\n",irank=omp_get_thread_num(),nrank);

    iseed = iseed - irank;
    //#pragma omp for
    for(i=irank;i<np;i+=nrank)
      {
	x1 = rand_r(&iseed);
	xp[i] = exp(-x1*x1);
      }

    j = nrank;
  }
  t1 = second();
  printf("TIME %d %.2f %.2f\n",j,timediff(t0,t1),timediff(t0,t1)/j);
  exit(0);
*/
  

  lnRhoMin = log(1.0E-3);
  lnRhoMax = log(100);
  rhohist = vector(0,NrhoBin-1);
  for(i=0;i<NrhoBin;++i)rhohist[i] = 0;

  ASCII_OUTPUT = 0;

  
  MASS_MIN = log(1.0E12);
  MASS_MAX = log(HOD.M_max); 
  NmassBin = (int)((MASS_MAX-MASS_MIN)/0.05+1);
  dlnMbin = (MASS_MAX-MASS_MIN)/NmassBin; 
  pnorm = vector(0,NmassBin-1);


  fp = fopen(Files.pmfile,"r");
  fprintf(stderr,"%s\n",Files.pmfile);

  // get header
  fread(&idum,sizeof(int),1,fp);
  fread(&idum,sizeof(int),1,fp);
  fread(&tpmhdr,1,idum,fp);
  np = tpmhdr.npart;
  printf("particles in TPM file: %d %d\n",tpmhdr.npart,np);
  printf("aexp: %f\n",tpmhdr.aa);
  fflush(stdout);
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
  nchunk = np/10;
  positions = malloc(((long)nchunk)*3*(long)sizeof(float));
  printf("done allocation of temp array. nchunk= %d %d\n",nchunk,np);
  fflush(stdout);

  nchunk_reads = np/nchunk + 1;
  if(np%nchunk==0)nchunk_reads--;

  ii = 1;
  for(ichunk=1;ichunk<=nchunk_reads;++ichunk)
    {
      nread = nchunk;
      if(ntotread + nread>np) nread = np-ntotread;
      fprintf(stderr,"%d %d %d %d %d\n",ichunk,ii,nread,ntotread,np);
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
  fprintf(stderr,"Read %d positions\n",ntotread);

  // now do the same thing for the velocities
  ii = 1;
  ntotread = 0;
  vfac = sqrt(OMEGA_M*pow(1+REDSHIFT,3.0)+(1-OMEGA_M))*100/(1+REDSHIFT);
  if(SIGV>0)
    fprintf(stdout,"Adding SIGV=%f random motions to all particles\n",SIGV);
  for(ichunk=1;ichunk<=nchunk_reads;++ichunk)
    {
      fprintf(stderr,"velread: chunk %d/%d\n",ichunk,nchunk_reads);
      nread = nchunk;
      if(ntotread + nread>np) nread = np-ntotread;
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


  fread(&density[1],sizeof(float),np,fp);
  fclose(fp);
  free(positions);

  if(SUBFRAC>0)
    {
      if(ARGC>4)
	fp = fopen(ARGV[4],"w");
      else 
	fp = fopen("pm_ascii.dat","w");
      n = 0;
      for(i=1;i<=np;++i)
	if(drand48()<SUBFRAC) { n++;
	  fprintf(fp,"%f %f %f %f %f %f %f %d\n",xp[i],yp[i],zp[i],vxp[i],vyp[i],vzp[i],density[i],n); }
      fprintf(stdout,"selected %d particles\n",n);
      exit(0);
    }

  //REDSHIFT = 0; // do we keep this? TODO
  RESET_COSMOLOGY++;

  // how many halos do we need?
  nexp = qromo(func_halo_density,MASS_MIN,MASS_MAX,midpnt)*BOX_SIZE*BOX_SIZE*BOX_SIZE;
  fraction = nexp/np;
  fprintf(stderr,"Fraction of particles to use: %e\n",fraction);

  // how many halos are we going to output
  x1 = qromo(func_halo_density,log(1.0E12),MASS_MAX,midpnt)*BOX_SIZE*BOX_SIZE*BOX_SIZE;
  nhalo = nexp*1.01; 
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
    
  // loop for histogramming the rhos over all sampled particles
  ii = 0;
  for (j=1; j<=np; j++) {
    if (density[j]>=0) {
      logprho = log(density[j]);
      if (logprho>lnRhoMin && logprho<lnRhoMax) {
	ibin=(int)( NrhoBin*(log(density[j])-lnRhoMin)/(lnRhoMax-lnRhoMin) );
	rhohist[ibin] += 1;
	ii++;
      }
    }
  }
  fprintf(stderr,"Done getting histogram of densities. %d %d\n",np,ii);

  // get the normalization for each bin

  // First compute the "average" number of particles selected per mass bin,
  // then the normalization is the ratio to what we want.
  for (jbin=0; jbin<NmassBin; jbin++) {
    logM = MASS_MIN+dlnMbin*(jbin+0.5);
    if(logM<log(1e12))continue;
    //ntar = exp(nofm->val(logM))*(MASS_MAX-MASS_MIN)/NmassBin * cc->vol; //density?
    ntar = qromo(func_halo_density, logM-0.5*dlnMbin, logM+0.5*dlnMbin, midpnt)*
      BOX_SIZE*BOX_SIZE*BOX_SIZE;
    ntot += ntar;
    avg  = 1e-30;
    for ( ibin=0; ibin<NrhoBin; ibin++) {
      logrho= lnRhoMin+(ibin+0.5)*(lnRhoMax-lnRhoMin)/NrhoBin;
      avg += PrhoMhalo(logrho,logM)*rhohist[ibin];
    }
    pnorm[jbin]=ntar/avg;
    printf("NORM: %d %e %e %e %e %e %e\n",jbin,exp(logM),pnorm[jbin],
	   PrhoMhalo(log(0.5),logM),PrhoMhalo(log(1),logM),
	   PrhoMhalo(log(2),logM),PrhoMhalo(log(4),logM));
    if (pnorm[jbin]>1.0) {
      fprintf(stderr,"Error in normalization for bin %d %e\n",jbin,pnorm[jbin]);
      fprintf(stderr,"%e %e %e\n",ntar, avg, exp(logM));
      //exit(1);
    }
  }
 //exit(0);
  fprintf(stderr,"Done getting histogram normalization. %e %e\n",nexp,ntot);

  /*----- Martin's code for converting the sampled particles into
   * masses via their densities.
   */
  fprintf(stderr,"%d\n",nrank);
  logMx = vector(0,nrank-1);
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

  // notes: when i put this inside the parallel region, it never terminates. when i put this outside the parallel region, it doesn't work. WTF. idea: try putting it inside but only executing for itask==0?
    mux = vector(0,NmassBin-1);
    for(i=0;i<NmassBin-1;++i)
      {
	lm = MASS_MIN+dlnMbin*(i+0.5);
	logMx[i] = lm;
	mux[i] = 0.4863-0.0007793*lm+0.07161*lm*lm; // martin's function
	// my function, z=1.5
	lm = exp(logMx[i]);
	m0 = pow(10,13.25);
	mux[i] = 1 + 0.13*log10(lm/m0) + pow(lm/m0,0.35)/(1+pow(lm/m0,-0.7)) - 0.5;
	//printf("%d %e %e %e\n",i,lm,logMx[i]/log(10),mux[i]);
      }


  nmax = 0;
  sig = 1/0.1;
#pragma omp parallel private(ii, logM, logprho, jbin, prob, irank, nrank, j, lm, mu,x1, iseed, sig, xran, drand_buf) \
  shared(np, xp, yp, zp, vxp, vyp, vzp, mvect, dvect, rvir, density, xxh, yxh, zxh, vxh, vyh, vzh,nhalo,pnorm,pnormx,mux)
  {
    nrank = omp_get_num_threads();
    pnormx = matrix(0,nrank-1,0,NmassBin-1);
    for(i=0;i<nrank;++i)
      for(j=0;j<NmassBin;++j)
	pnormx[i][j] = pnorm[j];

    #pragma omp master 
    {
      mux = vector(0,NmassBin-1);
      for(i=0;i<NmassBin-1;++i)
	{
	  lm = MASS_MIN+dlnMbin*(i+0.5);
	  logMx[i] = lm;
	  mux[i] = 0.4863-0.0007793*lm+0.07161*lm*lm; // martin's function
	  // my function, z=1.5
	  lm = exp(logMx[i]);
	  m0 = pow(10,13.25);
	  mux[i] = 1 + 0.13*log10(lm/m0) + pow(lm/m0,0.35)/(1+pow(lm/m0,-0.7)) - 0.5;
	  //printf("%d %e %e %e\n",i,lm,logMx[i]/log(10),mux[i]);
	}
    }
    #pragma omp barrier


    fprintf(stderr, "rank: %d %d\n",irank=omp_get_thread_num(),nrank);
    iseed = iseed - irank;
    srand48_r (iseed, &drand_buf);

    ii = irank+1;

    for(j=ii;j<=np;j+=nrank) {
      if( (j%1000000==0)&&(OUTPUT) ){ 
	drand48_r (&drand_buf, &xran);
	fprintf(stderr,"%d %e %e %d %f %e\n",j,nexp,j*1./np,jbin,ii*1./(j),xran); }
    if(density[j]<0){ printf("ZERO %d %f\n",j,density[j]); fflush(stdout); }
    if (density[j]>0) {
      logprho = log(density[j]);
      for ( jbin=NmassBin-1; jbin>=0; jbin--) {
	//logM = MASS_MIN+dlnMbin*(jbin+0.5);
	//prob = PrhoMhalo(logprho,logM) * pnorm[jbin];

	mu=mux[jbin];// 0.4863-0.0007793*lm+0.07161*lm*lm; //martin's version
	sig = 0.1;
	x1=(logprho-mu)/sig;
	prob=exp(-0.5*x1*x1)*pnorm[jbin];

	if (prob>=1.0) {
	  fprintf(stderr,"Out-of-bounds in create.\n");
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
	  //printf("BOO %d %d %d %e %e %e %d\n",irank, j, jbin, logprho, logM, prob, ii);
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
  indx=malloc(nhalo*sizeof(int));
  rsqr=malloc(nhalo*sizeof(float));
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

  for(ii=1;ii<=nhalo;++ii)
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
	  fflush(stdout);
	}
    }
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
