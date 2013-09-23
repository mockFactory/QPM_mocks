/* PROGRAM NBRSFIND2

   --- nbrsfind2(smin,smax,rmax,nmesh,xpos,ypos,zpos,nbrmax,indx,rsqr,x,y,z,meshparts,meshstart)
   --- nbrsfind2(x,y,z,smin,smax,rmax,meshparts,meshstart,nmesh,xpos,ypos,zpos,indx,rsqr,nbrmax)
   --- find all neighbours of a particle within a radius rmax.
   --- adapted from DHW's rfind.c
   --- version2 has all arrays pushed to the end

   * x,y,z = arrays of particle coordinates.
   * smin,smax = particles are located in a box running from 
                (smin,smin,smin) to (smax,smax,smax).
   * rmax = max. radius within which to search for neighbors.
   * meshparts = particle link list produced by `meshlink'
   * meshstart = lattice pointer array produced by `meshlink'
   * nmesh = dimension of mesh array in each dimension.
   * xpos,ypos,zpos = coordinates defining center of search.
   * indx = on return, contains list of particles within rmax of
       center. allocate as (int *)calloc(nbrmax,sizeof(int))
   * rmax = on return, rsqr[i] is the squared distance of particle
       indx[i]. allocate as (float *)calloc(nbrmax,sizeof(float)).
   * nbrmax = on input: dimension of indx[] and rsqr[], 
       for error checking on return: number of particles within rmax
       
   Notes: 09/10/96 - assumes periodic boundary condition.
   
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898

/*
void nbrsfind2(smin,smax,rmax,nmesh,xpos,ypos,zpos,nbrmax,indx,rsqr,x,y,z,meshparts,meshstart)

float smin,smax,rmax;
float *x,*y,*z;
int *meshparts,***meshstart;
int nmesh;
float xpos,ypos,zpos;
int *indx;
float *rsqr;
int *nbrmax;
*/

void nbrsfind2(float smin,float smax,float rmax,int nmesh,float xpos,float ypos,float zpos,
               int *nbrmax,int *indx,float *rsqr,float *x,float *y,float *z,
               int *meshparts,int ***meshstart,int ip)                         
{
   int i,j,k,ir;
   int ix,iy,iz,
       iix,iiy,iiz,
       iiix,iiiy,iiiz,
       nbr,p;
   float rmax2,r2,side,side2,sinv,dx,dy,dz;

   side=(smax-smin);
   side2=side/2.;
   sinv=1./side;
   
   /*   printf("nbrsfind2> %f %f %f %f %f %f %d %d\n",smin,smax,rmax,xpos,ypos,zpos,nmesh,*nbrmax);
    */

   /*   ix=(int)(nmesh*(xpos-smin)*sinv);
   iy=(int)(nmesh*(ypos-smin)*sinv);
   iz=(int)(nmesh*(zpos-smin)*sinv);

   if (ix>=nmesh||iy>=nmesh||iz>=nmesh ||ix<0||iy<0||iz<0){
     fprintf(stderr,"nbrsfind2> error in position or lattice parameters\n");
     fprintf(stderr,"nbrsfind2> nmesh, ix, iy, iz = %d %d %d %d\n",nmesh,ix,iy,iz) ;
     exit(-1) ;
   }
*/

   ix=(int)(nmesh*(xpos-smin)*sinv);
   if(ix>=nmesh)
     {
       fprintf(stderr,"meshlink> Warning: Particle at x = %f ix = %d\n",xpos,ix);
       ix=ix-1;
     }


   iy=(int)(nmesh*(ypos-smin)*sinv);
   if(iy>=nmesh)
     {
       fprintf(stderr,"meshlink> Warning: Particle at y = %f iy = %d\n",ypos,iy);
       iy=iy-1;
     }


   iz=(int)(nmesh*(zpos-smin)*sinv);
   if(iz>=nmesh)
     {
       fprintf(stderr,"meshlink> Warning: Particle at z = %f iz = %d\n",zpos,iz);
       iz=iz-1;
     }

   
   rmax2=rmax*rmax;
   nbr=0;
   ir=(int)(nmesh*rmax*sinv)+1;

   for(iix=-ir;iix<=ir;iix++)
   for(iiy=-ir;iiy<=ir;iiy++)
   for(iiz=-ir;iiz<=ir;iiz++){
     iiix=(ix+iix+nmesh)%nmesh ;
     iiiy=(iy+iiy+nmesh)%nmesh ;
     iiiz=(iz+iiz+nmesh)%nmesh ;
     p=meshstart[iiix][iiiy][iiiz] ;
     
     while(p>=0){
       dx=mabs(xpos-x[p]) ;
       dy=mabs(ypos-y[p]) ;
       dz=mabs(zpos-z[p]) ;
       if (dx>side2)  dx=side-dx ;
       if (dy>side2)  dy=side-dy ;
       if (dz>side2)  dz=side-dz ;
       r2=dx*dx+dy*dy+dz*dz ;
       
       if(r2<=rmax2) {
	 indx[nbr]=p;
	 rsqr[nbr]=r2;
	 nbr++;
       }
       if(nbr>*nbrmax){
	 fprintf(stderr,"nbrsfind2>too many particles in indx list\n");
         fprintf(stderr,"nbrsfind2>reset nbrmax and try again\n");
	 fprintf(stderr,"nbr = %d\n",nbr);
	 exit(-1) ;
       }
       p=meshparts[p];
     }
   }
   *nbrmax=nbr;
}
