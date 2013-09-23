/* PROGRAM MESHLINK2

   --- meshlink2(np,nmesh,smin,smax,rmax,x,y,z,meshparts,meshstart)
   --- creates a linked list for efficient neighbor searching.
   --- tag particles according to their mesh
   --- adapted from DHW's linklist.c
   --- version2 has all arrays pushed to the end
   
      * np = number of particles
      * x,y,z = arrays of particle coordinates.
      * smin,smax = particles are located in a box running from 
                    (smin,smin,smin) to (smax,smax,smax).
      * rmax = max. radius of a mesh. All neighbors closer than
               rmax. are included in a mesh.- to determine size of mesh
      * meshparts = on return, pointer to particle link list.
      * meshstart = on return, pointer to array of pointers to the
                    first particle in each mesh.
      * nmesh = on return,dimension of mesh array in each dimension.

   Notes: 09/09/96

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define sqr(x) ((x)*(x))
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define mabs(A) ((A) < 0.0 ? -(A) : (A))
#define pi 3.1415926535898
#define ind(a,b,c) (a)*n*n+(b)*n+(c)
#define ALLOC3D(type,dim0,dim1,dim2) \
        (type***)a3alloc((unsigned)dim0,(unsigned)dim1,(unsigned)dim2,(unsigned)sizeof(type))
#define NLATMAX 600      /* maximum lattice dimension */

int ***i3tensor_2(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

/*
void meshlink2(np1,nmesh,smin,smax,rmax,x1,y1,z1,meshparts,meshstart,meshfac)
int np1;
float smin,smax,rmax;
float *x1,*y1,*z1;
int **meshparts,****meshstart;
int *nmesh;
int meshfac;
*/
void meshlink2(int np1,int *nmesh,float smin,float smax,float rmax,float *x1,float *y1,float *z1,
	       int **meshparts,int ****meshstart,int meshfac)
{
   int nlatt; 
   int i,j,k;
   int ix,iy,iz;
   float sinv;

   nlatt=(int)((smax-smin)/rmax);
   if (nlatt>NLATMAX) nlatt=NLATMAX ;
   
   /*
   *meshstart=ALLOC3D(int,nlatt,nlatt,nlatt);
   */
   fprintf(stderr,"nlatt= %d %f %f %f\n",nlatt,smin,smax,rmax);
   *meshstart=(int ***)i3tensor_2(0,nlatt-1,0,nlatt-1,0,nlatt-1);
   fprintf(stderr,"done\n");
   *meshparts=(int *)calloc(np1,sizeof(int));

   for(i=0;i<np1;i++)(*meshparts)[i]=-1;
   for(i=0;i<nlatt;i++)
     for(j=0;j<nlatt;j++)
       for(k=0;k<nlatt;k++)
	 (*meshstart)[i][j][k]=-1;

   sinv=1./(smax-smin);
   for(i=0;i<np1;i++){
     ix=(int)(nlatt*(x1[i]-smin)*sinv);
     iy=(int)(nlatt*(y1[i]-smin)*sinv);
     iz=(int)(nlatt*(z1[i]-smin)*sinv);
     if (ix==nlatt) ix=nlatt-1;
     if (iy==nlatt) iy=nlatt-1;
     if (iz==nlatt) iz=nlatt-1;
     (*meshparts)[i]=(*meshstart)[ix][iy][iz];
     (*meshstart)[ix][iy][iz]=i;
   }

   *nmesh=nlatt;
   fprintf(stderr,"Done with meshlink.  nlatt=%d\n",nlatt);
}


