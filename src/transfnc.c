/* PROGRAM TRANSFERFUNCTION

   --- transfnc(xk)
   --- compute the transfer function T(k) of the power spectrum
   --- T(k) defined as P(k) = k^{xindx}*T^{2}(k)

       * itrans=type of transfer function
         0  -> no change (returns 1)
	 3  -> BBKS 
         4  -> Efstathiou, Bond & White transfer function with Gamma as
              specified (eqn. 7)
	 5  -> Eisnstein & Hu     
	 11 -> read in TF from file (usually CMBFAST)

   NOTE: xk is in h/Mpc and is defined as k=2pi/lambda (not k=1/lambda)
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "header.h"

/* EBW CDM parameters */

#define aa  6.4
#define bb  3.0
#define cc  1.7
#define xnu  1.13
#define twopi  6.283185
#define xnuinv  -0.884956

double transfnc(double xk) 
{
  double transf;
  double q,t1,t2;
 
  if(xk==0.)
    {
      transf=1.;
      return (double)transf;
    }

  switch(ITRANS)
    {
    case 0:
      transf=1.;
      break;
      
    case 3:
      q = xk/((OMEGA_M)*HUBBLE);
      transf = log(1+2.34*q)/2.34/q*pow(1 + 3.89*q + pow(16.1*q,2.0) + 
				    pow(5.46*q,3.0) + pow(6.71*q,4.0),-0.25);
      break;

    case 4:
      q = xk/GAMMA;
      t1=aa*q+pow((bb*q),1.5)+(cc*q)*(cc*q);
      t2=pow(t1,xnu);
      transf=pow((1.+t2),xnuinv);
      break;

    case 5:
      transf = tf_eisenstein_hu(xk);
      break;

    case 11:
      fprintf(stderr,"transfnc>Sorry, code cannot accept tabulated input %d \n",ITRANS);
      exit(-1);
      break;

    default:
      fprintf(stderr,"transfnc> Unrecognized transfer function %d \n",ITRANS);
      exit(-1);
      break;

    }
  return (double)transf;
}
