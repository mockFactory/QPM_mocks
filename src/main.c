#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

/* test file routines.
 */
void test(void);
void chi2_grid(int argc, char **argv);
void fit_scale_bias(int argc, char **argv);
void aspen_breakout(void);
void populate_simulation_clf(void);

int main(int argc, char **argv)
{
  double s1;
  int i, irank, nrank;

  ARGC = argc;
  ARGV = argv;

  OUTPUT=0;
  if(argc==1)
    endrun("./QPM.mock qpm.bat_file > output");

  read_parameter_file(argv[1]);

  SIGMA_8Z0 = SIGMA_8;
  SIGMA_8 = SIGMA_8*growthfactor(REDSHIFT);
  fprintf(stdout,"SIGMA_8(Z=%.2f)= %.3f\n",REDSHIFT,SIGMA_8);
  RESET_COSMOLOGY++;

  if(argc>2)
    {
      if(atoi(argv[2])==999)
	test();
      else
	SUBFRAC = atof(argv[2]);
    }
  
  if(Task.create_halos)
    create_lognormal_halos();
  if(Task.populate_simulation)
    populate_simulation_hod();
  exit(0);

}

