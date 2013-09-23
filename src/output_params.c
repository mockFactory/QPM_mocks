#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "header.h"

/* This function is taken from Volker Springel's GADGET code and
 * modified for the parameters used for the HOD.x code.
 */

/* Ths is the same as intput_params.c, but it is called after 
 * model fitting to data, and uses the new values whenever
 * applicable.
 */

/*
 *  This function parses the parameterfile in a simple way.
 *  Each paramater is defined by a keyword (`tag'), and can be
 *  either of type douple, int, or character string.
 *  The routine makes sure that each parameter appears 
 *  exactly once in the parameterfile.
 */
void output_parameter_file(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define CHAR 4
#define LONG 4
#define MAXTAGS 300

  FILE *fd,*fdout;

  char buf[200],buf1[200],buf2[200],buf3[200];
  int  i,j,nt;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][200];
  int  errorFlag=0;
  int IDUM_MCMC_TEMP=-555;

  nt=0;

  strcpy(tag[nt],"GAMMA");
  addr[nt]=&GAMMA;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"OMEGA_M");
  addr[nt]=&OMEGA_M;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"SIGMA_8");
  addr[nt]=&SIGMA_8;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"RHO_CRIT");
  addr[nt]=&RHO_CRIT;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"ITRANS");
  addr[nt]=&ITRANS;
  id[nt++]=INT;
  
  strcpy(tag[nt],"LINEAR_PSP");
  addr[nt]=&LINEAR_PSP;
  id[nt++]=INT;
  
  strcpy(tag[nt],"KAISER");
  addr[nt]=&KAISER;
  id[nt++]=INT;
  
  strcpy(tag[nt],"BETA");
  addr[nt]=&BETA;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"SIGV");
  addr[nt]=&SIGV;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"IDUM_MCMC");
  addr[nt]=&IDUM_MCMC_TEMP;
  id[nt++]=INT;
  
  strcpy(tag[nt],"SPECTRAL_INDX");
  addr[nt]=&SPECTRAL_INDX;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"DELTA_CRIT");
  addr[nt]=&DELTA_CRIT;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"DELTA_HALO");
  addr[nt]=&DELTA_HALO;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"BOX_SIZE");
  addr[nt]=&BOX_SIZE;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"RESOLUTION");
  addr[nt]=&RESOLUTION;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"VBIAS");
  addr[nt]=&VBIAS;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"VBIAS_C");
  addr[nt]=&VBIAS_C;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"TF_file"); 
  addr[nt]=Files.TF_file;
  id[nt++]=STRING;

  strcpy(tag[nt],"M1");
  addr[nt]=&HOD.M1;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"M_min");
  addr[nt]=&HOD.M_min;
  id[nt++]=DOUBLE;
        
  strcpy(tag[nt],"M_cen_max");
  addr[nt]=&HOD.M_cen_max;
  id[nt++]=DOUBLE;
        
  strcpy(tag[nt],"M_cut");
  addr[nt]=&HOD.M_cut;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"M_max");
  addr[nt]=&HOD.M_max;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"sigma_logM");
  addr[nt]=&HOD.sigma_logM;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"MaxCen");
  addr[nt]=&HOD.MaxCen;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"alpha");
  addr[nt]=&HOD.alpha;
  id[nt++]=DOUBLE;
    
  strcpy(tag[nt],"pdfc");
  addr[nt]=&HOD.pdfc;
  id[nt++]=INT;
    
  strcpy(tag[nt],"pdfs");
  addr[nt]=&HOD.pdfs;
  id[nt++]=INT;
    
  strcpy(tag[nt],"GALAXY_DENSITY");
  addr[nt]=&GALAXY_DENSITY;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"EXCLUSION");
  addr[nt]=&EXCLUSION;
  id[nt++]=INT;

  strcpy(tag[nt],"FIX_PARAM");
  addr[nt]=&FIX_PARAM;
  id[nt++]=INT;

  strcpy(tag[nt],"POWELL");
  addr[nt]=&POWELL;
  id[nt++]=INT;

  strcpy(tag[nt],"OUTPUT");
  addr[nt]=&OUTPUT;
  id[nt++]=INT;

  for(i=1;i<=11;++i)
    {
      sprintf(tag[nt],"free[%d]",i);
      addr[nt]=&HOD.free[i];
      id[nt++]=INT;
    }

  strcpy(tag[nt],"All");
  addr[nt]=&Task.All;
  id[nt++]=INT;
    
  strcpy(tag[nt],"real_space_xi");
  addr[nt]=&Task.real_space_xi;
  id[nt++]=INT;
    
  strcpy(tag[nt],"z_space_xi");
  addr[nt]=&Task.z_space_xi;
  id[nt++]=INT;
    
  strcpy(tag[nt],"kaiser_xi");
  addr[nt]=&Task.kaiser_xi;
  id[nt++]=INT;
    
  strcpy(tag[nt],"multipoles");
  addr[nt]=&Task.multipoles;
  id[nt++]=INT;
    
  strcpy(tag[nt],"r_half");
  addr[nt]=&Task.r_half;
  id[nt++]=INT;
    
  strcpy(tag[nt],"wp_minimize");
  addr[nt]=&Task.wp_minimize;
  id[nt++]=INT;
  Task.wp_minimize=0;

  strcpy(tag[nt],"zspace_minimize");
  addr[nt]=&Task.zspace_minimize;
  id[nt++]=INT;
  Task.zspace_minimize=0;

  strcpy(tag[nt],"MCMC");
  addr[nt]=&Task.MCMC;
  id[nt++]=INT;
  Task.MCMC=0;

  strcpy(tag[nt],"COVAR");
  addr[nt]=&COVAR;
  id[nt++]=INT;
  COVAR=1;
    
  strcpy(tag[nt],"DEPROJECTED");
  addr[nt]=&DEPROJECTED;
  id[nt++]=INT;
    
  strcpy(tag[nt],"fname_covar");
  addr[nt]=&wp.fname_covar;
  id[nt++]=STRING;

  strcpy(tag[nt],"pi_max");
  addr[nt]=&wp.pi_max;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"esys");
  addr[nt]=&wp.esys;
  id[nt++]=DOUBLE;
 
  strcpy(tag[nt],"fname_wp");
  addr[nt]=&wp.fname_wp;
  id[nt++]=STRING;

  strcpy(tag[nt],"wp_format");
  addr[nt]=&wp.format;
  id[nt++]=INT;

 strcpy(tag[nt],"root_filename");
  addr[nt]=&Task.root_filename;
  id[nt++]=STRING;
   
  strcpy(tag[nt],"CVIR_FAC");
  addr[nt]=&CVIR_FAC;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"JENKINS_A");
  addr[nt]=&JENKINS_A;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"JENKINS_C");
  addr[nt]=&JENKINS_C;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"JENKINS_B");
  addr[nt]=&JENKINS_B;
  id[nt++]=DOUBLE;

  strcpy(tag[nt],"BEST_FIT");
  addr[nt]=&BEST_FIT;
  id[nt++]=INT;


  if((fd=fopen(fname,"r")))
    {
      sprintf(buf,"%s.new",fname);
      if(!(fdout=fopen(buf,"w")))
	{
	  fprintf(stdout,"error opening file '%s' \n",buf);
	  errorFlag=1; 
	}
      else
	{
	  while(!feof(fd))
	    {
	      fgets(buf,200,fd);
	      if(sscanf(buf,"%s%s%s",buf1,buf2,buf3)<2)
		continue;
	      
	      if(buf1[0]=='%')
		continue;
	      
	      for(i=0,j=-1;i<nt;i++)
		if(strcmp(buf1,tag[i])==0)
		  {
		    j=i;
		    tag[i][0]=0;
		    break;
		  }
	      
	      if(j>=0)
		{
		  switch(id[j])
		    {
		    case DOUBLE:
		      if(!((double)atof(buf2)!=*((double*)addr[j])))
			*((double*)addr[j])=atof(buf2); 
		      fprintf(fdout,"%-35s%g\n",buf1,*((double*)addr[j]));
		      break;
		    case STRING:
		      strcpy(addr[j],buf2);
		      fprintf(fdout,"%-35s%s\n",buf1,buf2);
		      break;
		    case INT:
		      if(!(atoi(buf2)!=*((int*)addr[j])))
			*((int*)addr[j])=atoi(buf2);
		      fprintf(fdout,"%-35s%d\n",buf1,*((int*)addr[j]));
		      break;
		    case CHAR:
		      if(!(buf2[0]!=*((char*)addr[j])))
			*((char*)addr[j])=buf2[0];
		      fprintf(fdout,"%-35s%c\n",buf1,*((int*)addr[j]));
		      break;
		    }
		}
	    }
	}
    }
  else
    {
      fprintf(stderr,"Parameter file %s not found.\n", fname);
      exit(1);
    }

  fprintf(fdout,"\n\n");
  fclose(fd);
  fclose(fdout);


#undef DOUBLE 
#undef STRING 
#undef INT 
#undef MAXTAGS
}
 

