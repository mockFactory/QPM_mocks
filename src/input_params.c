#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "header.h"

/* This function is taken from Volker Springel's GADGET code and
 * modified for the parameters used for the HOD.x code.
 */

/*
 *  This function parses the parameterfile in a simple way.
 *  Each paramater is defined by a keyword (`tag'), and can be
 *  either of type douple, int, or character string.
 *  The routine makes sure that each parameter appears 
 *  exactly once in the parameterfile.
 */
void read_parameter_file(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define CHAR 4
#define LONG 4
#define MAXTAGS 300

  FILE *fd,*fdout;

  char buf[200],buf1[200],buf2[200],buf3[200],tempchar;
  int  i,j,nt,ii,nn,ctemp;
  int  id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][200];
  int  errorFlag=0;
  int IDUM_MCMC_TEMP=-555;

  nt=0;
  
  strcpy(tag[nt],"NO_FOF_HALOS");
  addr[nt]=&NO_FOF_HALOS;
  id[nt++]=INT;
  
  strcpy(tag[nt],"populate_simulation");
  addr[nt]=&Task.populate_simulation;
  id[nt++]=INT;
  
  strcpy(tag[nt],"create_halos");
  addr[nt]=&Task.create_halos;
  id[nt++]=INT;
  
  SIGV = 0;
  strcpy(tag[nt],"SIGV");
  addr[nt]=&SIGV;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"GAMMA");
  addr[nt]=&GAMMA;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"OMEGA_M");
  addr[nt]=&OMEGA_M;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"OMEGA_B");
  addr[nt]=&OMEGA_B;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"SIGMA_8");
  addr[nt]=&SIGMA_8;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"RHO_CRIT");
  addr[nt]=&RHO_CRIT;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"HUBBLE");
  addr[nt]=&HUBBLE;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"REDSHIFT");
  addr[nt]=&REDSHIFT;
  id[nt++]=DOUBLE;
  
  strcpy(tag[nt],"ITRANS");
  addr[nt]=&ITRANS;
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

  HOD.MaxCen=1;
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

  strcpy(tag[nt],"OUTPUT");
  addr[nt]=&OUTPUT;
  id[nt++]=INT;

  strcpy(tag[nt],"HaloFile");
  addr[nt]=&Files.HaloFile;
  id[nt++]=STRING;

  strcpy(tag[nt],"FOFHaloFile");
  addr[nt]=&Files.FOFHaloFile;
  id[nt++]=STRING;

  strcpy(tag[nt],"PBHaloFile");
  addr[nt]=&Files.PBHaloFile;
  id[nt++]=STRING;

  strcpy(tag[nt],"pmfile");
  addr[nt]=&Files.pmfile;
  id[nt++]=STRING;

  strcpy(tag[nt],"root_filename");
  addr[nt]=&Task.root_filename;
  id[nt++]=STRING;

  strcpy(tag[nt],"CVIR_FAC");
  addr[nt]=&CVIR_FAC;
  id[nt++]=DOUBLE;

  if((fd=fopen(fname,"r")))
    {
      nn=filesize(fd);
      sprintf(buf,"%s","hod-usedvalues");
      if(!(fdout=fopen(buf,"w")))
	{
	  fprintf(stdout,"error opening file '%s' \n",buf);
	  errorFlag=1; 
	}
      else
	{
	  /*while(!feof(fd))*/
	  for(ii=1;ii<=nn;++ii)
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
		      *((double*)addr[j])=atof(buf2); 
		      fprintf(fdout,"%-35s%g\n",buf1,*((double*)addr[j]));
		      break;
		    case STRING:
		      strcpy(addr[j],buf2);
		      fprintf(fdout,"%-35s%s\n",buf1,buf2);
		      break;
		    case INT:
		      *((int*)addr[j])=atoi(buf2);
		      fprintf(fdout,"%-35s%d\n",buf1,*((int*)addr[j]));
		      break;
		    case CHAR:
		      *((char*)addr[j])=buf2[0];
		      fprintf(fdout,"%-35s%c\n",buf1,*((int*)addr[j]));
		      break;
		    }
		}
	      else
		{
		  fprintf(stderr,"Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			  fname,buf1);
		  errorFlag=1;
		}
	    }
	}
      fclose(fd);
      fclose(fdout);

    }
  else
    {
      fprintf(stderr,"Parameter file %s not found.\n", fname);
      exit(1);
    }


  /* Check the params for some basic interpretation
   */
  MASS_PER_PARTICLE=pow(RESOLUTION,3.0)*RHO_CRIT*OMEGA_M;
  if(HOD.M_min<1.0E4 && HOD.M_min>0)
    {
      HOD.M_min=pow(RESOLUTION,3.0)*RHO_CRIT*OMEGA_M*HOD.M_min;
      fprintf(stderr,"HOD.M_min= %e\n",HOD.M_min);
    }
  if(HOD.M1<1.0E5)
    {
      HOD.M1=HOD.M1*MASS_PER_PARTICLE;
      fprintf(stderr,"HOD.M1= %e\n",HOD.M1);
    }

  /* SCale the M_max value by OMEGA_M=0.3
   */
  /*
  HOD.M_max*=(OMEGA_M/0.3);
  */

  for(i=0;i<nt;i++)
    {      
      if(!strcmp(tag[i],"FOFHaloFile"))continue;
      if(!strcmp(tag[i],"PBHaloFile"))continue;
      if(!strcmp(tag[i],"NO_FOF_HALOS"))continue;
      if(!strcmp(tag[i],"SMF_CORRECTION"))continue;
      if(!strcmp(tag[i],"FISHER"))continue;
      if(!strcmp(tag[i],"JPL_FLAG"))continue;
      if(!strcmp(tag[i],"DONT_FIT_CLUSTERING"))continue;
      if(!strcmp(tag[i],"DONT_FIT_SMF"))continue;
      if(!strcmp(tag[i],"DONT_FIT_LENSING"))continue;
      if(!strcmp(tag[i],"STEPFAC"))continue;
      if(!strcmp(tag[i],"LENSING_OUTPUT_FLAG"))continue;
      if(!strcmp(tag[i],"DENSITY_DEPENDENCE"))continue;
      if(!strcmp(tag[i],"DENSITY_THRESHOLD") && !DENSITY_DEPENDENCE)continue;
      if(!strcmp(tag[i],"M_min_fac") && !DENSITY_DEPENDENCE)continue;
      if(!strcmp(tag[i],"HaloDensityFile") && !DENSITY_DEPENDENCE)continue;
      
      if(!strcmp(tag[i],"HOD2.M_min") && !XCORR)continue;
      if(!strcmp(tag[i],"HOD2.M_max") && !XCORR)continue;
      if(!strcmp(tag[i],"HOD2.MaxCen") && !(HOD2.pdfc>=4 && HOD.pdfc<=6))continue;
      if(!strcmp(tag[i],"HOD2.M_cut") && HOD2.pdfs<2)continue;
      if(!strcmp(tag[i],"HOD2.M1") && !XCORR)continue;
      if(!strcmp(tag[i],"HOD2.M_cen_max") && HOD2.pdfc!=7)continue;
      if(!strcmp(tag[i],"HOD2.sigma_logM") && (!XCORR || (HOD2.pdfc==1 || HOD2.pdfc==7)))continue;
      if(!strcmp(tag[i],"HOD2.pdfs") && !XCORR)continue;
      if(!strcmp(tag[i],"HOD2.pdfc") && !XCORR)continue;
      if(!strcmp(tag[i],"HOD2.alpha") && !XCORR)continue;
      if(!strcmp(tag[i],"HOD2.alpha1") && HOD2.pdfs!=4 && HOD2.pdfs!=5)continue;
      if(!strcmp(tag[i],"HOD2.M_sat_break") && HOD2.pdfs!=4 && HOD2.pdfs!=5)continue;
      if(!strcmp(tag[i],"GALAXY_DENSITY2") && !XCORR)continue;
      if(!strcmp(tag[i],"XCORR"))continue;

      if(!strcmp(tag[i],"alpha1"))continue;
      if(!strcmp(tag[i],"M_sat_break"))continue;
      
      if(!strcmp(tag[i],"OMEGA_TEMP"))
	{
	  OMEGA_TEMP = OMEGA_M;
	  continue;
	}
      if(!strcmp(tag[i],"pdf"))continue;
      if(!strcmp(tag[i],"REDSHIFT"))continue;
      if(!strcmp(tag[i],"WP_ONLY"))continue;
      if(!strcmp(tag[i],"RESTART"))continue;
      if(!strcmp(tag[i],"RESTART_FILE"))continue;
      if(!strcmp(tag[i],"DEPROJECTED"))continue;
      if(!strcmp(tag[i],"POWELL"))continue;
      if(!strcmp(tag[i],"LINEAR_PSP"))continue;
      if(!strcmp(tag[i],"BOX_SIZE"))continue;
      if(!strcmp(tag[i],"RESOLUTION"))continue;
      if(!strcmp(tag[i],"DELTA_HALO"))continue;
      if(!strcmp(tag[i],"VBIAS"))continue;
      if(!strcmp(tag[i],"COVAR"))continue;
      if(!strcmp(tag[i],"VBIAS_C"))continue;
      if(!strcmp(tag[i],"CVIR_FAC"))continue;
      if(!strcmp(tag[i],"ITRANS"))continue;
      if(!strcmp(tag[i],"IDUM_MCMC"))continue;
      if(!strcmp(tag[i],"FIX_PARAM"))continue;
      if(!strcmp(tag[i],"DEPROJECTED"))continue;
      if(!strcmp(tag[i],"OUTPUT"))continue;
      if(!strcmp(tag[i],"JENKINS_A"))continue;
      if(!strcmp(tag[i],"JENKINS_B"))continue;
      if(!strcmp(tag[i],"JENKINS_C"))continue;
      if(!strcmp(tag[i],"BEST_FIT"))continue;
      if(!strcmp(tag[i],"KAISER"))continue;
      if(!strcmp(tag[i],"SIGV"))continue;
      if(!strcmp(tag[i],"BETA"))continue;
      if(!strcmp(tag[i],"HUBBLE"))continue;
      if(!strcmp(tag[i],"OMEGA_B"))continue;

      if(!strcmp(tag[i],"MCMC"))continue;
      if(!strcmp(tag[i],"wp_minimize"))continue;
      if(!strcmp(tag[i],"wp_format"))continue;
      if(!strcmp(tag[i],"n_wp"))continue;
      if(!strcmp(tag[i],"wp_npca"))continue;
      if(!strcmp(tag[i],"zspace_minimize"))continue;
      if(!strcmp(tag[i],"pi_max"))continue;
      if(!strcmp(tag[i],"esys"))continue;

      if(HOD.color)
	{
	  if(!strcmp(tag[i],"fblue0_cen") ||
	     !strcmp(tag[i],"sigma_fblue_cen") ||
	     !strcmp(tag[i],"fblue0_sat") ||
	     !strcmp(tag[i],"sigma_fblue_sat"))
	    {
	      fprintf(stderr,"Parameters for color HOD not specified.\n");
	      exit(0);
	    }
	  continue;
	}
	     
      if(!strcmp(tag[i],"color"))continue;
      if(!strcmp(tag[i],"fblue0_cen"))continue;
      if(!strcmp(tag[i],"sigma_fblue_cen"))continue;
      if(!strcmp(tag[i],"fblue0_sat"))continue;
      if(!strcmp(tag[i],"sigma_fblue_sat"))continue;


      if(!strcmp(tag[i],"free[0]"))continue;
      if(!strcmp(tag[i],"free[1]"))continue;
      if(!strcmp(tag[i],"free[2]"))continue;
      if(!strcmp(tag[i],"free[3]"))continue; 
      if(!strcmp(tag[i],"free[4]"))continue;
      if(!strcmp(tag[i],"free[5]"))continue;
      if(!strcmp(tag[i],"free[6]"))continue;
      if(!strcmp(tag[i],"free[7]"))continue;
      if(!strcmp(tag[i],"free[8]"))continue;
      if(!strcmp(tag[i],"free[9]"))continue;
      if(!strcmp(tag[i],"free[10]"))continue;
      if(!strcmp(tag[i],"free[11]"))continue;
      if(!strcmp(tag[i],"free[12]"))continue;
      if(!strcmp(tag[i],"free[13]"))continue;
      if(!strcmp(tag[i],"free[14]"))continue;
      if(!strcmp(tag[i],"free[15]"))continue;
      if(!strcmp(tag[i],"free[16]"))continue;

      if(!strcmp(tag[i],"All"))continue;
      if(!strcmp(tag[i],"populate_sim"))continue;
      //      if(!strcmp(tag[i],"HaloFile"))continue;
      if(!strcmp(tag[i],"HOD"))continue;
      if(!strcmp(tag[i],"PCA"))continue;
      if(!strcmp(tag[i],"PVD"))continue;
      if(!strcmp(tag[i],"matter_xi"))continue;
      if(!strcmp(tag[i],"matter_pk"))continue;
      if(!strcmp(tag[i],"sigma_r"))continue;
      if(!strcmp(tag[i],"massfunc"))continue;
      if(!strcmp(tag[i],"kaiser_xi"))continue;
      if(!strcmp(tag[i],"cvir"))continue;

      if(!strcmp(tag[i],"TF_file"))
	{
	  if(ITRANS==11) {
	    sprintf(Files.TF_file,"CMBFAST_trans.dat");
	    fprintf(stderr,"No transfer function file, using [%s]\n",Files.TF_file);
	  }
	  continue;
	}

      if(!strcmp(tag[i],"fname_covar"))
	{
	  if(Task.wp_minimize) {
	    fprintf(stderr,"No filename specified for covariance matrix.\n");
	    errorFlag=1;
	  }
	  continue;
	}
      if(!strcmp(tag[i],"fname_wp"))
	{
	  if(Task.wp_minimize) {
	    fprintf(stderr,"No filename specified for wp data.\n");
	    errorFlag=1;
	  }
	  continue;
	}
      if(!strcmp(tag[i],"M_cut"))
	{
	  if(HOD.pdfs==2 || HOD.pdfs==3){
	    fprintf(stderr,"No value for M_cut given for pdfs= 2/3\n");
	    errorFlag=1;
	  }
	  continue;
	}
      if(!strcmp(tag[i],"M_cen_max"))
	{
	  if(HOD.pdfc==7) {
	    fprintf(stderr,"No value for M_cen_max given for pdfc= 7\n");
	    errorFlag=1;
	  }
	  continue;
	}
      if(!strcmp(tag[i],"sigma_logM"))
	{
	  if(HOD.pdfc==2){
	    fprintf(stderr,"No value for sigma_logM given for pdfc=2\n");
	    errorFlag=1;
	  }
	  continue;
	}
      if(!strcmp(tag[i],"MaxCen"))
	{
	  if(HOD.pdfc==5){
	    fprintf(stderr,"No value for MaxCen given for pdfc=5\n");
	    errorFlag=1;
	  }
	  continue;
	}
      if(*tag[i])
	{
	  fprintf(stderr,"Error. I miss a value for tag '%s' in parameter file '%s'.\n",
		  tag[i],fname);
	  errorFlag=1;
	}
    }


  if(errorFlag)
    endrun("error in input_params ");


#undef DOUBLE 
#undef STRING 
#undef INT 
#undef MAXTAGS
}
 

