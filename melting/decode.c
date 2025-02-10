/******************************************************************************
 *                               MELTING v4.3                                 *
 * This program   computes for a nucleotide probe, the enthalpy, the entropy  *
 * and the melting temperature of the binding to its complementary template.  *
 * Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  *
 *          Copyright (C) Nicolas Le Novère and Marine Dumousseau  1997-2013  *
 *                                                                            *
 * File: decode.c                                                             *
 * Date: 01/APR/2009                                                          *
 * Aim : Read and treat the input parameters.                                 *
 ******************************************************************************/

/*    This program is free software; you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation; either version 2 of the License, or
      (at your option) any later version.

      This program is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
      GNU General Public License for more details.

      You should have received a copy of the GNU General Public License
      along with this program; if not, write to the Free Software
      Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

      Nicolas Le Novère
      Babraham Institute, Babraham Research Campus
      Babraham CB22 3AT Cambridge United-Kingdom.
      n.lenovere@gmail.com
       
      Marine Dumousseau
      EMBL-EBI, Wellcome-Trust Genome Campus
      Hinxton CB10 1SD Cambridge United-Kingdom. 
      marine@ebi.ac.uk 
      
*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>PREPROCESSOR INFORMATIONS<<<<<<<<<<<<<<<<<<<<<<<<*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include "common.h"
#include "decode.h"


/*******************************************************
 * Decode a string containing configuration parameters *
 *******************************************************/

struct param *decode_input(struct param *pst_in_param, char *ps_input, char *ps_path){
  
/* those four variables are used to construct the name of the outfile */
  time_t  universal_time;        
  struct tm *now;                
  int i_year, i_month, i_day, i_hour, i_min;
  char s_nmonth[4]; /* three first letters of the month (+eos)*/

  char *ps_line;
  char *ps_inputline;
  FILE *pF_INFILE;
  
  switch (ps_input[1]){
  case 'A':         /* an alternative NN set is required */
      if ( strlen(&ps_input[2]) != 0 && isalnum((int)ps_input[2]) ){
	  i_alt_nn = TRUE;
	  if (pst_in_param->pst_present_nn != NULL)
	      pst_in_param->pst_present_nn = NULL; /* Reset the NN set */
	  pst_in_param->pst_present_nn = read_nn(&ps_input[2],ps_path);
	  strncpy(pst_in_param->pst_present_nn->s_nnfile,&ps_input[2],FILE_MAX);
	  pst_in_param->pst_present_nn->s_nnfile[FILE_MAX-1] = '\0'; /* security check */
	  i_hybridtype = TRUE;	/* The entry of a NN set is equivalent to define an hybrid style */
      } else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
      }
	break;
  case 'C':	    /* a complement is furnished (seems to mean mismatches or dangling ends or inosine mismatches) */
      if ( strlen(&ps_input[2]) != 0 ){
	  i_complement = TRUE;
	  i_mismatchesneed = TRUE;
	  i_inosineneed = TRUE;
	  i_dangendsneed = TRUE;
	  if ( ( pst_in_param->ps_complement = (char *)malloc(strlen(&ps_input[2])+1) ) == NULL){
	      fprintf(ERROR," Function decode_input, line __LINE__:"
		      " Unable to allocate memory to register the complement\n");
	      exit(EXIT_FAILURE);
	  }
	  strncpy(pst_in_param->ps_complement,&ps_input[2],strlen(&ps_input[2])+1);
	  pst_in_param->ps_complement[strlen(&ps_input[2])] = '\0'; /* security check */
      } else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
    }
    break;
  case 'D':         /* an alternative dangling ends set is required */
      if ( strlen(&ps_input[2]) != 0 && isalnum((int)ps_input[2]) ){
	  i_alt_de = TRUE;
	  if (pst_in_param->pst_present_de != NULL)
	      pst_in_param->pst_present_de = NULL; /* Reset the NN set */
	  pst_in_param->pst_present_de = read_dangends(&ps_input[2],ps_path);
	  strncpy(pst_in_param->pst_present_de->s_defile,&ps_input[2],FILE_MAX); 
	  pst_in_param->pst_present_de->s_defile[FILE_MAX-1] = '\0'; /* security check */
      } else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
    }
    break;
  case 'F':        /* change correction factor for nucleic acid concentration */
      if ( strlen(&ps_input[2]) != 0 && isdigit((int)ps_input[2]) )
	  pst_in_param->d_gnat = strtod(&ps_input[2],NULL);
      else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
      }
      break;
  case 'h':       /* help required */
      usage();
      exit(EXIT_SUCCESS);
  case 'H':       /* hybridisation type, max 6 characters*/
      /* CAUTION strncmp sends '0' when identical */
      /* All the single character cases are here for compatibility */
      /* with version < 4,  However they are deprecated. */
      i_hybridtype = TRUE; 
      if( strncmp(&ps_input[2],"dnadna",6) == 0 
	  || strncmp(&ps_input[2],"A",6) == 0){
	  if (pst_in_param->pst_present_nn == NULL){
	      pst_in_param->pst_present_nn = read_nn(DEFAULT_DNADNA_NN,ps_path);
	      strncpy(pst_in_param->pst_present_nn->s_nnfile,DEFAULT_DNADNA_NN,FILE_MAX);
	  } 
	  i_dnadna = TRUE;
	  i_dnarna = FALSE;
	  i_rnarna = FALSE;
      } else if( strncmp(&ps_input[2],"dnarna",6) == 0 
		 || strncmp(&ps_input[2],"rnadna",6) == 0 
		 || strncmp(&ps_input[2],"B",6) == 0){
	  if (pst_in_param->pst_present_nn == NULL){
	      pst_in_param->pst_present_nn = read_nn(DEFAULT_DNARNA_NN,ps_path);
	      strncpy(pst_in_param->pst_present_nn->s_nnfile,DEFAULT_DNARNA_NN,FILE_MAX);
	  } 
	  i_dnadna = FALSE;
	  i_dnarna = TRUE;
	  i_rnarna = FALSE;
      }
      else if( strncmp(&ps_input[2],"rnarna",6) == 0 
	       ||strncmp(&ps_input[2],"C",6) == 0){
	  if (pst_in_param->pst_present_nn == NULL){
	      pst_in_param->pst_present_nn = read_nn(DEFAULT_RNARNA_NN,ps_path);
	      strncpy(pst_in_param->pst_present_nn->s_nnfile,DEFAULT_RNARNA_NN,FILE_MAX); 
	  }
	  i_dnadna = FALSE;
	  i_dnarna = FALSE;
	  i_rnarna = TRUE;
    }
      else if( strncmp(&ps_input[2],"F",2) == 0){ 
	  /* compare 2 letters because EOS make sure it's not just the first letters of a word */
	  if (pst_in_param->pst_present_nn == NULL){
	      pst_in_param->pst_present_nn = read_nn("fre86a.nn",ps_path);
	      strncpy(pst_in_param->pst_present_nn->s_nnfile,"fre86a.nn",FILE_MAX);
	  } 
	  i_dnadna = FALSE;
	  i_dnarna = FALSE;
	  i_rnarna = TRUE;
      }
      else if( strncmp(&ps_input[2],"R",2) == 0){
	  if (pst_in_param->pst_present_nn == NULL){
	      pst_in_param->pst_present_nn = read_nn("bre86a.nn",ps_path);
	      strncpy(pst_in_param->pst_present_nn->s_nnfile,"bre86a.nn",FILE_MAX); 
	  }
	  i_dnadna = TRUE;
	  i_dnarna = FALSE;
	  i_rnarna = FALSE;
      }
      else if( strncmp(&ps_input[2],"S",2) == 0){
	  if (pst_in_param->pst_present_nn == NULL){
	      pst_in_param->pst_present_nn = read_nn("sug96a.nn",ps_path);
	      strncpy(pst_in_param->pst_present_nn->s_nnfile,"sug96a.nn",FILE_MAX);
	  } 
	  i_dnadna = TRUE;
	  i_dnarna = FALSE;
	  i_rnarna = FALSE;
      }
      else if( strncmp(&ps_input[2],"T",2) == 0){
	  if (pst_in_param->pst_present_nn == NULL){
	      pst_in_param->pst_present_nn = read_nn("san96a.nn",ps_path);
	      strncpy(pst_in_param->pst_present_nn->s_nnfile,"san96a.nn",FILE_MAX); 
	  }
	  i_dnadna = TRUE;
	  i_dnarna = FALSE;
	  i_rnarna = FALSE;
      }
      else if( strncmp(&ps_input[2],"U",2) == 0){
	  if (pst_in_param->pst_present_nn == NULL){
	      pst_in_param->pst_present_nn = read_nn("sug95a.nn",ps_path);
	      strncpy(pst_in_param->pst_present_nn->s_nnfile,"sug95a.nn",FILE_MAX); 
	  }
	  i_dnadna = FALSE;
	  i_dnarna = TRUE;
	  i_rnarna = FALSE;
    }
      else if( strncmp(&ps_input[2],"W",2) == 0){
	  if (pst_in_param->pst_present_nn == NULL){
	      pst_in_param->pst_present_nn = read_nn("all97a.nn",ps_path);
	      strncpy(pst_in_param->pst_present_nn->s_nnfile,"all97a.nn",FILE_MAX);
	  } 
	  i_dnadna = TRUE;
	  i_dnarna = FALSE;
	  i_rnarna = FALSE;
      }
      else{
	  fprintf(ERROR," I did not understand the hybridisation type %s\n",&ps_input[2]);
	  usage();
	  exit(EXIT_FAILURE);
      }
      break;
  case 'I':       /* An input file is provided */
      if ( strlen(&ps_input[2]) != 0){
	  i_infile = TRUE;
	  if ( (pF_INFILE = fopen(&ps_input[2],"r")) == NULL){
	      fprintf(ERROR," I was not able to open the file %s\n",&ps_input[2]);
	      usage();
	      exit(EXIT_FAILURE);
	  } else {
	      /* HERE THERE IS A BUG: DOES NOT READ THE VERY LAST LINE */
	      while(!feof(pF_INFILE)){
		  ps_line = read_string(pF_INFILE);
		  if (ps_line == NULL || strlen(ps_line) == 0)
		      continue;
		  if ( ( ps_inputline = (char *)malloc(strlen(ps_line)+1) ) == 0){
		      fprintf(ERROR," function decode_input, line __LINE__:"
			      " Unable to allocate memory to get the name of inputfile\n");
		      exit(EXIT_FAILURE);
		  }
		  sscanf(ps_line,"%s",ps_inputline);
		  if (strncmp(ps_inputline,"-",1) == 0){
		      pst_in_param = decode_input(pst_in_param,ps_inputline,ps_path);
		  } else {
		      fprintf(ERROR," I did not understand this line of input file: %s\n",ps_inputline);
		      usage();
		      exit(EXIT_FAILURE);
		  }
		  free(ps_inputline);
	      }
	  }
      } else {
	  fprintf(ERROR," I did not understand the file %s\n",&ps_input[2]);
	  usage();
	  exit(EXIT_FAILURE);
      }
      break;
  case 'K':       /* Enter another correction for salt concentration */
      if (strncmp(&ps_input[2],"san96a",6) == 0
	  || strncmp(&ps_input[2],"san98a",6) == 0
	  || strncmp(&ps_input[2],"nak99a",6) == 0
	  || strncmp(&ps_input[2],"wet91a",6) == 0){
	  strncpy(pst_in_param->s_sodium_correction,&ps_input[2],6);
	  pst_in_param->s_sodium_correction[6]='\0';
      } else {
	  fprintf(ERROR," I did not understand your salt correction\n"
		  " Please read the manual to find the available corrections\n");
	  usage();
	  exit(EXIT_FAILURE);
      }
      break;
  case 'L':       /* please give me the legal notice */
      legal();
      exit(EXIT_SUCCESS);
  case 'M':       /* Alternative Nearest-neighbor set for mismatches */
      if ( strlen(&ps_input[2]) != 0 && isalnum((int)ps_input[2]) ){
	  i_alt_mm = TRUE;
	  if (pst_in_param->pst_present_mm != NULL)
	      pst_in_param->pst_present_mm = NULL; /* Reset the NN set */
	  pst_in_param->pst_present_mm = read_mismatches(&ps_input[2],ps_path);
	  strncpy(pst_in_param->pst_present_mm->s_mmfile,&ps_input[2],FILE_MAX); 
	  pst_in_param->pst_present_mm->s_mmfile[FILE_MAX-1] = '\0'; /* security check */
      } else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
      }
    break;
  case 'i':       /* Alternative Nearest-neighbor set for inosine base pairs */
      if ( strlen(&ps_input[2]) != 0 && isalnum((int)ps_input[2]) ){
	  i_alt_inosine = TRUE;
	  i_inosineneed = TRUE;
	  if (pst_in_param->pst_present_inosine != NULL)
	      pst_in_param->pst_present_inosine = NULL; /* Reset the NN set */
	  pst_in_param->pst_present_inosine = read_inosine(&ps_input[2],ps_path);
	  strncpy(pst_in_param->pst_present_inosine->s_inosinefile,&ps_input[2],FILE_MAX); 
	  pst_in_param->pst_present_inosine->s_inosinefile[FILE_MAX-1] = '\0'; /* security check */
      } else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
      }
    break;
   case 'N':
      /* sodium concentration */
      if ( strlen(&ps_input[2]) != 0 && isdigit((int)ps_input[2]) ){
      pst_in_param->d_conc_salt = strtod(&ps_input[2],NULL);
      i_salt = TRUE;
      } else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
      }
      break;
  case 'k':
      /* potassium concentration */
      if ( strlen(&ps_input[2]) != 0 && isdigit((int)ps_input[2]) ){
      pst_in_param->d_conc_potassium = strtod(&ps_input[2],NULL);
      i_magnesium = TRUE;
      i_salt = TRUE;
      } else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
      }
      break;
  case 't':
      /* tris concentration */
      if ( strlen(&ps_input[2]) != 0 && isdigit((int)ps_input[2]) ){
      pst_in_param->d_conc_tris = strtod(&ps_input[2],NULL);
      i_magnesium = TRUE;
      i_salt = TRUE;
      } else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
      }
      break;
  case 'G':
      /* magnesium concentration */
      if ( strlen(&ps_input[2]) != 0 && isdigit((int)ps_input[2]) ){
      pst_in_param->d_conc_magnesium = strtod(&ps_input[2],NULL);
      i_magnesium = TRUE;
      i_salt = TRUE;
      } else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
      }
      break;
  case 'O':
      /* An output file is required */
      i_outfile = TRUE;
      if (strlen(&ps_input[2]) != 0 && strlen(&ps_input[2]) <= FILE_MAX)
	  strncpy(pst_in_param->s_outfile,&ps_input[2],FILE_MAX);
      else {
	  if (strlen(&ps_input[2]) == 0){
	      /*+-----------------------------------------------+ 
		| compute the a code based on the year/day/time |
		| to complement the file names                  |
		+-----------------------------------------------+*/
	      time(&universal_time);
	      now = localtime(&universal_time);
	      i_year = 1900 + now->tm_year;
	      i_month = now->tm_mon;
	      switch(i_month){
	      case  0: strncpy(s_nmonth,"JAN",sizeof(s_nmonth)); break;
	      case  1: strncpy(s_nmonth,"FEB",sizeof(s_nmonth)); break;
	      case  2: strncpy(s_nmonth,"MAR",sizeof(s_nmonth)); break;
	      case  3: strncpy(s_nmonth,"APR",sizeof(s_nmonth)); break;
	      case  4: strncpy(s_nmonth,"MAY",sizeof(s_nmonth)); break;
	      case  5: strncpy(s_nmonth,"JUN",sizeof(s_nmonth)); break;
	      case  6: strncpy(s_nmonth,"JUL",sizeof(s_nmonth)); break;
	      case  7: strncpy(s_nmonth,"AUG",sizeof(s_nmonth)); break;
	      case  8: strncpy(s_nmonth,"SEP",sizeof(s_nmonth)); break;
	      case  9: strncpy(s_nmonth,"OCT",sizeof(s_nmonth)); break;
	      case 10: strncpy(s_nmonth,"NOV",sizeof(s_nmonth)); break;
	      case 11: strncpy(s_nmonth,"DEC",sizeof(s_nmonth)); break;
	      default: break;
	      }
	      i_day = now->tm_mday;
	      i_hour = now->tm_hour;
	      i_min = now->tm_min;
	      sprintf(pst_in_param->s_outfile,"melting%d%s%d_%dh%dm.out",i_year,s_nmonth,i_day,i_hour,i_min);
	  } else {
	      fprintf(ERROR," I was not able to treat the option %s\n",ps_input);
	      usage();
	      exit(EXIT_FAILURE);
	  }
      }
      break;
  case 'P':
      /* concentration of the strand in excess (P for "probe") */
      if ( strlen(&ps_input[2]) != 0 && isdigit((int)ps_input[2]) ){
	  pst_in_param->d_conc_probe = strtod(&ps_input[2],NULL);
	  i_probe = TRUE;
      } else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
      }
      break;
  case 'p':
      /* displays the path where to look for the set of parameters and quit */
      fprintf(OUTPUT,"path: %s\n",ps_path);
      exit(EXIT_SUCCESS);
      break;
  case 'q':
      /* no interactive correction of parameters */
      if (i_quiet == FALSE) 
	  i_quiet = TRUE;
      else i_quiet = FALSE;
      break;
  case 'S':
      /* Sequence */
      if ( strlen(&ps_input[2]) != 0 ){
	  if ( ( pst_in_param->ps_sequence = (char *)malloc(strlen(&ps_input[2])+1) ) == NULL){
	      fprintf(ERROR," Function decode_input, line __LINE__:"
		      "Unable to allocate memory to register the sequence\n");
	      exit(EXIT_FAILURE);
	  }
	  strncpy(pst_in_param->ps_sequence,&ps_input[2],strlen(&ps_input[2])+1);
	  pst_in_param->ps_sequence[strlen(&ps_input[2])] = '\0'; /* security check */
	  i_seq = TRUE;
      } else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
      }
      break;
    case 'T':
      /* max length before approximative calculus */
      if ( strlen(&ps_input[2]) != 0 && isdigit((int)ps_input[2]) ){
	  i_threshold = strtol(&ps_input[2],NULL,0); /* humm, check length here, no? */
      } else {
	  fprintf(ERROR," I did not understand the option %s\n",ps_input);
	  usage();
	  exit(EXIT_FAILURE);
      }      
      break;
    case 'v':
    /* Verbose mode */
      if (i_verbose == FALSE) 
	  i_verbose = TRUE;
      else i_verbose = FALSE;
      break;
  case 'V':
      /* Displays version and quit */
      fprintf(OUTPUT,"Version: %3.1f\n",VERSION);
      exit(EXIT_SUCCESS);
  case 'x':
      /* Force approximative tm computation */
      i_approx = TRUE;
      break;
  default:
      fprintf(ERROR," I did not understand the option %s\n",ps_input);
      usage();
      exit(EXIT_FAILURE);
  }
  return pst_in_param;
}

/***********************************
 * read a file containing a nn set *
 ***********************************/

struct nnset *read_nn(char *ps_nn_set, char *ps_path){
    struct nnset *pst_current_nn; /* pointer on a structure containing a set of nn_param */
    FILE *pF_nn_file;		  /* handle of file containing a set of nn param */
    char s_line[MAX_LINE];	  /* contains a line of a file or of stdin */
    char *pc_line_ptr;		  /* pointer moving along an input line */
    char *ps_nn_path;		  /* contains the address of the nn file */
    int i_crickcount = 0;	  /* counter of recorded crick's pairs */
    int i_count;
  
       /*+-----------------------------------------------------+
         | initialise a structure containing the nn parameters |
         +-----------------------------------------------------+*/

    pst_current_nn = (struct nnset *)calloc(1,sizeof(struct nnset));
    /* I AM SUPPOSED TO FREE THAT INSIDE THIS FUNCTION */

    if ((ps_nn_path = (char *)malloc(FILE_MAX)) == NULL){
      fprintf(ERROR," Function decode_input, line __LINE__:"
	      " Unable to allocate memory for the path of the NN file.\n");
      exit(EXIT_FAILURE);
    }

    sprintf(ps_nn_path,"%s/%s",ps_path,ps_nn_set);
    /* construct the complete name of the nn set file */
    
    if ( (pF_nn_file = fopen(ps_nn_path,"r")) == NULL){
      /* cannot open file containing alternative nn set in this path */
      fprintf(ERROR," I was not able to open the file %s,\n"
	      " supposed to contain the set of nearest_neighbor parameters.\n",ps_nn_path);
      
      sprintf(ps_nn_path,"%s/%s",NN_BASE,ps_nn_set);
      
      if ( (pF_nn_file = fopen(ps_nn_path,"r")) == NULL){
	/* permits to check default directory if the env defined does not work */
	/* for instance -Hdnadna defined in default and -Afoo97a.nn define by env var*/
	fprintf(ERROR," I was not able to open the file %s,\n"
		" supposed to contain the set of nearest_neighbor parameters.\n",ps_nn_path);
	usage();
	exit(EXIT_FAILURE);
      }
      fprintf(ERROR," I am going to use the file %s instead.\n",ps_nn_path);
    }

       /*+-------------------------------------+
         | read the file containing the nn set |
         +-------------------------------------+*/
    i_count = 0;		/* initialise reference counter */
    while(!feof(pF_nn_file)){	/* read the file nn_file until it ends */
	fgets(s_line,sizeof(s_line),pF_nn_file);	/* read a line the file */
	pc_line_ptr = s_line;
	if(*pc_line_ptr == ' ')	/* skip uninformative spaces */
	    pc_line_ptr++;
	if(*pc_line_ptr == '/' || *pc_line_ptr == '\n')	/* skip empty line */
	    continue;
	if(*pc_line_ptr == 'R'){
	    if (i_count >= NUM_REF) /* In case there are too many references */
		continue;
/*	    pst_current_nn->s_reference[i_count] = (char *)malloc(MAX_REF); */
/* FIXME: has to allocate place for a new reference */
	    strncpy(pst_current_nn->s_reference[i_count],s_line,MAX_REF); /* enter the reference */
	    pst_current_nn->s_reference[i_count][MAX_REF - 1] = '\0'; /* security lock */
	    i_count++;		/* ready for next reference */
	}
	if(*pc_line_ptr == 'A'||*pc_line_ptr == 'a'
	   ||*pc_line_ptr == 'G'||*pc_line_ptr == 'g'
	   ||*pc_line_ptr == 'C'||*pc_line_ptr == 'c'
	   ||*pc_line_ptr == 'T'||*pc_line_ptr == 't'
	   ||*pc_line_ptr == 'U'||*pc_line_ptr == 'u'
	   ||*pc_line_ptr == 'I'||*pc_line_ptr == 'i'){
	    if (i_crickcount <= NBNN){
		sscanf(s_line,"%3s %lf %lf",pst_current_nn->ast_nndata[i_crickcount].s_crick_pair,
		       &(pst_current_nn->ast_nndata[i_crickcount].d_enthalpy),
		       &(pst_current_nn->ast_nndata[i_crickcount].d_entropy));
		i_crickcount++;
	    }else {
		fprintf(ERROR," I detected too many Crick's pairs in that file.\n"
			      " Only 16 pairs and two initiation factors are allowed.\n");
		exit(EXIT_FAILURE);
	    }
	}
    }
    fclose(pF_nn_file);
    return pst_current_nn;
}

struct mmset *read_mismatches(char *ps_mm_set, char *ps_path){
    struct mmset *pst_current_mm; /* pointer on a structure containing a set of mismatches NN param */
    FILE *pF_mm_file;		  /* handle of file containing a set of mismatches NN param */
    char s_line[MAX_LINE];	  /* contains a line of a file or of stdin */
    char *pc_line_ptr;		  /* pointer moving along an input line */
    int i_crickcount = 0;	  /* counter of recorded crick's pairs */
    char *ps_mm_path;		  /* contains the address of the mismatches NN file */
    int i_count;

       /*+-----------------------------------------------------------------+
         | initialise a structure containing the parameters for mismatches |
         +-----------------------------------------------------------------+*/

    pst_current_mm = (struct mmset *)calloc(1,sizeof(struct mmset));
    /* I AM SUPPOSED TO FREE THAT INSIDE THIS FUNCTION */
  
       /*+-------------------------------------+
         | read the file containing the mm set |
         +-------------------------------------+*/

    if ((ps_mm_path = (char *)malloc(FILE_MAX)) == NULL){
	fprintf(ERROR," function read_mismatches, line __LINE__:\n"
		" Memory allocation error.\n");
	exit(EXIT_FAILURE);
    }
    sprintf(ps_mm_path,"%s/%s",ps_path,ps_mm_set);
    /* construct the complete name of the mm set file */

    if ( (pF_mm_file = fopen(ps_mm_path,"r")) == NULL){
    /* cannot open file containing alternative mm set in this path */

      fprintf(ERROR," I was not able to open the file %s,\n"
	      " supposed to contain the set of parameters for mismatches.\n",ps_mm_path);
      
      sprintf(ps_mm_path,"%s/%s",NN_BASE,ps_mm_set);

      if ( (pF_mm_file = fopen(ps_mm_path,"r")) == NULL){
	/* permits to check default directory if the env defined does not work */
	/* for instance -Hdnadna defined in default and -Afoo97a.nn define by env var*/
	fprintf(ERROR," I was not able to open the file %s,\n"
		" supposed to contain the set of parameters for mismatches.\n",ps_mm_path);
	usage();
	exit(EXIT_FAILURE);
      }
      fprintf(ERROR," I am going to use the file %s instead.\n",ps_mm_path);
    }
    i_count = 0;		/* initialise reference counter */
    while(!feof(pF_mm_file)){	/* read the file mm_file until it ends */
	fgets(s_line,sizeof(s_line),pF_mm_file);	/* read a line the file */
	pc_line_ptr = s_line;
	if(*pc_line_ptr == ' ')	/* skip uninformative spaces */
	    pc_line_ptr++;
	if(*pc_line_ptr == '/' || *pc_line_ptr == '\n')	/* skip empty line*/
	    continue;
	if(*pc_line_ptr == 'R'){
	    if (i_count >= NUM_REF)       /* In case there are too many references */
		continue;
	    /*	    pst_current_mm->s_reference[i_count] = (char *)malloc(MAX_REF); */
	    /* FIXME: has to allocate place for a new reference */
	    strncpy(pst_current_mm->s_reference[i_count],s_line,MAX_REF); /* enter the reference */
	    pst_current_mm->s_reference[i_count][MAX_REF - 1] = '\0'; /* security lock */
	    i_count++;		/* ready for next reference */
	}
	if(*pc_line_ptr == 'A'||*pc_line_ptr == 'a'
	   ||*pc_line_ptr == 'G'||*pc_line_ptr == 'g'
	   ||*pc_line_ptr == 'C'||*pc_line_ptr == 'c'
	   ||*pc_line_ptr == 'T'||*pc_line_ptr == 't'
	   ||*pc_line_ptr == 'U'||*pc_line_ptr == 'u'
	   ||*pc_line_ptr == 'I'||*pc_line_ptr == 'i'){
	    if (i_crickcount <= NBMM){
	      sscanf(s_line,"%6s %lf %lf",pst_current_mm->ast_mmdata[i_crickcount].s_crick_pair,
		     &(pst_current_mm->ast_mmdata[i_crickcount].d_enthalpy),
		     &(pst_current_mm->ast_mmdata[i_crickcount].d_entropy));
	      i_crickcount++;
	    } else {
		fprintf(ERROR," I detected too many Crick's pairs in that file.\n"
			      " Only %d mismatch pairs are allowed.\n",NBMM);
		exit(EXIT_FAILURE);
	    }
	}
    }
    fclose(pF_mm_file);
    return pst_current_mm;
}

struct inosineset *read_inosine(char *ps_inosine_set, char *ps_path){
    struct inosineset *pst_current_inosine; /* pointer on a structure containing a set of inosine mismatches NN param */
    FILE *pF_inosine_file;		  /* handle of file containing a set of inosine mismatches NN param */
    char s_line[MAX_LINE];	  /* contains a line of a file or of stdin */
    char *pc_line_ptr;		  /* pointer moving along an input line */
    int i_crickcount = 0;	  /* counter of recorded crick's pairs */
    char *ps_inosine_path;		  /* contains the address of the inosine mismatches NN file */
    int i_count;
    
           /*+-----------------------------------------------------------------+
         | initialise a structure containing the parameters for inosine base pairs |
         +-----------------------------------------------------------------+*/

    pst_current_inosine = (struct inosineset *)calloc(1,sizeof(struct inosineset));
    /* I AM SUPPOSED TO FREE THAT INSIDE THIS FUNCTION */
  
       /*+-------------------------------------+
         | read the file containing the inosine set |
         +-------------------------------------+*/

    if ((ps_inosine_path = (char *)malloc(FILE_MAX)) == NULL){
	fprintf(ERROR," function read_inosine, line __LINE__:\n"
		" Memory allocation error.\n");
	exit(EXIT_FAILURE);
    }
    sprintf(ps_inosine_path,"%s/%s",ps_path,ps_inosine_set);
    /* construct the complete name of the inosine set file */

    if ( (pF_inosine_file = fopen(ps_inosine_path,"r")) == NULL){
    /* cannot open file containing alternative inosine set in this path */

      fprintf(ERROR," I was not able to open the file %s,\n"
	      " supposed to contain the set of parameters for inosine mismatches.\n",ps_inosine_path);
      
      sprintf(ps_inosine_path,"%s/%s",NN_BASE,ps_inosine_set);

      if ( (pF_inosine_file = fopen(ps_inosine_path,"r")) == NULL){
	/* permits to check default directory if the env defined does not work */
	/* for instance -Hdnadna defined in default and -Afoo97a.nn define by env var*/
	fprintf(ERROR," I was not able to open the file %s,\n"
		" supposed to contain the set of parameters for inosine mismatches.\n",ps_inosine_path);
	usage();
	exit(EXIT_FAILURE);
      }
      fprintf(ERROR," I am going to use the file %s instead.\n",ps_inosine_path);
    }
    i_count = 0;		/* initialise reference counter */
    while(!feof(pF_inosine_file)){	/* read the file inosine_file until it ends */
	fgets(s_line,sizeof(s_line),pF_inosine_file);	/* read a line the file */
	pc_line_ptr = s_line;
	if(*pc_line_ptr == ' ')	/* skip uninformative spaces */
	    pc_line_ptr++;
	if(*pc_line_ptr == '/' || *pc_line_ptr == '\n')	/* skip empty line*/
	    continue;
	if(*pc_line_ptr == 'R'){
	    if (i_count >= NUM_REF)       /* In case there are too many references */
		continue;
	    /*	    pst_current_inosine->s_reference[i_count] = (char *)malloc(MAX_REF); */
	    /* FIXME: has to allocate place for a new reference */
	    strncpy(pst_current_inosine->s_reference[i_count],s_line,MAX_REF); /* enter the reference */
	    pst_current_inosine->s_reference[i_count][MAX_REF - 1] = '\0'; /* security lock */
	    i_count++;		/* ready for next reference */
	}
	if(*pc_line_ptr == 'A'||*pc_line_ptr == 'a'
	   ||*pc_line_ptr == 'G'||*pc_line_ptr == 'g'
	   ||*pc_line_ptr == 'C'||*pc_line_ptr == 'c'
	   ||*pc_line_ptr == 'T'||*pc_line_ptr == 't'
	   ||*pc_line_ptr == 'U'||*pc_line_ptr == 'u'
	   ||*pc_line_ptr == 'I'||*pc_line_ptr == 'i'){
	    if (i_crickcount <= NBIN){
	      sscanf(s_line,"%6s %lf %lf",pst_current_inosine->ast_inosinedata[i_crickcount].s_crick_pair,
		     &(pst_current_inosine->ast_inosinedata[i_crickcount].d_enthalpy),
		     &(pst_current_inosine->ast_inosinedata[i_crickcount].d_entropy));
	      i_crickcount++;
	    } else {
		fprintf(ERROR," I detected too many Crick's pairs in that file.\n"
			      " Only %d inosine mismatch pairs are allowed.\n",NBIN);
		exit(EXIT_FAILURE);
	    }
	}
    }
    fclose(pF_inosine_file);
    return pst_current_inosine;
}

struct deset *read_dangends(char *ps_de_set, char *ps_path){
    struct deset *pst_current_de; /* pointer on a structure containing a set of dangling ends NN param */
    FILE *pF_de_file;		  /* handle of file containing a set of dangling ends  NN param */
    char s_line[MAX_LINE];	  /* contains a line of a file or of stdin */
    char *pc_line_ptr;		  /* pointer moving along an input line */
    int i_crickcount = 0;	  /* counter of recorded crick's pairs */
    char *ps_de_path;		  /* contains the address of the mismatches NN file */
    int i_count;



       /*+--------------------------------------------------------------------+
         | initialise a structure containing the parameters for dangling ends |
         +--------------------------------------------------------------------+*/

    pst_current_de = (struct deset *)calloc(1,sizeof(struct deset));
    /* I AM SUPPOSED TO FREE THAT INSIDE THIS FUNCTION */
  
       /*+-------------------------------------+
         | read the file containing the de set |
         +-------------------------------------+*/

    if ((ps_de_path = (char *)malloc(FILE_MAX)) == NULL){
	fprintf(ERROR," function read_dangends, line __LINE__:\n"
		" Memory allocation error.\n");
	exit(EXIT_FAILURE);
    }

    sprintf(ps_de_path,"%s/%s",ps_path,ps_de_set);
    /* construct the complete name of the de set file */

    if ( (pF_de_file = fopen(ps_de_path,"r")) == NULL){
    /* cannot open file containing alternative nn set in this path */

      fprintf(ERROR," I was not able to open the file %s,\n"
	      " supposed to contain the set of parameters for dangling ends.\n",ps_de_path);

      sprintf(ps_de_path,"%s/%s",NN_BASE,ps_de_set);

      if ( (pF_de_file = fopen(ps_de_path,"r")) == NULL){
	/* permits to check default directory if the env defined does not work */
	/* for instance -Hdnadna defined in default and -Afoo97a.nn define by env var*/
	fprintf(ERROR," I was not able to open the file %s,\n"
		" supposed to contain the set of parameters for dangling ends.\n",ps_de_path);
	usage();
	exit(EXIT_FAILURE);
      }
      fprintf(ERROR," I am going to use the file %s instead.\n",ps_de_path);
    }
    
    i_count = 0;		/* initialise reference counter */
    while(!feof(pF_de_file)){	/* read the file de_file until it ends */
	fgets(s_line,sizeof(s_line),pF_de_file);	/* read a line the file */
	pc_line_ptr = s_line;
	if(*pc_line_ptr == ' ')	/* skip uninformative spaces */
	    pc_line_ptr++;
	if(*pc_line_ptr == '/' || *pc_line_ptr == '\n')	/* skip empty line */
	    continue;
	if(*pc_line_ptr == 'R'){
	    if (i_count >= NUM_REF)       /* In case there are too many references */
		continue;
	    /*	    pst_current_mm->s_reference[i_count] = (char *)malloc(MAX_REF); *//* FIXME: has to allocate place for a new reference */
	    strncpy(pst_current_de->s_reference[i_count],s_line,MAX_REF); /* enter the reference */
	    pst_current_de->s_reference[i_count][MAX_REF - 1] = '\0'; /* security lock */
	    i_count++;		/* ready for next reference */
	}
	if(*pc_line_ptr == 'A'||*pc_line_ptr == 'a'
	   ||*pc_line_ptr == 'G'||*pc_line_ptr == 'g'
	   ||*pc_line_ptr == 'C'||*pc_line_ptr == 'c'
	   ||*pc_line_ptr == 'T'||*pc_line_ptr == 't'
	   ||*pc_line_ptr == 'U'||*pc_line_ptr == 'u'
	   ||*pc_line_ptr == 'I'||*pc_line_ptr == 'i'
	   ||*pc_line_ptr == '-' ){
	    if (i_crickcount <= NBDE){
		sscanf(s_line,"%6s %lf %lf",pst_current_de->ast_dedata[i_crickcount].s_crick_pair,
		       &(pst_current_de->ast_dedata[i_crickcount].d_enthalpy),
		       &(pst_current_de->ast_dedata[i_crickcount].d_entropy));
		i_crickcount++;
		printf("%s",pst_current_de->ast_dedata[i_crickcount].s_crick_pair);
	    }else {
		fprintf(ERROR," I detected too many Crick's pairs in that file.\n"
			      " Only %d dangling end pairs are allowed.\n",NBDE);
		exit(EXIT_FAILURE);
	    }
	}
    }
    fclose(pF_de_file);
    return pst_current_de;
}

/*****************************************************************************
 * This function has been adapted from the book of Jean-Pierre Braquelaire : *
 * << Methodologie de la programmation en C >> Masson (1998) pp200-201       *
 * ISBN : 2-225-83269-2                                                      *
 * It reads a string of variable size. If newline are present within the     *
 * string, they have to be preceded by backslash.                            * 
 *****************************************************************************/

char *read_string(FILE *stream){
    int i_size = 8;		 /* initial size of the buffer */
    char *ps_buffer;             /* reading buffer */
    char *ps_string;             /* contain the complete string, e.g. sequence */

    if ( (ps_buffer = (char *)malloc(i_size)) == NULL ){
	fprintf(ERROR," function read_string: \n"
                      " Unable to allocate memory for the reading_buffer\n");  
	exit(EXIT_FAILURE);
    }
    ps_string = ps_buffer;
    
    for(;;){
	int c = getc(stream);	/* reading a character */
	if (c == EOF){
	    free(ps_buffer);
	    return NULL;
	}
	*ps_string = c;
				/* newline */
	if (*ps_string == '\n'){
				/* backslash before */
	    if (ps_string > ps_buffer && *(ps_string - 1) == '\\'){
		ps_string--;
		continue;
	    }
	    break;		/* end of string to read */
	}
	if(++ps_string == ps_buffer + i_size){
				/* increase reading buffer size */
	    if( (ps_buffer = realloc(ps_buffer,2*i_size)) == NULL){
		fprintf(ERROR," function read_string: \n"
			      " Unable to re-allocate memory for the reading_buffer\n");  
		exit(EXIT_FAILURE);
	    }
	    ps_string = ps_buffer + i_size;
	    i_size*=2;		/* we increase the buffer geometrically */
	} 
    }
    *ps_string = '\0';
    if ( (ps_string = (char *)malloc(ps_string - ps_buffer + 1)) == NULL ){
	fprintf(ERROR," function read_string: \n"
		" Unable to allocate memory for the reading_buffer\n");  
	exit(EXIT_FAILURE);
    };
    /* I am scared by this strcpy. But the strncpy doesn't work. */
    strcpy(ps_string,ps_buffer);
    free(ps_buffer);
    return ps_string;
}

/**************************
 * Print the legal notice *
 **************************/

void legal(void){
    fprintf(OUTPUT,"   Melting is copyright (C) 1997, 2009 by Nicolas Le Novère\n\n");
    fprintf(OUTPUT,"   and Marine Dumousseau\n\n");
    fprintf(OUTPUT,"   This  program  is  free  software; you can redistribute it\n");
    fprintf(OUTPUT,"   and/or modify it under the terms of the GNU General Public\n");
    fprintf(OUTPUT,"   License  as  published  by  the  Free Software Foundation;\n");
    fprintf(OUTPUT,"   either version 2 of the License, or (at your  option)  any\n");
    fprintf(OUTPUT,"   later version.\n\n");
    fprintf(OUTPUT,"   This  program  is  distributed in the hope that it will be\n");
    fprintf(OUTPUT,"   useful, but WITHOUT ANY WARRANTY; without even the implied\n");
    fprintf(OUTPUT,"   warranty  of  MERCHANTABILITY  or FITNESS FOR A PARTICULAR\n");
    fprintf(OUTPUT,"   PURPOSE.  See the GNU  General  Public  License  for  more\n");
    fprintf(OUTPUT,"   details.\n\n");
    fprintf(OUTPUT,"   You  should have received a copy of the GNU General Public\n");
    fprintf(OUTPUT,"   License along with this program; if not, write to the Free\n");
    fprintf(OUTPUT,"   Software  Foundation,  Inc.,  59  Temple Place, Suite 330,\n");
    fprintf(OUTPUT,"   Boston, MA  02111-1307 USA\n\n");
    fprintf(OUTPUT,"   Nicolas Le Novère, Computational Biology, EMBL-EBI        \n");
    fprintf(OUTPUT,"   Hinxton CB10 1SD United-Kingdom                           \n");
    fprintf(OUTPUT,"   lenov@ebi.ac.uk\n");
}

