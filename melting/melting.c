/******************************************************************************
 *                               MELTING v4.3                                 *
 * This program   computes for a nucleotide probe, the enthalpy, the entropy  *
 * and the melting temperature of the binding to its complementary template.  *
 * Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  *
 *          Copyright (C) Nicolas Le Novère and Marine Dumousseau  1997-2013  *
 *                                                                            *
 * File: melting.c                                                            *
 * Date: 01/APR/2009                                                          *
 * Aim : program entry. Output results                                        *
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
      Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

      Nicolas Le Novère
      Babraham Institute, Babraham Research Campus
      Babraham CB22 3AT Cambridge United-Kingdom.
      n.lenovere@gmail.com
       
      Marine Dumousseau
      EMBL-EBI, Wellcome-Trust Genome Campus
      Hinxton CB10 1SD Cambridge United-Kingdom. 
      marine@ebi.ac.uk  
      
*/

/*-----------------------------------------------------------------------*
 | Priority for options is depend on the order of argument. On the       |
 | command line, in case of conflict, the last defined win. The          |
 | interactive entry exists only to correct errors of command line and   | 
 | infiles                                                               |
 |                                                                       |
 | Command line arguments:                                               |
 |        -A[Alternative NN set]                                         |
 |        -C[Complement]                                                 |
 |        -D[Alternative Dangling ends NN set]                           |
 |        -F[Factor to correct the concentration of nucleic acid]        |
 |        -G[magnesium]                                                  |
 |        -h     displays Help                                           |
 |        -H[Hybridation type]                                           |
 |        -I[Infile]                                                     |
 |        -i[Alternative inosine set]                                    |
 |        -K[salt Korrection]                                            |
 |        -k[potassium]                                                  |
 |        -L     displays Legal information                              |
 |        -M[Alternative Mismaches NN set]                               |
 |        -N[salt (N states for Na)]                                     |
 |        -G[magnesium]                                                  |
 |        -O[Outfile] (the name can be omitted)                          |
 |        -P[concentration of the strand in excess (P states for Probe)] |
 |        -p     displays the path where to seek the parameters and quit |
 |        -q     Quiet. Switch off interactive correction of parameters  |
 |        -S[Sequence]                                                   |
 |        -T[Threshold for approximative computation]                    |
 |        -t[tris]                                                       |
 |        -v     Verbose mode                                            |
 |        -V     displays Version and quit                               |
 |        -x     force approXimative calculus                            |
 |                                                                       |
 | here describe the structure of input file                             |
 |                                                                       |
 | IF SEQUENCE < i_threshold NT, THE CALCUL IS EXACT (METHOD OF NEAREST- |
 | NEIGHBORS). IF > i_threshold, THE CALCUL IS AN APPROXIMATION, BASED   |
 | ON %(G+C)                                                             |
 *-----------------------------------------------------------------------*/


/*>>>>>>>>>>>>>>>>>>>>>>>>>>>PREPROCESSOR INFORMATIONS<<<<<<<<<<<<<<<<<<<<<<<<*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "common.h"
#include "melting.h"

/*****************
 * main function *
 *****************/

int main(int argc, char *argv[]){
    
    int i_count;			/* loop counter */
    int i_seq_errors;		        /* used to count mistakes in sequence */
    char *pc_scan;		        /* scan a line */
    char c_answer;		        /* single-letter answer */
    char *ps_inputstring;	        /* set of configuration strings */
    char s_line[MAX_LINE];		/* Just to read a small line of input */
    struct param *pst_param;	        /* contains the parameters of the current run */
    struct thermodynamic *pst_results;  /* contains the results of the computation */
    char *ps_getenv;	 	        /* content of the NN_PATH variable */
    FILE *OUTFILE;

  	/* THE HANDLING OF  *pst_present_nn IS COMPLETELY SILLY */
        /* HAS TO BE ALLOCATED HERE RATHER THAN IN READ_NN */

    /*-----------------------------------------*
     | Initialisation of a parameter structure |
     *-----------------------------------------*/
        
    if ( (pst_param = (struct param *)malloc(sizeof(struct param))) == NULL){
	 fprintf(ERROR," Function main, line __LINE__:\n"
		 " Unable to allocate memory for the parameter structure\n");
	 return EXIT_FAILURE;
    }
    pst_param->d_conc_salt = 0.0;
    pst_param->d_conc_potassium = 0.0;
    pst_param->d_conc_tris = 0.0;
    pst_param->d_conc_magnesium = 0.0;
    pst_param->d_conc_probe = 0.0;       
/* Note however that the probe concentration has to be > 0, because of a logarithm, and 
   because an hybridation without nucleic acid isn't that much interresting ...*/
    if ( (pst_param->ps_sequence = (char *)malloc(1)) == NULL){
	fprintf(ERROR," Function main, line __LINE__:\n"
		" Unable to allocate memory for the sequence\n");
	return EXIT_FAILURE;
    }
    if ( (pst_param->ps_complement = (char *)malloc(1)) == NULL){
	fprintf(ERROR," Function main, line __LINE__:\n"
		" Unable to allocate memory for the complement\n");
	return EXIT_FAILURE;
    }
    pst_param->ps_sequence[0] = '\0';
    pst_param->ps_complement[0] = '\0';
    pst_param->d_gnat = DEFAULT_NUC_CORR;
    /* the following three lines are necessary under Win32 */
    pst_param->pst_present_nn = NULL;
    pst_param->pst_present_mm = NULL;
    pst_param->pst_present_inosine = NULL;
    pst_param->pst_present_de = NULL;
    strncpy(pst_param->s_sodium_correction, DEFAULT_SALT_CORR,sizeof(pst_param->s_sodium_correction));
    pst_param->s_sodium_correction[sizeof(pst_param->s_sodium_correction)]='\0'; 
    /* the length of the correction has to be only 6 characters + eos */

       /*+----------------------------+
         | check the path of nn files |
         +----------------------------+*/

    /* read the environment variable specifying the repository directory */
    if ( (ps_getenv = getenv("NN_PATH")) == NULL){ 
	
	ps_getenv = NN_BASE;
    }

     /*-------------------------------------*
      | sequential reading of the arguments |
      *-------------------------------------*/
/* I copy the arguments because I read in fr.comp.lang.c that argv 
   and argc are only guaranteed only within the main function. */

    for (i_count = 1; i_count < argc; i_count++){	
	if ( (ps_inputstring = (char *)malloc(strlen(argv[i_count])+1)) == NULL){
	fprintf(ERROR," Function main, line __LINE__:\n"
		" Unable to allocate memory for the parsing of the %dth argument\n",i_count);
	return EXIT_FAILURE;
    }
	strncpy(ps_inputstring,argv[i_count],strlen(argv[i_count]));
	ps_inputstring[strlen(argv[i_count])] = '\0';
	pst_param = decode_input(pst_param,ps_inputstring,ps_getenv);
	free(ps_inputstring);
    }

/* All the following is redundant. Recode to call decode_input with the adequat
   argument. Maybe separate parsing of arguments from fullfilling the
   instructions: arguments -> parsing -> call read sequence, read salt etc.
                 STDIN                -> call the adequate function
 */

    /*-----------------------------------*
     | The hybridation type is mandatory |
     *-----------------------------------*/
    while(i_hybridtype == FALSE){
	if (i_quiet == FALSE){
	    fprintf(MENU,"  No type of hybridation has been properly entered.\n"
		         "  The specification of this parameter is mandatory.\n"
		         "      [A]-default DNA/DNA\n"
		         "      [B]-default DNA/RNA\n"
		         "      [C]-default RNA/RNA\n"
		         "      [Q]-Quit the program\n");
	    fgets(s_line,sizeof(s_line),INPUT); 
/* FIXME: Try to see what is happening if the line 
   is over sizeof(s_line). Maybe suck the remaining with while( getchar() != EOF) */
	    pc_scan = s_line;
	    while(*pc_scan == ' ')
		pc_scan++;
	    c_answer = toupper((int)*pc_scan);
	    switch (c_answer){
		case 'A': 
		    i_hybridtype = TRUE; 
		    if (pst_param->pst_present_nn == NULL){
			pst_param->pst_present_nn = read_nn(DEFAULT_DNADNA_NN,ps_getenv);
			strncpy(pst_param->pst_present_nn->s_nnfile,DEFAULT_DNADNA_NN,FILE_MAX); 
		    }
		    i_dnadna = TRUE;
		    i_dnarna = FALSE;
		    i_rnarna = FALSE;
		    break;
		case 'B': 
		    i_hybridtype = TRUE; 
		    if (pst_param->pst_present_nn == NULL){
			pst_param->pst_present_nn = read_nn(DEFAULT_DNARNA_NN,ps_getenv);
			strncpy(pst_param->pst_present_nn->s_nnfile,DEFAULT_DNARNA_NN,FILE_MAX);
		    } 
		    i_dnadna = FALSE;
		    i_dnarna = TRUE;
		    i_rnarna = FALSE;
		    break;
		case 'C': 
		    i_hybridtype = TRUE; 
		    if (pst_param->pst_present_nn == NULL){
		    pst_param->pst_present_nn = read_nn(DEFAULT_RNARNA_NN,ps_getenv);
		    strncpy(pst_param->pst_present_nn->s_nnfile,DEFAULT_RNARNA_NN,FILE_MAX); 
		    }
		    i_dnadna = FALSE;
		    i_dnarna = FALSE;
		    i_rnarna = TRUE;
		    break;
		case 'Q': return EXIT_SUCCESS; 
		default: break; /* nothing */		    
	    }
	} else {
	    fprintf(ERROR," No proper type of hybridation has been entered.\n");
	    return EXIT_FAILURE;
	}
    }
    
    /*-------------------------------------*
     | The different ion concentrations are mandatory |
     *-------------------------------------*/
    if (pst_param->d_conc_salt < MIN_SALT || pst_param->d_conc_salt >= MAX_SALT )
	i_salt = FALSE;
	if (pst_param->d_conc_salt == 0 && i_magnesium == FALSE){
	i_salt = FALSE;
	}
    while(i_salt == FALSE){
	if (i_quiet == FALSE){
	    fprintf(MENU,"  No salt concentration has been properly entered.\n"
		         "  The specification of this parameter is mandatory.\n"
		         "  This concentration has to belong to ]%4.2f,%5.2f[\n"
		         "  Enter it now (Q to quit)                         \n",MIN_SALT,MAX_SALT);
	    fgets(s_line,sizeof(s_line),INPUT); /* Try to see what is happening if the line 
		      is over sizeof(ac_line). Maybe suck the remaining with while( getchar() != EOF) */
	    pc_scan = s_line;
	    while(*pc_scan == ' ')
		pc_scan++;
	    if (*pc_scan == 'Q' || *pc_scan == 'q')
		return EXIT_SUCCESS;
	    if ( isdigit((int)*pc_scan) ){
		pst_param->d_conc_salt = strtod(pc_scan,NULL);
		if (pst_param->d_conc_salt > MIN_SALT && pst_param->d_conc_salt < MAX_SALT)
		    i_salt = TRUE;
	    }
	} else{
	    fprintf(ERROR," No proper salt concentration has been entered.\n");
	    return EXIT_FAILURE;
	}
    }
    
   if (pst_param->d_conc_potassium < MIN_SALT)
	i_salt = FALSE;
    while(i_salt == FALSE){
	if (i_quiet == FALSE){
	    fprintf(MENU,"  No potassium concentration has been properly entered.\n"
		         "  The specification of this parameter is mandatory.\n"
		         "  This concentration has to be positive\n"
		         "  Enter it now (Q to quit)                         \n");
	    fgets(s_line,sizeof(s_line),INPUT); /* Try to see what is happening if the line 
		      is over sizeof(ac_line). Maybe suck the remaining with while( getchar() != EOF) */
	    pc_scan = s_line;
	    while(*pc_scan == ' ')
		pc_scan++;
	    if (*pc_scan == 'Q' || *pc_scan == 'q')
		return EXIT_SUCCESS;
	    if ( isdigit((int)*pc_scan) ){
		pst_param->d_conc_potassium = strtod(pc_scan,NULL);
		if (pst_param->d_conc_potassium >= MIN_SALT)
		    i_salt = TRUE;
	    }
	} else{
	    fprintf(ERROR," No proper potassium concentration has been entered.\n");
	    return EXIT_FAILURE;
	}
    }
    
    if (pst_param->d_conc_tris < MIN_SALT)
    	i_salt = FALSE;
    while(i_salt == FALSE){
	if (i_quiet == FALSE){
	    fprintf(MENU,"  No tris concentration has been properly entered.\n"
		         "  The specification of this parameter is mandatory.\n"
		         "  This concentration has to be positive\n"
		         "  Enter it now (Q to quit)                         \n");
	    fgets(s_line,sizeof(s_line),INPUT); /* Try to see what is happening if the line 
		      is over sizeof(ac_line). Maybe suck the remaining with while( getchar() != EOF) */
	    pc_scan = s_line;
	    while(*pc_scan == ' ')
		pc_scan++;
	    if (*pc_scan == 'Q' || *pc_scan == 'q')
		return EXIT_SUCCESS;
	    if ( isdigit((int)*pc_scan) ){
		pst_param->d_conc_tris = strtod(pc_scan,NULL);
		if (pst_param->d_conc_tris > MIN_SALT)
		    i_salt = TRUE;
	    }
	} else{
	    fprintf(ERROR," No proper tris concentration has been entered.\n");
	    return EXIT_FAILURE;
	}
    }
    
    if (pst_param->d_conc_magnesium < MIN_SALT)
    	i_salt = FALSE;
    while(i_salt == FALSE){
	if (i_quiet == FALSE){
	    fprintf(MENU,"  No magnesium concentration has been properly entered.\n"
		         "  The specification of this parameter is mandatory.\n"
		         "  This concentration has to be positive\n"
		         "  Enter it now (Q to quit)                         \n");
	    fgets(s_line,sizeof(s_line),INPUT); /* Try to see what is happening if the line 
		      is over sizeof(ac_line). Maybe suck the remaining with while( getchar() != EOF) */
	    pc_scan = s_line;
	    while(*pc_scan == ' ')
		pc_scan++;
	    if (*pc_scan == 'Q' || *pc_scan == 'q')
		return EXIT_SUCCESS;
	    if ( isdigit((int)*pc_scan) ){
		pst_param->d_conc_magnesium = strtod(pc_scan,NULL);
		if (pst_param->d_conc_tris >= MIN_SALT)
		    i_salt = TRUE;
	    }
	} else{
	    fprintf(ERROR," No proper ion concentration has been entered.\n");
	    return EXIT_FAILURE;
	}
     }
    if (i_approx == FALSE){ /* The approximative mode do not need the concentration of nucleic acid 
                               A good indication of how accurate it is ...*/
      /*-----------------------------------------------------------------*
	| The nucleic acid  concentration (strand in excess) is mandatory |
	*-----------------------------------------------------------------*/
      
      if (pst_param->d_conc_probe <= MIN_PROBE || pst_param->d_conc_probe >= MAX_PROBE )
	i_probe = FALSE;
      while(i_probe == FALSE){
	if (i_quiet == FALSE){
	  fprintf(MENU,"  No nucleic acid concentration has been properly entered.\n"
		  "  The specification of this parameter is mandatory.\n"
		  "  This concentration has to belong to ]%4.2f,%4.2f[\n"
		  "  Enter it now (Q to quit)                         \n",MIN_PROBE,MAX_PROBE);
	  fgets(s_line,sizeof(s_line),INPUT); /* Try to see what is happening if the line 
						 is over sizeof(ac_line). Maybe suck the remaining with while( getchar() != EOF) */
	  pc_scan = s_line;
	  while(*pc_scan == ' ')
	    pc_scan++;
	  if (*pc_scan == 'Q' || *pc_scan == 'q')
	    return EXIT_SUCCESS;
	  if ( isdigit((int)*pc_scan) ){
	    pst_param->d_conc_probe = strtod(pc_scan,NULL);
	    if (pst_param->d_conc_probe > MIN_PROBE && pst_param->d_conc_probe < MAX_PROBE )
	      i_probe = TRUE;
	  } 
	} else{
	  fprintf(ERROR," No proper nucleic acid concentration has been entered.\n");
	  return EXIT_FAILURE;
	}
      }
    }

    /*---------------------------*
     | The sequence is mandatory |
     *---------------------------*/

    if (i_seq == TRUE){
	if ( (i_seq_errors = check_sequence(pst_param->ps_sequence)) != 0)
	i_seq = FALSE;
    }
    while(i_seq == FALSE){
	if (i_quiet == FALSE){
	    fprintf(MENU,"  No nucleic acid sequence has been properly entered.\n"
		         "  The specification of this parameter is mandatory.\n"
		         "  Enter the sequence (Q to quit)\n"
                         "  (if there are newlines in sequence, precede them with a \n"
                         "  backslash \\)\n");
	    pst_param->ps_sequence = read_string(INPUT); /* read the sequence from INPUT */
	    if (pst_param->ps_sequence[0] == 'Q' || pst_param->ps_sequence[0] == 'q')
		return EXIT_SUCCESS;                   /* user wants to quit */
	    else {
		i_seq_errors = check_sequence(pst_param->ps_sequence);
		if ( i_seq_errors != 0 ){                  /*The sequence contains illegal characters*/
		    fprintf(ERROR," Your sequence contain %d non legal character(s)\n",i_seq_errors);
		} else i_seq = TRUE;	                   /* the sequence is acceptable */
	    }
	}else{
	    fprintf(ERROR," No proper sequence has been entered.\n");
	    return EXIT_FAILURE;
	}
    }
    
    /*------------------*
     | Check complement |
     *------------------*/
    
    if (i_complement == TRUE){
	if ( (i_seq_errors = check_sequence(pst_param->ps_complement)) != 0
	     || strlen(pst_param->ps_complement) != strlen(pst_param->ps_sequence) ){
	    /* The complement contains illegal characters or has not the right size */
	    fprintf(ERROR," Your complement contain illegal characters or has not \n"
                          " same length than the sequence\n");
	    i_complement = FALSE;
	    while(i_complement == FALSE){
		if (i_quiet == FALSE){
		    fprintf(MENU,"  Enter the sequence of the complement now (Q to quit)\n"
			         "  (if there are newlines in sequence, precede them with a \n"
                                 "  backslash \\)");
		    pst_param->ps_complement=read_string(INPUT); /* read the sequence from INPUT */
		    if (pst_param->ps_complement[0] == 'Q' || pst_param->ps_complement[0] == 'q')
			return EXIT_SUCCESS;                   /* user wants to quit */
		    else if ( (i_seq_errors = check_sequence(pst_param->ps_complement)) != 0
			      || strlen(pst_param->ps_complement) != strlen(pst_param->ps_sequence) ){
			/* The complement contains illegal characters or has not the right size */
			fprintf(ERROR," Your complement contain illegal characters or has not \n"
				" same length than the sequence\n");
		    } else i_complement = TRUE;	                   /* the sequence is acceptable */
		} else{
		    fprintf(ERROR," No proper sequence complement has been entered.\n");
		    return EXIT_FAILURE;
		}
	    }		
	} 
    } else pst_param->ps_complement = make_complement(pst_param->ps_sequence);
    
    /*+------------------------------------------------------------+
      | If we need mismatches parameters but none were entered ... |
      +------------------------------------------------------------+*/
    if (i_mismatchesneed == TRUE && i_alt_mm == FALSE){
	if (pst_param->pst_present_mm != NULL)
	    pst_param->pst_present_mm = NULL; 
	if (i_dnadna){
	  pst_param->pst_present_mm = read_mismatches(DEFAULT_DNADNA_MISMATCHES,ps_getenv);
	  strncpy(pst_param->pst_present_mm->s_mmfile,DEFAULT_DNADNA_MISMATCHES,FILE_MAX); 
	} else if (i_dnarna){
	  pst_param->pst_present_mm = read_mismatches(DEFAULT_DNARNA_MISMATCHES,ps_getenv);
	  strncpy(pst_param->pst_present_mm->s_mmfile,DEFAULT_DNARNA_MISMATCHES,FILE_MAX); 
	} else if (i_rnarna){
	  pst_param->pst_present_mm = read_mismatches(DEFAULT_RNARNA_MISMATCHES,ps_getenv);
	  strncpy(pst_param->pst_present_mm->s_mmfile,DEFAULT_RNARNA_MISMATCHES,FILE_MAX); 
	}
/*	i_alt_mm = TRUE;*/
    }
    
        /*+------------------------------------------------------------+
      | If we need inosine mismatches parameters but none were entered ... |
      +------------------------------------------------------------+*/
    if (i_inosineneed == TRUE && i_alt_inosine == FALSE){
	if (pst_param->pst_present_inosine != NULL)
	    pst_param->pst_present_inosine = NULL; 
	if (i_dnadna){
	  pst_param->pst_present_inosine = read_inosine(DEFAULT_DNADNA_INOSINE_MISMATCHES,ps_getenv);
	  strncpy(pst_param->pst_present_inosine->s_inosinefile,DEFAULT_DNADNA_INOSINE_MISMATCHES,FILE_MAX); 
	} else if (i_dnarna){
	  pst_param->pst_present_inosine = read_inosine(DEFAULT_DNARNA_INOSINE_MISMATCHES,ps_getenv);
	  strncpy(pst_param->pst_present_inosine->s_inosinefile,DEFAULT_DNARNA_INOSINE_MISMATCHES,FILE_MAX); 
	} else if (i_rnarna){
	  pst_param->pst_present_inosine = read_inosine(DEFAULT_RNARNA_INOSINE_MISMATCHES,ps_getenv);
	  strncpy(pst_param->pst_present_inosine->s_inosinefile,DEFAULT_RNARNA_INOSINE_MISMATCHES,FILE_MAX); 
	}
/*	i_alt_inosine = TRUE;*/
    }

    /*+---------------------------------------------------------------+
      | If we need dangling ends parameters but none were entered ... |
      +---------------------------------------------------------------+*/
    if (i_dangendsneed == TRUE && i_alt_de == FALSE){
	if (pst_param->pst_present_de != NULL)
	    pst_param->pst_present_de = NULL; 
	if (i_dnadna){
	  pst_param->pst_present_de = read_dangends(DEFAULT_DNADNA_DANGENDS,ps_getenv);
	  strncpy(pst_param->pst_present_de->s_defile,DEFAULT_DNADNA_DANGENDS,FILE_MAX); 
	} else if (i_dnarna){
	  pst_param->pst_present_de = read_dangends(DEFAULT_DNARNA_DANGENDS,ps_getenv);
	  strncpy(pst_param->pst_present_de->s_defile,DEFAULT_DNARNA_DANGENDS,FILE_MAX); 
	} else if (i_rnarna){
	  pst_param->pst_present_de = read_dangends(DEFAULT_RNARNA_DANGENDS,ps_getenv);
	  strncpy(pst_param->pst_present_de->s_defile,DEFAULT_RNARNA_DANGENDS,FILE_MAX); 
	}
/*	i_alt_de = TRUE;*/
      }

    /*+-------------------------------------+
      | Let's launch the actual computation |
      +-------------------------------------+*/
    pst_results = get_results(pst_param);
    
    if (i_outfile == TRUE){	/* REDIRECTION IN OUTFILE */
	OUTFILE = fopen(pst_param->s_outfile,"w");
	/*+-----------------------------------------+
	  | printing verbose information in outfile |
	  +-----------------------------------------+*/
	if (i_verbose){
	    fprintf(OUTFILE,"\n");
	    fprintf(OUTFILE,"******************************************************************************\n");
	    fprintf(OUTFILE,"*                                MELTING %3.1f                                  *\n",VERSION);
            fprintf(OUTFILE,"* This program   computes for a nucleotide probe, the enthalpy, the entropy  *\n");
            fprintf(OUTFILE,"* and the melting temperature of the binding to its complementary template.  *\n");
            fprintf(OUTFILE,"* Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  *\n");
            fprintf(OUTFILE,"*    Copyright (C) Nicolas Le Novère and Marine Dumousseau 1997-2009        *\n");
            fprintf(OUTFILE,"******************************************************************************\n");
	    fprintf(OUTFILE,"\n");
	    fprintf(OUTFILE,"sequence  : %s\n",pst_param->ps_sequence);
	    fprintf(OUTFILE,"complement: %s\n",pst_param->ps_complement);
	    fprintf(OUTFILE,"\n");
	    if (i_dnarna == TRUE || i_rnarna == TRUE)
		fprintf(OUTFILE,"(Note that uridine is changed into thymidine for sake of simplification. The\n"
		                "computation has been nevertheless performed with the specified hybridisation\n"
                                "type. There is also not magnesium correction for dnarna and rnarna hybridization.)\n");
	    fprintf(OUTFILE,"Sodium concentration: %5.2e M\n",pst_param->d_conc_salt);
	    fprintf(OUTFILE,"Potassium concentration: %5.2e M\n",pst_param->d_conc_potassium);
	    fprintf(OUTFILE,"Tris concentration: %5.2e M\n",pst_param->d_conc_tris);
	    fprintf(OUTFILE,"Magnesium concentration: %5.2e M\n",pst_param->d_conc_magnesium);
	    fprintf(OUTFILE,"Nucleic acid concentration (strand in excess): %5.2e M\n",pst_param->d_conc_probe);
	    if (i_approx == FALSE){
		fprintf(OUTFILE,"File containing the nearest_neighbor parameters is %s.\n\n",pst_param->pst_present_nn->s_nnfile);
		for (i_count = 0; i_count < NUM_REF;i_count++){
		    if (pst_param->pst_present_nn->s_reference[i_count][0] == 'R')
			fprintf(OUTFILE,"%s",pst_param->pst_present_nn->s_reference[i_count]);
		    fprintf(OUTFILE,"\n");
		}
		fprintf(OUTFILE,"NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n"
			        "--------------------------------\n");
		for (i_count = 0; i_count < NBNN; i_count++)
		  fprintf(OUTFILE,"%s\t%8.1f\t%6.2f\n",pst_param->pst_present_nn->ast_nndata[i_count].s_crick_pair,
			  pst_param->pst_present_nn->ast_nndata[i_count].d_enthalpy * 4.18,
			  pst_param->pst_present_nn->ast_nndata[i_count].d_entropy * 4.18);
		
		if (i_mismatchesneed){
		    fprintf(OUTFILE,"File containing the nearest_neighbor parameters for mismatches is %s.\n\n",pst_param->pst_present_mm->s_mmfile);
		    for (i_count = 0; i_count < NUM_REF;i_count++){
			if (pst_param->pst_present_mm->s_reference[i_count][0] == 'R')
			    fprintf(OUTFILE,"%s",pst_param->pst_present_mm->s_reference[i_count]);
			fprintf(OUTFILE,"\n");
		    }
		    fprintf(OUTFILE,"NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n"
		                    "--------------------------------\n");
		    for (i_count = 0; i_count < NBMM; i_count++)
			if (strncmp(pst_param->pst_present_mm->ast_mmdata[i_count].s_crick_pair,"",2) != 0){
			    fprintf(OUTFILE,"%s\t%8.1f\t%6.2f\n",pst_param->pst_present_mm->ast_mmdata[i_count].s_crick_pair,
				    pst_param->pst_present_mm->ast_mmdata[i_count].d_enthalpy * 4.18,
				    pst_param->pst_present_mm->ast_mmdata[i_count].d_entropy * 4.18);
			}
		}
		if (i_inosineneed){
		    fprintf(OUTFILE,"File containing the nearest_neighbor parameters for inosine mismatches is %s.\n\n",pst_param->pst_present_inosine->s_inosinefile);
		    for (i_count = 0; i_count < NUM_REF;i_count++){
			if (pst_param->pst_present_inosine->s_reference[i_count][0] == 'R')
			    fprintf(OUTFILE,"%s",pst_param->pst_present_inosine->s_reference[i_count]);
			fprintf(OUTFILE,"\n");
		    }
		    fprintf(OUTFILE,"NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n"
		                    "--------------------------------\n");
		    for (i_count = 0; i_count < NBIN; i_count++)
			if (strncmp(pst_param->pst_present_inosine->ast_inosinedata[i_count].s_crick_pair,"",2) != 0){
			    fprintf(OUTFILE,"%s\t%8.1f\t%6.2f\n",pst_param->pst_present_inosine->ast_inosinedata[i_count].s_crick_pair,
				    pst_param->pst_present_inosine->ast_inosinedata[i_count].d_enthalpy * 4.18,
				    pst_param->pst_present_inosine->ast_inosinedata[i_count].d_entropy * 4.18);
			}
		}
		if (i_dangendsneed){
		    fprintf(OUTFILE,"File containing the nearest_neighbor parameters for dangling ends is %s.\n\n",pst_param->pst_present_de->s_defile);
		    for (i_count = 0; i_count < NUM_REF;i_count++){
			if (pst_param->pst_present_de->s_reference[i_count][0] == 'R')
			    fprintf(OUTFILE,"%s",pst_param->pst_present_de->s_reference[i_count]);
			fprintf(OUTFILE,"\n");
		    }
		    fprintf(OUTFILE,"NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n"
		                    "--------------------------------\n");
		    for (i_count = 0; i_count < NBDE; i_count++)
			if (strncmp(pst_param->pst_present_de->ast_dedata[i_count].s_crick_pair,"",2) != 0){
			    fprintf(OUTFILE,"%s\t%8.1f\t%6.2f\n",pst_param->pst_present_de->ast_dedata[i_count].s_crick_pair,
				    pst_param->pst_present_de->ast_dedata[i_count].d_enthalpy * 4.18,
				    pst_param->pst_present_de->ast_dedata[i_count].d_entropy * 4.18);
			}
		}
		
		if (i_magnesium == TRUE){
		    fprintf(OUTFILE,"\nThe monovalent ion correction and divalent ion correction is from owczarzy (2008), i.e,\n");
		    if (pst_param->d_conc_salt + pst_param->d_conc_potassium + pst_param->d_conc_tris/2 == 0) {
		    	fprintf(OUTFILE,"1/Tm(Mg2+) = 1/Tm(1M Na+) + a - b x ln([Mg2+]) + Fgc x (c + d x ln([Mg2+]) + 1/(2 x (Nbp - 1)) x (- e + f x ln([Mg2+]) + g x ln([Mg2+]) x\n"
			"ln([Mg2+]))\n"
			"where : a = 3.92/100000, b = 9.11/1000000, c = 6.26/100000,d = 1.42/100000,e = 4.82/10000;f = 5.25/10000, g = 8.31/100000.\n");
		    }
		    else {
		    	if (sqrt(pst_param->d_conc_magnesium)/(pst_param->d_conc_salt + pst_param->d_conc_potassium + pst_param->d_conc_tris/2) < 0.22) {
				fprintf(OUTFILE,"1/Tm(Mg2+) = 1/Tm(1M Na+) + (4.29 x Fgc - 3.95) x 1/100000 x ln([monovalent+]) + 9.40 x 1/1000000 x ln([monoalent+]) x \n"
				"ln([monovalent+])\n"
				"where : [Monovalent+] = [Na+] + [K+] + [tris+].\n");
		    	}
			else {
		    		if (sqrt(pst_param->d_conc_magnesium)/(pst_param->d_conc_salt + pst_param->d_conc_potassium + pst_param->d_conc_tris/2) < 6) {
					fprintf(OUTFILE,"1/Tm(Mg2+) = 1/Tm(1M Na+) + a - b x ln([Mg2+]) + Fgc x (c + d x ln([Mg2+]) + 1/(2 x (Nbp - 1)) x (- e + f x ln([Mg2+]) + g x ln([Mg2+]) x\n"
					"ln([Mg2+]))\n"
					"where : a = 3.92/100000 x (0.843 - 0.352 x [Monovalent+]0.5 x ln([Monovalent+])), b = 9.11/1000000, c = 6.26/100000, d = 1.42/100000 x\n"
					"(1.279 - 4.03/1000 x ln([monovalent+]) - 8.03/1000 x ln([monovalent+] x ln([monovalent+]),e = 4.82/10000;f = 5.25/10000, g = 8.31/100000\n"
					"x (0.486 - 0.258 x ln([monovalent+]) +\n"
					"5.25/1000 x ln([monovalent+] x\n");
	                                fprintf(OUTFILE,"ln([monovalent+] x ln([monovalent+]).\n"
				        "and [Monovalent+] = [Na+] + [K+] + [tris+]\n");
		    		}
				else {
					fprintf(OUTFILE,"1/Tm(Mg2+) = 1/Tm(1M Na+) + a - b x ln([Mg2+]) + Fgc x (c + d x ln([Mg2+]) + 1/(2 x (Nbp - 1)) x (- e + f x ln([Mg2+]) + g x ln([Mg2+]) x\n"
					"ln([Mg2+]))\n"
					"where : a = 3.92/100000, b = 9.11/1000000, c = 6.26/100000,d = 1.42/100000,e = 4.82/10000;f = 5.25/10000, g = 8.31/100000.\n");
				}
		    	}
		    }
		}
		else {
			if (strncmp(pst_param->s_sodium_correction,"wet91a",sizeof(pst_param->s_sodium_correction)) == 0){
		    		fprintf(OUTFILE,"\nThe salt correction is from Wetmur (1991), i.e,\n"
			            "16.6 x log([Na+] / (1 + 0.7 x [Na+])) + 3.85\n");
			}
			else if (strncmp(pst_param->s_sodium_correction,"san96a",sizeof(pst_param->s_sodium_correction)) == 0){
		    		fprintf(OUTFILE,"\nThe salt correction is from SantaLucia et al. (1996), i.e,\n"
                                    "12.5 x log[Na+]\n");
			}
			else if (strncmp(pst_param->s_sodium_correction,"san98a",sizeof(pst_param->s_sodium_correction)) == 0){
		    		fprintf(OUTFILE,"\nThe salt correction is from SantaLucia (1998), i.e,\n"
                                    "DeltaS = DeltaS([Na+]=1M) + 0.368 x (N-1) x ln[Na+]\n");
			}
		}
			
		fprintf(OUTFILE,"\nThe correction of the nucleic acid concentration is %3.1f,\n"
                                "i.e. the Tm for [Na+]=1M is DeltaH / [DeltaS + R x ln c/%3.1f]\n",pst_param->d_gnat,pst_param->d_gnat);
 
		fprintf(OUTFILE,"\nCrick's pairs contained in your sequence:\n");
		for(i_count = 0; i_count < NBNN; i_count++)
		    if (pst_results->i_crick[i_count] != 0)
			fprintf(OUTFILE,"%s\t%d\n",pst_param->pst_present_nn->ast_nndata[i_count].s_crick_pair,pst_results->i_crick[i_count]);	
		fprintf(OUTFILE,"\nMismatched pairs contained in your sequence:\n");
		for(i_count = 0; i_count< NBMM; i_count++)
		    if (pst_results->i_mismatch[i_count] != 0)
			fprintf(OUTFILE,"%s\t%d\n",pst_param->pst_present_mm->ast_mmdata[i_count].s_crick_pair,pst_results->i_mismatch[i_count]);
		fprintf(OUTFILE,"\nInosine mismatched pairs contained in your sequence:\n");
		for(i_count = 0; i_count< NBIN; i_count++)
		    if (pst_results->i_inosine[i_count] != 0)
			fprintf(OUTFILE,"%s\t%d\n",pst_param->pst_present_inosine->ast_inosinedata[i_count].s_crick_pair,pst_results->i_inosine[i_count]);		
		fprintf(OUTFILE,"\nDangling ends contained in your sequence:\n");
		for(i_count = 0; i_count< NBDE; i_count++)
		  if (pst_results->i_dangends[i_count] != 0)
			fprintf(OUTFILE,"%s\t%d\n",pst_param->pst_present_de->ast_dedata[i_count].s_crick_pair,pst_results->i_dangends[i_count]);		
		fprintf(OUTFILE,"\n");	
	    }
	}
      /*+------------------------------------+
        | print essential results in outfile |
        +------------------------------------+*/
	if (i_approx == FALSE){
	    fprintf(OUTFILE,"  Enthalpy: %7.0f J.mol-1\n", pst_results->d_total_enthalpy * 4.18);
	    fprintf(OUTFILE,"  Entropy: %7.2f J.mol-1.K-1\n", pst_results->d_total_entropy * 4.18);
	} else {
	  fprintf(OUTFILE,"  Sequence length above threshold: approximative mode\n");
	}
	fprintf(OUTFILE,"  Melting temperature: %5.2f deg C\n", pst_results->d_tm);
	fclose(OUTFILE);
    } else {			/* OUTPUT ON STDIN */
      /*+-------------------------------------------------+
        | printing verbose information BEFORE computation |
        +-------------------------------------------------+*/
	if (i_verbose){
	    fprintf(VERBOSE,"\n");
	    fprintf(VERBOSE,"******************************************************************************\n");
	    fprintf(VERBOSE,"*                                MELTING %3.1f                                  *\n",VERSION);
            fprintf(VERBOSE,"* This program   computes for a nucleotide probe, the enthalpy, the entropy  *\n");
            fprintf(VERBOSE,"* and the melting temperature of the binding to its complementary template.  *\n");
            fprintf(VERBOSE,"* Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  *\n");
            fprintf(VERBOSE,"*    Copyright (C) Nicolas Le Novère and Marine Dumousseau 1997-2009        *\n");
            fprintf(VERBOSE,"******************************************************************************\n");
	    fprintf(VERBOSE,"\n");
	    fprintf(VERBOSE,"sequence  : %s\n",pst_param->ps_sequence);
	    fprintf(VERBOSE,"complement: %s\n",pst_param->ps_complement);
	    fprintf(VERBOSE,"\n");
	    if (i_dnarna == TRUE || i_rnarna == TRUE)
		fprintf(VERBOSE,"(Note that uridine is changed into thymidine for sake of simplification. The\n"
                                "computation has been nevertheless performed with the specified hybridisation\n"
                                "type. There is also not magnesium correction for dnarna and rnarna hybridization.\n");
	    fprintf(VERBOSE,"Sodium concentration: %5.2e M\n",pst_param->d_conc_salt);
	    fprintf(VERBOSE,"Potassium concentration: %5.2e M\n",pst_param->d_conc_potassium);
	    fprintf(VERBOSE,"Tris concentration: %5.2e M\n",pst_param->d_conc_tris);
	    fprintf(VERBOSE,"Magnesium concentration: %5.2e M\n",pst_param->d_conc_magnesium);
	    fprintf(VERBOSE,"Nucleic acid concentration (strand in excess): %5.2e M\n",pst_param->d_conc_probe);
	    if (i_approx == FALSE){
		fprintf(VERBOSE,"File containing the nearest_neighbor parameters is %s.\n\n",pst_param->pst_present_nn->s_nnfile);
		for (i_count = 0; i_count < NUM_REF;i_count++){
		    if (pst_param->pst_present_nn->s_reference[i_count][0] == 'R')
			fprintf(VERBOSE,"%s",pst_param->pst_present_nn->s_reference[i_count]);
		}
		fprintf(VERBOSE,"\n");
		fprintf(VERBOSE,"NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n"
		                "-------------------------------\n");
		for (i_count = 0; i_count < NBNN; i_count++){
		    fprintf(VERBOSE,"%s\t%8.1f\t%6.2f\n",pst_param->pst_present_nn->ast_nndata[i_count].s_crick_pair,
			    pst_param->pst_present_nn->ast_nndata[i_count].d_enthalpy * 4.18,
			    pst_param->pst_present_nn->ast_nndata[i_count].d_entropy * 4.18);
		}
		if (i_mismatchesneed){
		    fprintf(VERBOSE,"File containing the nearest_neighbor parameters for mismatches is %s.\n\n",pst_param->pst_present_mm->s_mmfile);
		    for (i_count = 0; i_count < NUM_REF;i_count++){
			if (pst_param->pst_present_mm->s_reference[i_count][0] == 'R')
			    fprintf(VERBOSE,"%s",pst_param->pst_present_mm->s_reference[i_count]);
		    }
		    fprintf(VERBOSE,"\n");
		    fprintf(VERBOSE,"NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n"
		                    "-------------------------------\n");
		    for (i_count = 0; i_count < NBMM; i_count++)
			if ( pst_param->pst_present_mm->ast_mmdata[i_count].d_enthalpy != 99999){
			    fprintf(VERBOSE,"%s\t%8.1f\t%6.2f\n",pst_param->pst_present_mm->ast_mmdata[i_count].s_crick_pair,
				    pst_param->pst_present_mm->ast_mmdata[i_count].d_enthalpy * 4.18,
				    pst_param->pst_present_mm->ast_mmdata[i_count].d_entropy * 4.18);
			}
		}
		if (i_inosineneed){
		    fprintf(VERBOSE,"File containing the nearest_neighbor parameters for inosine mismatches is %s.\n\n",pst_param->pst_present_inosine->s_inosinefile);
		    for (i_count = 0; i_count < NUM_REF;i_count++){
			if (pst_param->pst_present_inosine->s_reference[i_count][0] == 'R')
			    fprintf(VERBOSE,"%s",pst_param->pst_present_inosine->s_reference[i_count]);
		    }
		    fprintf(VERBOSE,"\n");
		    fprintf(VERBOSE,"NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n"
		                    "-------------------------------\n");
		    for (i_count = 0; i_count < NBIN; i_count++)
			if ( pst_param->pst_present_inosine->ast_inosinedata[i_count].d_enthalpy != 99999){
			    fprintf(VERBOSE,"%s\t%8.1f\t%6.2f\n",pst_param->pst_present_inosine->ast_inosinedata[i_count].s_crick_pair,
				    pst_param->pst_present_inosine->ast_inosinedata[i_count].d_enthalpy * 4.18,
				    pst_param->pst_present_inosine->ast_inosinedata[i_count].d_entropy * 4.18);
			}
		}
		if (i_dangendsneed){
		    fprintf(VERBOSE,"File containing the nearest_neighbor parameters for dangling ends is %s.\n\n",pst_param->pst_present_de->s_defile);
		    for (i_count = 0; i_count < NUM_REF;i_count++){
			if (pst_param->pst_present_de->s_reference[i_count][0] == 'R')
			    fprintf(VERBOSE,"%s",pst_param->pst_present_de->s_reference[i_count]);
			fprintf(VERBOSE,"\n");
		    }
		    fprintf(VERBOSE,"NN\tenthalpy\tentropy\n\t(J.mol-1)\t(J.mol-1.K-1)\n"
		                    "--------------------------------\n");
		    for (i_count = 0; i_count < NBDE; i_count++)
			if (strncmp(pst_param->pst_present_de->ast_dedata[i_count].s_crick_pair,"",2) != 0){
			    fprintf(VERBOSE,"%s\t%8.1f\t%6.2f\n",pst_param->pst_present_de->ast_dedata[i_count].s_crick_pair,
				    pst_param->pst_present_de->ast_dedata[i_count].d_enthalpy * 4.18,
				    pst_param->pst_present_de->ast_dedata[i_count].d_entropy * 4.18);
			}
		}
		if (i_magnesium == TRUE){
		    fprintf(VERBOSE,"\nThe monovalent and bivalent ions correction is from Owczarzy (2008), i.e,\n");
		          	    
		   if (pst_param->d_conc_salt + pst_param->d_conc_potassium + pst_param->d_conc_tris/2 == 0) {
		    	fprintf(VERBOSE,"1/Tm(Mg2+) = 1/Tm(1M Na+) + a - b x ln([Mg2+]) + Fgc x (c + d x ln([Mg2+]) + 1/(2 x (Nbp - 1)) x (- e + f x ln([Mg2+]) + g x ln([Mg2+]) x\n"
			"ln([Mg2+]))\n"
			"where : a = 3.92/100000, b = 9.11/1000000, c = 6.26/100000,d = 1.42/100000,e = 4.82/10000;f = 5.25/10000, g = 8.31/100000.\n");
		    	}
		    	else {
		    		if (sqrt(pst_param->d_conc_magnesium)/(pst_param->d_conc_salt + pst_param->d_conc_potassium + pst_param->d_conc_tris/2) < 0.22) {
					fprintf(VERBOSE,"1/Tm(Mg2+) = 1/Tm(1M Na+) + (4.29 x Fgc - 3.95) x 1/100000 x ln([monovalent+]) + 9.40 x 1/1000000 x ln([monovalent+]) x\n"
					"ln([monovalent+])\n"
					"where : [Monovalent+] = [Na+] + [K+] + [Tris+].\n");
		    		}
				else {
		    			if (sqrt(pst_param->d_conc_magnesium)/(pst_param->d_conc_salt + pst_param->d_conc_potassium + pst_param->d_conc_tris/2) < 6) {
						fprintf(VERBOSE,"1/Tm(Mg2+) = 1/Tm(1M Na+) + a - b x ln([Mg2+]) + Fgc x (c + d x ln([Mg2+]) + 1/(2 x (Nbp - 1)) x (- e + f x ln([Mg2+]) + g x ln([Mg2+]) x\n"
						"ln([Mg2+]))\n");
						fprintf(VERBOSE,"where : a = 3.92/100000 x (0.843 - 0.352 x [Monovalent+]0.5 x ln([Monovalent+])), b = 9.11/1000000, c = 6.26/100000, d = 1.42/100000 x\n"
						"(1.279 - 4.03/1000 x ln([monovalent+]) - 8.03/1000 x ln([monovalent+] x ln([monovalent+]),e = 4.82/10000;f = 5.25/10000, g = 8.31/100000\n"
						"x (0.486 - 0.258 x ln([monovalent+]) + 5.25/1000 x ln([monovalent+] x ln([monovalent+] x ln([monovalent+]).\n"
						"and [Monovalent+] = [Na+] + [K+] + [Tris+]\n");
		    			}
					else {
						fprintf(VERBOSE,"1/Tm(Mg2+) = 1/Tm(1M Na+) + a - b x ln([Mg2+]) + Fgc x (c + d x ln([Mg2+]) + 1/(2 x (Nbp - 1)) x (- e + f x ln([Mg2+]) + g x ln([Mg2+]) x\n"
						"ln([Mg2+]))\n"
						"where : a = 3.92/100000, b = 9.11/1000000, c = 6.26/100000,d = 1.42/100000,e = 4.82/10000;f = 5.25/10000, g = 8.31/100000.\n");
					}
		    		}
		    	}
		}
		else {
			if (strncmp(pst_param->s_sodium_correction,"wet91a",sizeof(pst_param->s_sodium_correction)) == 0){
		    		fprintf(VERBOSE,"\nThe salt correction is from Wetmur (1991), i.e,\n"
			            "16.6 x log([Na+] / (1 + 0.7 x [Na+])) + 3.85\n");
			}
			else if (strncmp(pst_param->s_sodium_correction,"san96a",sizeof(pst_param->s_sodium_correction)) == 0){
		    		fprintf(VERBOSE,"\nThe salt correction is from SantaLucia et al. (1996), i.e,\n"
			            "12.5 x log[Na+]\n");
			}
			else if (strncmp(pst_param->s_sodium_correction,"san98a",sizeof(pst_param->s_sodium_correction)) == 0){
		    		fprintf(VERBOSE,"\nThe salt correction is from SantaLucia (1998), i.e,\n"
			            "DeltaS = DeltaS([Na+]=1M) + 0.368 x (N-1) x ln[Na+]\n");
			}
		}	
		fprintf(VERBOSE,"\nThe correction of the nucleic acid concentration is %3.1f,\n"
			        "i.e. the Tm for [Na+]=1M is DeltaH / [DeltaS + R x ln c/%3.1f]\n",pst_param->d_gnat,pst_param->d_gnat);
		fprintf(VERBOSE,"\nCrick's pairs contained in your sequence:\n");
		for(i_count = 0; i_count < NBNN; i_count++)
		    if (pst_results->i_crick[i_count] != 0)
			fprintf(VERBOSE,"%s\t%d\n",pst_param->pst_present_nn->ast_nndata[i_count].s_crick_pair,pst_results->i_crick[i_count]);	
		fprintf(VERBOSE,"\nMismatched pairs contained in your sequence:\n");
		for(i_count = 0; i_count< NBMM; i_count++)
		    if (pst_results->i_mismatch[i_count] != 0)
			fprintf(VERBOSE,"%s\t%d\n",pst_param->pst_present_mm->ast_mmdata[i_count].s_crick_pair,pst_results->i_mismatch[i_count]);
		fprintf(VERBOSE,"\nInosine mismatched pairs contained in your sequence:\n");
		for(i_count = 0; i_count< NBIN; i_count++)
		    if (pst_results->i_inosine[i_count] != 0)
			fprintf(VERBOSE,"%s\t%d\n",pst_param->pst_present_inosine->ast_inosinedata[i_count].s_crick_pair,pst_results->i_inosine[i_count]);
		fprintf(VERBOSE,"\nDangling ends contained in your sequence:\n");
		for(i_count = 0; i_count< NBDE; i_count++)
		  if (pst_results->i_dangends[i_count] != 0)
			fprintf(VERBOSE,"%s\t%d\n",pst_param->pst_present_de->ast_dedata[i_count].s_crick_pair,pst_results->i_dangends[i_count]);		
		fprintf(VERBOSE,"\n");	
	    }
	}
      /*+------------------------------------+
        |  print essential results on stdout |
        +------------------------------------+*/
	if (i_approx == FALSE){
	    fprintf(OUTPUT,"  Enthalpy: %7.0f J.mol-1\n", pst_results->d_total_enthalpy * 4.18);
	    fprintf(OUTPUT,"  Entropy: %7.2f J.mol-1.K-1\n", pst_results->d_total_entropy * 4.18);
	} else {
	  fprintf(OUTPUT,"  Sequence length above threshold: approximative mode\n");
	}

	fprintf(OUTPUT,"  Melting temperature: %5.2f deg C\n", pst_results->d_tm);
    }    	
			   /* This way of output the results is heavy and not satisfying ... */

    free(pst_param->ps_complement);
    free(pst_param->ps_sequence);
    free(pst_results);

    return EXIT_SUCCESS;
}



/**************************************
 * Precise the way to use the program *
 **************************************/

void usage(void){
    fprintf(OUTPUT,"  Usage is 'melting OPTIONS' where OPTIONS are:                        \n");
    fprintf(OUTPUT,"     -A[xxxxxx.nn]  Name of a file containing alternative nn parameters\n");
    fprintf(OUTPUT,"                    Defaults are: DNA/DNA: "DEFAULT_DNADNA_NN"         \n");
    fprintf(OUTPUT,"                                  DNA/RNA: "DEFAULT_DNARNA_NN"         \n");
    fprintf(OUTPUT,"                                  RNA/RNA: "DEFAULT_RNARNA_NN"         \n");
    fprintf(OUTPUT,"     -D[xxxxxx.nn]  Name of a file containing nn parameters for dangling ends\n");
    fprintf(OUTPUT,"                    Default is "DEFAULT_DNADNA_DANGENDS"             \n"); 
    fprintf(OUTPUT,"     -C[XXXXXXXXXX] Complementary sequence, mandatory if mismaches     \n");
    fprintf(OUTPUT,"     -F[x.xx]       Correction for the concentration of nucleic acid   \n");
    fprintf(OUTPUT,"                    Default is DEFAULT_NUC_CORR                       \n"); 
    fprintf(OUTPUT,"     -h             Displays this help and quit                        \n");
    fprintf(OUTPUT,"    -H[xxxxxx]     Type of hybridisation (exemple dnadna), mandatory  \n");
    fprintf(OUTPUT,"     -I[XXXXXX]     Name of an input file setting up the options       \n");
    fprintf(OUTPUT,"     -K             Salt correction. Default is "DEFAULT_SALT_CORR"    \n" );
    fprintf(OUTPUT,"    -L             Displays legal information and quit                \n");
    fprintf(OUTPUT,"     -M[xxxxxx.nn]  Name of a file containing nn parameters for mismatches\n");
    fprintf(OUTPUT,"                    Default is "DEFAULT_DNADNA_MISMATCHES"             \n");
    fprintf(OUTPUT,"     -i[xxxxxx.nn]  Name of a file containing nn parameters for inosine mismatches\n"); 
    fprintf(OUTPUT,"                    Defaults are: DNA/DNA: "DEFAULT_DNADNA_INOSINE_MISMATCHES"         \n");
    fprintf(OUTPUT,"                                  DNA/RNA: "DEFAULT_DNARNA_INOSINE_MISMATCHES"         \n");
    fprintf(OUTPUT,"                                  RNA/RNA: "DEFAULT_RNARNA_INOSINE_MISMATCHES"         \n");
    fprintf(OUTPUT,"     -N[x.xe-x]     Sodium concentration in mol.l-1. Mandatory         \n");
    fprintf(OUTPUT,"     -k[x.xe-x]     Potassium concentration in mol.l-1. Mandatory         \n");
    fprintf(OUTPUT,"     -t[x.xe-x]     Tris concentration in mol.l-1. The Tri+ concentration is about \n");
    fprintf(OUTPUT,"                    half of total Tris concentration Mandatory         \n");
    fprintf(OUTPUT,"     -G[x.xe-x]     Magnesium concentration in mol.l-1. Mandatory         \n");  
    fprintf(OUTPUT,"     -O[XXXXXX]     Name of an output file (the name can be omitted)   \n");
    fprintf(OUTPUT,"     -P[x.xe-x]     Concentration of single strand nucleic acid in mol.l-1. Mandatory\n");
    fprintf(OUTPUT,"     -p             Return path where to find the calorimetric tables\n");
    fprintf(OUTPUT,"     -q             Quiet. Switch off interactive correction of parameters\n");
    fprintf(OUTPUT,"     -S[XXXXXXXXXX] Nucleic acid sequence, mandatory                   \n");
    fprintf(OUTPUT,"     -T[XXX]        Threshold for approximative computation            \n");
    fprintf(OUTPUT,"     -v             Switch ON the verbose mode, issuing lot more info  \n");
    fprintf(OUTPUT,"                    (if already ON, switch if OFF). Default is OFF     \n");
    fprintf(OUTPUT,"     -V             Print the version number                           \n");
    fprintf(OUTPUT,"     -x             Force to compute an approximative tm               \n");
    fprintf(OUTPUT,"  More information is available in the user-guide. Type `man melting'  \n"
	           "  to access it, or consult one of the melting.xxx files, where xxx     \n"
                   "  states for lat1 (isolatin1 text), ps (postscript), pdf or html.\n");
}

/************************************
 * Check the legality of a sequence *
 ************************************/

int check_sequence(char *ps_sequence){
  char *pc_base;                /*moving pointer on base*/
  int i_mistakes = 0;           /*error counter*/

  pc_base = ps_sequence;

  while (*pc_base != '\0'){
    /* capitalise the bases */
      *pc_base = toupper((int)*pc_base);

				/* change uridine into thymidine */
      if (*pc_base == 'U')
	  *pc_base = 'T';
      if (*pc_base != 'A' && *pc_base != 'G' && *pc_base != 'C' && *pc_base != 'T' && *pc_base != '-' && *pc_base != 'I')
	  i_mistakes++;
      pc_base++;
  }
  return(i_mistakes);
}

/******************************************
 * Construct the complement of a sequence *
 ******************************************/

char *make_complement(char *ps_sequence){
  char *ps_complement;	       /* computed complement of the sequence input */
  char *pc_base_seq;           /* moving pointer on sequence */
  char *pc_base_comp;	       /* moving pointer on complement */
  
  if ( (ps_complement = (char *)malloc(strlen(ps_sequence)+1)) == NULL){
    fprintf(ERROR," function make_complement, line __LINE__:\n"
	    "Unable to allocate memory to register the sequence\n");
    exit(EXIT_FAILURE);
  }
  pc_base_seq = ps_sequence;	/* pointer on the beginning of sequence */
  pc_base_comp = ps_complement;
  while (*pc_base_seq != '\0'){
    switch (*pc_base_seq){
    case 'A':
      *pc_base_comp = 'T';
      break;
    case 'G':
      *pc_base_comp = 'C';
      break;
    case 'C':
      *pc_base_comp = 'G';
      break;
    case 'T':
      *pc_base_comp = 'A';
      break;
    case '-':
      *pc_base_comp = '-';
      break;
    case 'I':
       fprintf(ERROR," There are inosine base pairs in your sequence and no complementary sequence has been entered.\n");
       exit(EXIT_FAILURE);
       break;
	
     default:
      fprintf(ERROR," It seems that one base of sequence is illegal.\n"
	      " I cannot compute the complement.\n");
      exit(EXIT_FAILURE);
    }
    pc_base_seq++;
    pc_base_comp++;
  }
  *pc_base_comp = '\0';
  return ps_complement;
}





