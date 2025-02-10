

/******************************************************************************
 *                               MELTING v4.3                                 *
 * This program   computes for a nucleotide probe, the enthalpy, the entropy  *
 * and the melting temperature of the binding to its complementary template.  *
 * Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  *
 *          Copyright (C) Nicolas Le Novère and Marine Dumousseau  1997-2013  *
 *                                                                            *
 * File: calcul.c                                                             *
 * Date: 01/APR/2009                                                          *
 * Aim : Computes enthalpy, entropy and melting temperature.                  *
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
#include <math.h>
#include "common.h"
#include "calcul.h"

struct thermodynamic *get_results(struct param *pst_param){
    int i,j;			/* loop counters */
    int i_mismatch;             /* mismatche detector */
    int i_inosine;              /* inosine mismatche detector */
    int i_dangend;              /* dangling end detector */
    int i_length = 0;		/* length of the sequence */
    int i_proxoffset = 0;       /* offset due to dangling end on the proximal side*/
    int i_distoffset = 0;       /* offset due to dangling end on the distal side*/
    struct thermodynamic *pst_results; /* contains the results ... */

				/* initialisation of result variables */

    if ( (pst_results = (struct thermodynamic *)malloc(sizeof(struct thermodynamic))) == NULL){
      	 fprintf(ERROR," function get_results, line __LINE__:"
		 " Unable to allocate memory for the result structure\n");
	 exit(EXIT_FAILURE);
    }
    pst_results->d_total_enthalpy = 0.0;
    pst_results->d_total_entropy = 0.0;
    for ( i = 0; i < NBNN ; i++)
	pst_results->i_crick[i] = 0;
    for ( i = 0; i < NBMM ; i++)
	pst_results->i_mismatch[i] = 0;
    for ( i = 0; i < NBIN ; i++)
	pst_results->i_inosine[i] = 0;
    for ( i = 0; i < NBDE ; i++)
	pst_results->i_dangends[i] = 0;
    pst_results->d_tm = 0.0;

    /*+------------------------------------------------------------------+
      | The length is too important. approximative computation performed |
      +------------------------------------------------------------------+*/

    if (strlen(pst_param->ps_sequence) > i_threshold)
	i_approx = TRUE;

    if (i_approx == TRUE)
	pst_results->d_tm = tm_approx(pst_param);

    /*+------------------------------+
      | nearest-neighbor computation |
      +------------------------------+*/

    else {
	/* The algorithm of screening is heavy, not general enough and does not offer room for evolution. To be changed! */

      if ( *(pst_param->ps_sequence) == '-' || *(pst_param->ps_complement) == '-'){
	i_dangend = TRUE;
	if (i_dnadna == FALSE && i_alt_de == FALSE){
	  fprintf(OUTPUT,"  WARNING: The default dangling ends parameters can efficiently\n"
		      "  account only for the DNA/DNA hybridisation. You can enter an\n"
		      "  alternative set of parameters with the option -D\n");
	    }
	    i_proxoffset++;
	    for (i = 0 ; i < NBDE ; i++){ /* seek the dangling-end term */
	      if ( (strncmp(pst_param->ps_sequence,pst_param->pst_present_de->ast_dedata[i].s_crick_pair,2) == 0)
		   && (strncmp(pst_param->ps_complement,&pst_param->pst_present_de->ast_dedata[i].s_crick_pair[3],2) == 0)){
		pst_results->d_total_enthalpy += pst_param->pst_present_de->ast_dedata[i].d_enthalpy;
		pst_results->d_total_entropy += pst_param->pst_present_de->ast_dedata[i].d_entropy;
		pst_results->i_dangends[i]++;
		i_dangend = FALSE;
	      }
	    }
	    if (i_dangend == TRUE){
	      fprintf(ERROR," NN parameters for %c%c/%c%c not found. Check the file containing\n"
		      " the information on dangling ends\n",pst_param->ps_sequence[0],pst_param->ps_sequence[1],pst_param->ps_complement[0],pst_param->ps_complement[1]);
	      exit(EXIT_FAILURE);
	    }
	}	
	
	if ( *(pst_param->ps_sequence + strlen(pst_param->ps_sequence)-1) == '-' || *(pst_param->ps_complement + strlen(pst_param->ps_complement)-1) == '-'){
	  i_dangend = TRUE;
	  if (i_dnadna == FALSE && i_alt_de == FALSE){
	    fprintf(OUTPUT,"  WARNING: The default dangling ends parameters can efficiently\n"
		    "  account only for the DNA/DNA hybridisation. You can enter an\n"
		    "  alternative set of parameters with the option -D\n");
	  }
	  i_distoffset++;
	  for (i = 0 ; i < NBDE ; i++){ /* seek the dangling-end term */
	    if ( (strncmp(pst_param->ps_sequence+strlen(pst_param->ps_sequence)-2,pst_param->pst_present_de->ast_dedata[i].s_crick_pair,2) == 0)
		 && (strncmp(pst_param->ps_complement+strlen(pst_param->ps_sequence)-2,&pst_param->pst_present_de->ast_dedata[i].s_crick_pair[3],2) == 0)){
	      pst_results->d_total_enthalpy += pst_param->pst_present_de->ast_dedata[i].d_enthalpy;
	      pst_results->d_total_entropy += pst_param->pst_present_de->ast_dedata[i].d_entropy;
	      pst_results->i_dangends[i]++;
	      i_dangend = FALSE;
	    }
	  }
	  if (i_dangend == TRUE){
	    fprintf(ERROR," NN parameters for %c%c/%c%c not found. Check the file containing\n"
		    " the information on dangling ends\n",*(pst_param->ps_sequence+strlen(pst_param->ps_sequence)-2),*(pst_param->ps_sequence+strlen(pst_param->ps_sequence)-1),*(pst_param->ps_complement+strlen(pst_param->ps_sequence)-2),*(pst_param->ps_complement+strlen(pst_param->ps_sequence)-1));
	    exit(EXIT_FAILURE);
	  }
	}
				/* determination of initiation terms for proximal extremity*/
	if ( *(pst_param->ps_sequence + i_proxoffset) == 'A' || *(pst_param->ps_sequence + i_proxoffset) == 'T')
	    for(i = 0 ; i < NBNN ; i++) /* seek the initiation term */
		if ( strncmp(pst_param->pst_present_nn->ast_nndata[i].s_crick_pair,"IA",2) == 0){
		    pst_results->d_total_enthalpy += pst_param->pst_present_nn->ast_nndata[i].d_enthalpy;
		    pst_results->d_total_entropy += pst_param->pst_present_nn->ast_nndata[i].d_entropy;
		}
	if ( *(pst_param->ps_sequence + i_proxoffset) == 'G' || *(pst_param->ps_sequence + i_proxoffset) == 'C')
	    for(i = 0 ; i < NBNN ; i++) /* seek the initiation term */
		if ( strncmp(pst_param->pst_present_nn->ast_nndata[i].s_crick_pair,"IG",2) == 0){
		    pst_results->d_total_enthalpy += pst_param->pst_present_nn->ast_nndata[i].d_enthalpy;
		    pst_results->d_total_entropy += pst_param->pst_present_nn->ast_nndata[i].d_entropy;
		}
				/* determination of initiation terms for distal extremity*/
	if ( *(pst_param->ps_sequence + strlen(pst_param->ps_sequence)-1-i_distoffset) == 'A' || *(pst_param->ps_sequence + strlen(pst_param->ps_sequence)-1-i_distoffset) == 'T')
	    for(i = 0 ; i < NBNN ; i++) /* seek the initiation term */
		if ( strncmp(pst_param->pst_present_nn->ast_nndata[i].s_crick_pair,"IA",2) == 0){
		    pst_results->d_total_enthalpy += pst_param->pst_present_nn->ast_nndata[i].d_enthalpy;
		    pst_results->d_total_entropy += pst_param->pst_present_nn->ast_nndata[i].d_entropy;
		}
	if ( *(pst_param->ps_sequence + strlen(pst_param->ps_sequence)-1-i_distoffset) == 'G' || *(pst_param->ps_sequence + strlen(pst_param->ps_sequence)-1-i_distoffset) == 'C')
	    for(i = 0 ; i < NBNN ; i++) /* seek the initiation term */
		if ( strncmp(pst_param->pst_present_nn->ast_nndata[i].s_crick_pair,"IG",2) == 0){
		    pst_results->d_total_enthalpy += pst_param->pst_present_nn->ast_nndata[i].d_enthalpy;
		    pst_results->d_total_entropy += pst_param->pst_present_nn->ast_nndata[i].d_entropy;
		}
				/* Travel through the sequence */
	if (strlen(pst_param->ps_sequence)-1 <= 0){
	    fprintf(ERROR," Oups, the lengh of the sequence seems zero or less ...\n");
	    exit(1);
	}

	i_length = strlen(pst_param->ps_sequence)-1 -i_proxoffset -i_distoffset;
	for (i = 0 + i_proxoffset ; i < i_length; i++){
	    i_mismatch = FALSE;
	    i_inosine = FALSE;
	    if (pst_param->ps_complement[i] == 'I' || pst_param->ps_complement[i+1] == 'I'){
			i_inosine = TRUE;
			}
	    switch (pst_param->ps_sequence[i]) {
		case 'A':
		    if (pst_param->ps_complement[i] != 'T' && pst_param->ps_complement[i] != 'I')
			i_mismatch = TRUE;
		    break;
		case 'G':
		    if (pst_param->ps_complement[i] != 'C' && pst_param->ps_complement[i] != 'I')
			i_mismatch = TRUE;
		    break;
		case 'C':
		    if (pst_param->ps_complement[i] != 'G' && pst_param->ps_complement[i] != 'I')
			i_mismatch = TRUE;
		    break;
		case 'T':
		    if (pst_param->ps_complement[i] != 'A' && pst_param->ps_complement[i] != 'I')
			i_mismatch = TRUE;
		    break;
		case 'I':
		    i_inosine = TRUE;
		    break;
		default:
		    fprintf(ERROR," I do not recognize the base %c\n",pst_param->ps_sequence[i]);
		    exit(EXIT_FAILURE);
	    }
	    switch (pst_param->ps_sequence[i+1]) {
		case 'A':
		    if (pst_param->ps_complement[i+1] != 'T' && pst_param->ps_complement[i+1] != 'I')
			i_mismatch = TRUE;
		    break;
		case 'G':
		    if (pst_param->ps_complement[i+1] != 'C' && pst_param->ps_complement[i+1] != 'I')
			i_mismatch = TRUE;
		    break;
		case 'C':
		    if (pst_param->ps_complement[i+1] != 'G' && pst_param->ps_complement[i+1] != 'I')
			i_mismatch = TRUE;
		    break;
		case 'T':
		    if (pst_param->ps_complement[i+1] != 'A' && pst_param->ps_complement[i+1] != 'I'){
			i_mismatch = TRUE;
		    }
		    break;
		case 'I':
		    i_inosine = TRUE;
		    break;
		default:
		    fprintf(ERROR," I do not recognize the base %c\n",pst_param->ps_sequence[i+1]);
	    }
	    if (i_mismatch == TRUE || i_inosine == TRUE){
		if (i == (0 + i_proxoffset) || i == (i_length - 1)){
		    fprintf(ERROR," The effect of mismatches (inosine mismatches included) located on the two extreme positions\n"
                                  " of a duplex are unpredictable (i.e. each case has to be \n"
                                  " considered separately).\n");
		    exit(EXIT_FAILURE);
		}
		if (i_dnadna == FALSE && i_alt_mm == FALSE && i_mismatch == TRUE){
		  fprintf(OUTPUT,"  WARNING: The default mismatches parameters can efficiently\n"
			  "  account only for the DNA/DNA hybridisation. You can enter an\n"
			  "  alternative set of parameters with the option -M\n");
		}
		if (i_dnarna == TRUE && i_alt_inosine == FALSE){
		  fprintf(OUTPUT,"  WARNING: The default inosine mismatches parameters can efficiently\n"
			  "  account only for the DNA/DNA hybridisation or RNA/RNA hybridization (however not completed yet). You can enter an\n"
			  "  alternative set of parameters with the option -M\n");
		}
		if (i_rnarna == TRUE && i_alt_inosine == FALSE){
		  fprintf(OUTPUT,"  WARNING: The only default inosine mismatches parameters available\n"
			  "  are the I.U bas pairs. You can enter an\n"
			  "  alternative set of parameters with the option -i\n");
		}
		if (i_mismatch == TRUE){
				/*compare with each possible mismatched pair*/
		for (j = 0 ; j < NBMM ; j++){
		  if ( (strncmp(&pst_param->ps_sequence[i],pst_param->pst_present_mm->ast_mmdata[j].s_crick_pair,2) == 0)
		       && (strncmp(&pst_param->ps_complement[i],&pst_param->pst_present_mm->ast_mmdata[j].s_crick_pair[3],2) == 0)){
		    pst_results->i_mismatch[j]++;
		    if (pst_param->pst_present_mm->ast_mmdata[j].d_enthalpy != 99999){
		      pst_results->d_total_enthalpy += pst_param->pst_present_mm->ast_mmdata[j].d_enthalpy;
		      pst_results->d_total_entropy += pst_param->pst_present_mm->ast_mmdata[j].d_entropy;
		      i_mismatch = FALSE; /* NN for mismatch identified, return to normality */
		    }
		    break;
		  }
		}
		}
		
		if (i_inosine == TRUE){
				/*compare with each possible inosine mismatched pair*/
		for (j = 0 ; j < NBIN ; j++){
		  if ( (strncmp(&pst_param->ps_sequence[i],pst_param->pst_present_inosine->ast_inosinedata[j].s_crick_pair,2) == 0)
		       && (strncmp(&pst_param->ps_complement[i],&pst_param->pst_present_inosine->ast_inosinedata[j].s_crick_pair[3],2) == 0)){
		    pst_results->i_inosine[j]++;
		    if (pst_param->pst_present_inosine->ast_inosinedata[j].d_enthalpy != 99999){
		      pst_results->d_total_enthalpy += pst_param->pst_present_inosine->ast_inosinedata[j].d_enthalpy;
		      pst_results->d_total_entropy += pst_param->pst_present_inosine->ast_inosinedata[j].d_entropy;
		      i_inosine = FALSE; /* NN for inosine mismatch identified, return to normality */
		    }
		  break;
		  }
		}
		}
		
		if (i_mismatch == TRUE || i_inosine == TRUE){
		  fprintf(ERROR," NN parameters for %c%c/%c%c not found. Check the file containing\n"
			  " the information on mismatches\n",pst_param->ps_sequence[i],pst_param->ps_sequence[i+1],pst_param->ps_complement[i],pst_param->ps_complement[i+1]);
		  exit(EXIT_FAILURE);
		}
	    } else 
		/*compare with each possible regular pair*/
		for ( j = 0; j < NBNN; j++)
		    if ( strncmp(&pst_param->ps_sequence[i],pst_param->pst_present_nn->ast_nndata[j].s_crick_pair,2) == 0){
			pst_results->i_crick[j]++;
			pst_results->d_total_enthalpy += pst_param->pst_present_nn->ast_nndata[j].d_enthalpy;
			pst_results->d_total_entropy += pst_param->pst_present_nn->ast_nndata[j].d_entropy;
		    }
	}
	pst_results->d_tm = tm_exact(pst_param,pst_results);
    }
    return pst_results;
}

/********************************************************************
 * The length is too important. approximative computation performed *
 ********************************************************************/

double tm_approx(struct param *pst_param){
    double d_temp;		/* melting temperature */
    int i_size;			/* size of the duplex */
    int i_numbergc;		/* ... */
    double d_percentgc;		/* need an explanation? */
    char *pc_screen;     	/* screen the sequence */

    /*+--------------------+
      | Size of the duplex |
      +--------------------+*/

    i_size = strlen(pst_param->ps_sequence);
    if (i_size == 0){
	fprintf(ERROR," The size of the duplex appears to be null. Therefore I\n"
                      " cannot compute approximation of the melting temperature.\n");
	exit(EXIT_FAILURE);
    }

    /*+----------------+
      | percent of G+C |
      +----------------+*/
	
    pc_screen = pst_param->ps_sequence;
    i_numbergc = 0;
    while(*pc_screen != '\0'){
	if (*pc_screen == 'G' || *pc_screen == 'C')
	    i_numbergc++;
	pc_screen++;
    }
    d_percentgc = ( (double)i_numbergc / (double)i_size ) * 100;
    

    /*+---------------------+
      | melting temperature |
      +---------------------+*/

    if(i_dnadna == TRUE)
	d_temp = 81.5 
	    + 16.6 * log10(pst_param->d_conc_salt / (1.0 + 0.7 * pst_param->d_conc_salt)) 
	    + 0.41 * d_percentgc
	    - 500.0 / (double)i_size;
    else 
	if(i_dnarna == TRUE)
	    d_temp = 67 
		+ 16.6 * log10(pst_param->d_conc_salt / (1.0 + 0.7 * pst_param->d_conc_salt)) 
		+ 0.8 * d_percentgc
		- 500.0 / (double)i_size;
	else 
	    if(i_rnarna == TRUE)
		d_temp = 78 
		    + 16.6 * log10(pst_param->d_conc_salt / (1.0 + 0.7 * pst_param->d_conc_salt)) 
		    + 0.8 * d_percentgc
		    - 500.0 / (double)i_size;
	    else{
		fprintf(ERROR," I do not find any hybridisation type and therefore\n"
                              " I cannot compute the approximative melting temperature\n");
		exit(EXIT_FAILURE);
	    }
    
    return d_temp;
}

/******************************************
 * Nearest-neighbor computation performed *
 ******************************************/

double tm_exact(struct param *pst_param, struct thermodynamic *pst_results){

    int i;			/* loop counters */
    double d_temp;		/* melting temperature */
    double d_temp_na;           /*melting temperature in 1M na+ */
    double d_salt_corr_value = 0.0; /* ... */
    double d_magn_corr_value = 0.0;
    int i_size;			/* size of the duplex */
    int i_numbergc;		/* ... */
    double d_fgc;		/* need an explanation? */
    double d_conc_monovalents = pst_param->d_conc_salt + pst_param->d_conc_potassium + pst_param->d_conc_tris/2;
    double d_ratio_ions = sqrt(pst_param->d_conc_magnesium)/d_conc_monovalents;
    double d_a = 3.92/100000.0; /*Parameters from the article of Owczarzy*/
    double d_b = 9.11/1000000;  /*Parameters from the article of Owczarzy*/
    double d_c = 6.26/100000;   /*Parameters from the article of Owczarzy*/
    double d_d = 1.42/100000;  /*Parameters from the article of Owczarzy*/
    double d_e = 4.82/10000;   /*Parameters from the article of Owczarzy*/
    double d_f = 5.25/10000;   /*Parameters from the article of Owczarzy*/
    double d_g = 8.31/100000;   /*Parameters from the article of Owczarzy*/  
    
    /*+--------------------+
      | Size of the duplex |
      +--------------------+*/

    i_size = strlen(pst_param->ps_sequence);
    if (i_size == 0){
	fprintf(ERROR," The size of the duplex appears to be null. Therefore I\n"
                      " cannot compute approximation of the melting temperature.\n");
	exit(EXIT_FAILURE);
    }
    
      /*+----------------+
      | fraction of G+C |
      +----------------+*/

    i_numbergc = 0;
    for (i = 0; i < strlen(pst_param->ps_sequence); i++){
    switch (pst_param->ps_sequence[i]) {
	case 'G':
	    if (pst_param->ps_complement[i] == 'C')
		i_numbergc++;
		break;
        case 'C':
	    if (pst_param->ps_complement[i] == 'G')
		i_numbergc++;
		break;
	    }
}
    d_fgc = ( (double)i_numbergc / (double)i_size );
    


    /*+-----------------+
      | ion correction |
      +-----------------+*/
      
    if (i_magnesium == FALSE){ /*if [Na+] != 0 and the other ions = 0, we can use the approximations with sodium correction*/
    
    	if (strncmp(pst_param->s_sodium_correction,"wet91a",6) == 0)
		d_salt_corr_value = 16.6 * log10 (pst_param->d_conc_salt / (1.0 + 0.7 * pst_param->d_conc_salt)) - 269.32;
    	else 
	    if (strncmp(pst_param->s_sodium_correction,"san96a",6) == 0)
	    	d_salt_corr_value = 12.5 * log10 (pst_param->d_conc_salt) - 273.15;
	    else 
	        if (strncmp(pst_param->s_sodium_correction,"san98a",6) == 0){
			d_salt_corr_value = -273.15;
			pst_results->d_total_entropy += 0.368 * (strlen(pst_param->ps_sequence)-1) * log (pst_param->d_conc_salt);
	    	} else 
			if (strncmp(pst_param->s_sodium_correction,"nak99a",6) == 0){
		   		 fprintf(ERROR," Sorry, not implemented yet\n");
		   		 exit(EXIT_FAILURE);
			}
    				/* thermodynamic term */
    d_temp = pst_results->d_total_enthalpy / (pst_results->d_total_entropy + 1.987 * log (pst_param->d_conc_probe/pst_param->d_gnat))
				/* salt correction */
	+ d_salt_corr_value;
    }
    else if (i_magnesium == TRUE && i_dnadna == TRUE){ /*The following algorithm is from the article of Owczarzy*/
    	if (d_conc_monovalents == 0) {
		d_magn_corr_value = d_a - d_b * log (pst_param->d_conc_magnesium) + d_fgc * (d_c + d_d * log
		(pst_param->d_conc_magnesium)) + 1/(2 * ((double)i_size - 1)) *
		(- d_e + d_f * log (pst_param->d_conc_magnesium) + d_g *
		pow(log (pst_param->d_conc_magnesium),2));
	}
	else {		
		if (d_ratio_ions < 0.22) {		
			d_magn_corr_value = (4.29 * d_fgc - 3.95) * 1/100000 * log (d_conc_monovalents) + 9.40 * 1/1000000 * pow(log (d_conc_monovalents),2);
		}
		else {
			if (d_ratio_ions < 6) {
			
				d_a = 3.92/100000 * (0.843 - 0.352 * sqrt(d_conc_monovalents) * log (d_conc_monovalents));
				d_d = 1.42/100000 * (1.279 - 4.03/1000 * log (d_conc_monovalents) - 8.03/1000 * pow(log (d_conc_monovalents),2));
				d_g = 8.31/100000 * (0.486 - 0.258 * log (d_conc_monovalents) + 5.25/1000 * pow(log (d_conc_monovalents),3));
				
				d_magn_corr_value = d_a - d_b * log
				(pst_param->d_conc_magnesium) + d_fgc *
				(d_c + d_d * log (pst_param->d_conc_magnesium)) + 1/(2 * ((double)i_size - 1)) * (- d_e + d_f * log (pst_param->d_conc_magnesium) + d_g *
				pow(log (pst_param->d_conc_magnesium),2));
			}
			else {			
				d_magn_corr_value = d_a - d_b * log (pst_param->d_conc_magnesium) + d_fgc * (d_c + d_d * log (pst_param->d_conc_magnesium)) + 1/(2 *
				((double)i_size - 1)) * (- d_e + d_f * log (pst_param->d_conc_magnesium) + d_g *
				pow(log(pst_param->d_conc_magnesium),2));
			
			}
		}
	}
    d_temp_na = pst_results->d_total_enthalpy / (pst_results->d_total_entropy +
    1.987 * log (pst_param->d_conc_probe/pst_param->d_gnat));	
    				/* thermodynamic term with magnesium correction */
    d_temp = 1/(1/d_temp_na + d_magn_corr_value) - 273.15;
    
    }
    else {
    	fprintf(OUTPUT,"  WARNING: The magnesium correction can efficiently\n"
			  "  account only for the DNA/DNA hybridisation. So we can't take in account the magnesium, potassium ant tris concentration for the melting temperature\n" 
			  "computation of RNA or hybrids RNA/DNA duplexes.\n");
    d_temp = 0.0;
    }
    return d_temp;
}
