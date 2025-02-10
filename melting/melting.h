/******************************************************************************
 *                               MELTING v4.3                                 *
 * This program   computes for a nucleotide probe, the enthalpy, the entropy  *
 * and the melting temperature of the binding to its complementary template.  *
 * Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  *
 *          Copyright (C) Nicolas Le Novère and Marine Dumousseau  1997-2013  *
 *                                                                            *
 * File: melting.h                                                            *
 * Date: 01/APR/2009                                                          *
 * Aim : This file contains the definitions of MACRO and variables as well as *
 *       function prototypes for the module melting.c                         *
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

#ifndef MELTING_H
#define MELTING_H

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>MACRO DEFINITIONS<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>VARIABLE DEFINITIONS<<<<<<<<<<<<<<<<<<<<<<<<<*/

extern int i_alt_nn;		/* an alternative set of nn parameters? */
extern int i_alt_mm;	        /* an alternative set of mismatches parameters is required */
extern int i_alt_inosine;	/* an alternative set of inosine mismatches parameters is required */
extern int i_alt_de;	        /* an alternative set of dangling ends parameters is required */
extern int i_verbose;		/* is verbose mode on? */
extern int i_complement;	/* correct complementary sequence? */
extern int i_infile;		/* infile firnished? */
extern int i_outfile;		/* outfile requested? */
extern int i_hybridtype;	/* correct hybridisation type? */
extern int i_salt;		/* correct sodium concentration? */
extern int i_magnesium;		/* can we use the magnesium correction algorithm? */
extern int i_probe;		/* correct nucleic acid concentration? */
extern int i_seq;		/* correct sequence? */
extern int i_dnadna;		/* those flags specify the type of hybridisation */
extern int i_dnarna;		/* (useful fo the approximative computations) */
extern int i_rnarna;
extern int i_approx;		/* approximative tm computation? */
extern int i_quiet;		/* stay quiet, i.e. no interactive correction of parameters */
extern int i_mismatchesneed;	/* We need mismaches parameters */
extern int i_inosineneed;	/* We need mismaches parameters */
extern int i_dangendsneed;	/* We need dangling end parameters */

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>FUNCTION PROTOTYPES<<<<<<<<<<<<<<<<<<<<<<<<<*/

extern struct nnset *read_nn(char *ps_nn_set, char *ps_path);         /* read a file containing the nn set */
extern struct mmset *read_mismatches(char *ps_mm_set, char *ps_path); /* read a file containing a mismatch set */
extern struct inosineset *read_inosine(char *ps_inosine_set, char *ps_path); /* read a file containing a inosine mismatch set */
extern struct deset *read_dangends(char *ps_de_set, char *ps_path);   /* read a file containing a dangling ends set */

extern struct thermodynamic *get_results(struct param *pst_param);
extern struct param *decode_input(struct param *pst_in_param, char *ps_input, char *ps_path);	
                                 /* decodes input line (command-line or inputfile) */
extern char *read_string(FILE *stream); /* read a line of input of unknown size */
extern void legal(void);		/* precises the license under which melting is released */

void usage(void);		/* precises the command line parameters*/

int check_sequence(char *ps_sequence); /* check the legality of every sequence base */
char *make_complement(char *ps_sequence); /* construct the reverse complement from a sequence */

#endif /* MELTING_H */
