/******************************************************************************
 *                               MELTING v4.3                                 *
 * This program   computes for a nucleotide probe, the enthalpy, the entropy  *
 * and the melting temperature of the binding to its complementary template.  *
 * Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  *
 *          Copyright (C) Nicolas Le Novère and Marine Dumousseau  1997-2013  *
 *                                                                            *
 * File: decode.h                                                             *
 * Date: 01/APR/2009                                                          *
 * Aim : Variable definitions for decode.c                                    *
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

#ifndef DECODE_H
#define DECODE_H

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>MACRO DEFINITIONS<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>VARIABLE DEFINITIONS<<<<<<<<<<<<<<<<<<<<<<<<<*/

int i_alt_nn = FALSE;		 /* an alternative set of nn parameters is required */
int i_alt_mm = FALSE;	         /* an alternative set of mismatches parameters is required */
int i_alt_inosine = FALSE;	 /* an alternative set of inosine mismatches parameters is required */
int i_alt_de = FALSE;	         /* an alternative set of dangling ends parameters is required */
int i_approx = FALSE;		 /* approximative tm computation? */
int i_complement = FALSE;	 /* correct complementary sequence? */
int i_dnadna = FALSE;		 /* those flags specify the type of hybridisation */
int i_dnarna = FALSE;		 /* (useful fo the approximative computations) */
int i_rnarna = FALSE;
int i_hybridtype = FALSE;	 /* correct hybridisation type? */
int i_infile = FALSE;		 /* infile furnished? */
int i_mismatchesneed = FALSE;	 /* We need mismaches parameters */
int i_inosineneed = FALSE;	 /* We need inosine mismaches parameters */
int i_dangendsneed = FALSE;	 /* We need dangling ends parameters */
int i_outfile = FALSE;		 /* outfile requested? */
int i_probe = FALSE;		 /* correct nucleic acid concentration? */
int i_quiet = FALSE;		 /* stay quiet, i.e. no interactive correction of parameters */
int i_salt = FALSE;		 /* correct sodium concentration? */
int i_magnesium = FALSE;          /* can we use magnesium correction algorithm? */
int i_seq = FALSE;		 /* correct sequence? */
int i_verbose = FALSE;		 /* is verbose mode on? */
int i_threshold = MAX_SIZE_NN;   /* threshold before approximative calculus */

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>FUNCTION PROTOTYPES<<<<<<<<<<<<<<<<<<<<<<<<<*/

extern void usage(void);         /*precises the command line parameters*/

struct param *decode_input(struct param *pst_in_param, char *ps_input, char *ps_path);
struct nnset *read_nn(char *ps_nn_set, char *ps_path);         /* read a file containing a nn set */
struct mmset *read_mismatches(char *ps_mm_set, char *ps_path); /* read a file containing a mismatch set */
struct inosineset *read_inosine(char *ps_inosine_set, char *ps_path); /* read a file containing a inosine mismatch set */
struct deset *read_dangends(char *ps_de_set, char *ps_path);   /* read a file containing a dangling ends set */
char *read_string(FILE *stream); /* read a line of input of unknown size */
void legal(void);		 /* precises the copyright under which melting is released*/

#endif /* DECODE_H */



