/******************************************************************************
 *                               MELTING v4.3                                 *
 * This program   computes for a nucleotide probe, the enthalpy, the entropy  *
 * and the melting temperature of the binding to its complementary template.  *
 * Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  *
 *          Copyright (C) Nicolas Le Novère and Marine Dumousseau  1997-2013  *
 *                                                                            *
 * File: calcul.h                                                             *
 * Date: 01/APR/2009                                                          *
 * Aim : Variable definitions for calcul.c                                    *
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

#ifndef CALCUL_H
#define CALCUL_H

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>MACRO DEFINITIONS<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>VARIABLE DEFINITIONS<<<<<<<<<<<<<<<<<<<<<<<<<*/

extern int i_approx;		/* True if we use an approximative method instead of nearest-neighbors */
extern int i_dnadna;		/* those flags specify the type of hybridisation */
extern int i_dnarna;		/* (useful fo the approximative computations) */
extern int i_rnarna;
extern int i_complement;	/* correct complementary sequence?  */
extern int i_alt_de;
extern int i_alt_mm;
extern int i_alt_inosine;
extern int i_threshold;         /* threshold before approximative calculus */
extern int i_magnesium;		/* can we use the magnesium correction algorithm? */

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>FUNCTION PROTOTYPES<<<<<<<<<<<<<<<<<<<<<<<<<*/
struct thermodynamic *get_results(struct param *pst_param);
double tm_approx(struct param *pst_param);
double tm_exact(struct param *pst_param, struct thermodynamic *pst_results);

#endif /* CALCUL_H */




