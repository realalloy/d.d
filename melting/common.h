/******************************************************************************
 *                               MELTING v4.3                                 *
 * This program   computes for a nucleotide probe, the enthalpy, the entropy  *
 * and the melting temperature of the binding to its complementary template.  *
 * Three types of hybridisation are possible: DNA/DNA, DNA/RNA, and RNA/RNA.  *
 *          Copyright (C) Nicolas Le Novère and Marine Dumousseau  1997-2009 *
 *                                                                            *
 * File: common.h                                                             *
 * Date: 01/APR/2009                                                          *
 * Aim : This file contains the definitions of MACRO and variables common to  *
 *       several modules of melting.                                          *
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

/*    Hungarian abbreviations: 

ast_  array of structures
d_    double precision float
i_    integer
pc_   pointer on character
ps_   pointer to a character string
s_    character string

*/

#ifndef COMMON_H
#define COMMON_H

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>MACRO DEFINITIONS<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#define VERSION 4.3         /* current version */
/* The base directory of nn files is set up during the compilation with the 
   instruction: -DNN_BASE=$(NN_DIR) with NN_DIR is specified by the user */
#ifndef NN_BASE		    /* default definition, in case of */
#define NN_BASE "/usr/local/share/melting/NNFILES"
#endif /*NN_BASE*/
/* environment variable targeting the directory containing the nn sets. Used at
   runtime to supersede the path defined during the compilation */
#define NN_PATH NN_PATH	

/* CHECK THAT CAREFULLY, it is clearly needed, particularly 
for a MacOS port. maybe directly in the code rather than as MACRO.  */
/* #define DEFAULT_NN . */	

#define TRUE 1		    /* definition of the truth */
#define FALSE 0		    /* definition of the lie */

/* FILENAME_MAX already defined in stdio.h but dangerous because on certain 
   system it is arbitrarily large, hence the new definition. */
#if FILENAME_MAX > 256          
#define FILE_MAX 256	        
#else                           
#define FILE_MAX FILENAME_MAX
#endif /* FILENAME_MAX */

#define ERROR    stderr	    /* where to send error messages */
#define VERBOSE  stdout	    /* where to send verbose details */
#define MENU     stdout	    /* where to write interactive menus */
#define INPUT     stdin	    /* where to acquire the infos */
#define OUTPUT   stdout	    /* where to send results */
#define MAX_REF     256	    /* max size of a reference to a paper */
#define NUM_REF      16	    /* max number of references for a parameter file. Unelegant TO BE FIXED */
#define MAX_LINE    128	    /* maximum size of an input line */
#define DEFAULT_DNADNA_NN "all97a.nn"           /* default nearest-neighbor set used for DNA/DNA hybridisation */
#define DEFAULT_DNARNA_NN "sug95a.nn"           /* default nearest-neighbor set used for DNA/RNA hybridisation */
#define DEFAULT_RNARNA_NN "xia98a.nn"           /* default nearest-neighbor set used for DNA/RNA hybridisation */
#define DEFAULT_DNADNA_MISMATCHES "dnadnamm.nn" /* default nearest-neighbor set used for DNA/DNA mismatches */ 
#define DEFAULT_DNARNA_MISMATCHES "dnadnamm.nn" /* default nearest-neighbor set used for DNA/RNA mismatches */ 
#define DEFAULT_RNARNA_MISMATCHES "dnadnamm.nn" /* default nearest-neighbor set used for RNA/RNA mismatches */ 
#define DEFAULT_DNADNA_INOSINE_MISMATCHES "san05a.nn" /* default nearest-neighbor set used for DNA/DNA inosine mismatches */ 
#define DEFAULT_DNARNA_INOSINE_MISMATCHES "san05a.nn" /* default nearest-neighbor set used for DNA/RNA inosine mismatches */ 
#define DEFAULT_RNARNA_INOSINE_MISMATCHES "bre07a.nn" /* default nearest-neighbor set used for RNA/RNA inosine mismatches */ 
#define DEFAULT_DNADNA_DANGENDS "dnadnade.nn"   /* default nearest-neighbor set used for DNA/DNA dangling ends */ 
#define DEFAULT_DNARNA_DANGENDS "dnadnade.nn"   /* default nearest-neighbor set used for DNA/RNA dangling ends */ 
#define DEFAULT_RNARNA_DANGENDS "dnadnade.nn"   /* default nearest-neighbor set used for RNA/RNA dangling ends */ 
#define DEFAULT_NUC_CORR 4                      /* default correction for nucleic acid concentration */
#define DEFAULT_SALT_CORR "san98a"              /* default method for the correction of salt concentration */
#define MIN_SALT    0.0	    /* minimal sodium concentration */
#define MAX_SALT   10.0	    /* maximal salt concentration */
#define MIN_PROBE   0.0	    /* minimal nucleic acid concentration */
#define MAX_PROBE   0.1	    /* maximal nucleic acid concentration */
#define MAX_SIZE_NN  60	    /* Maximal size of a duplex being analised by the nearest-neighbor approach */
#define NBNN         18     /* number of regular nn parameters per set */
#define NBMM        240     /* number of mismatch parameters per set */
#define NBIN        109     /* number of inosine mismatch parameters per set */
#define NBDE         64     /* number of dangling end parameters per set */

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>VARIABLE DEFINITIONS<<<<<<<<<<<<<<<<<<<<<<<<<*/

/* Will contain a Crick's pair and the associated calorimetric parameters */
struct calor_const{         
    char     s_crick_pair[6]; /* Six because mismaches are designed as XX/XX plus end of string */
    double   d_enthalpy;
    double   d_entropy;
};

/* contains the parameters for the regular hybridisations */
struct nnset{
    char   s_reference[NUM_REF][MAX_REF]; /* contains the references to the articles */
/*  char   *s_reference[]; */   /* FIXME: A construction like this one would be better */
    struct calor_const ast_nndata[NBNN];  /* parameters for present hybridization*/
    char s_nnfile[FILE_MAX];              /* name of the file containing the nn params */
};

/* contains the parameters for the mismatches*/
struct mmset{
    char   s_reference[NUM_REF][MAX_REF]; /* contains the references to the articles */
/*  char   *s_reference[]; */   /* FIXME: A construction like this one would be better */
    struct calor_const ast_mmdata[NBMM];  /* parameters for present hybridization: 4x4x4x4 possibilities - 16 regular pairs*/
    char s_mmfile[FILE_MAX];              /* name of the file containing the mm params */
};

/* contains the parameters for the inosine mismatches*/
struct inosineset{
    char   s_reference[NUM_REF][MAX_REF]; /* contains the references to the articles */
/*  char   *s_reference[]; */   /* FIXME: A construction like this one would be better */
    struct calor_const ast_inosinedata[NBIN];  /* parameters for present hybridization*/
    char s_inosinefile[FILE_MAX];              /* name of the file containing the insoine params */
};

/* contains the parameters for the dangling ends*/
struct deset{
    char   s_reference[NUM_REF][MAX_REF]; /* contains the references to the articles */
/*  char   *s_reference[]; */   /* FIXME: A construction like this one would be better */
    struct calor_const ast_dedata[NBDE];  /* parameters for present hybridization */
    char s_defile[FILE_MAX];              /* name of the file containing the de params */
};

/* Contains the parameters of the present computation */
struct param {			  
    char *ps_sequence;	          /* sequence of the nucleic acid to test */
    char *ps_complement;          /* sequence of the (supposed) reverse complement */
    double d_conc_probe;	  /* concentration of the strand in excess */
    double d_conc_salt;	          /* concentration in sodium */
    double d_conc_potassium;	  /* concentration in potassium */
    double d_conc_tris;	          /* concentration in tris */
    double d_conc_magnesium;	  /* concentration in magnesium */
    double d_gnat;	          /* correction facteur for the probe concentration */
    struct nnset *pst_present_nn; /* Contains the current nearest-neighbor parameters set */
    struct mmset *pst_present_mm; /* Contains the current parameters for mismatches */
    struct inosineset *pst_present_inosine; /* Contains the current parameters for inosine mismatches */
    struct deset *pst_present_de; /* Contains the current parameters for dangling ends */
    char s_sodium_correction[7];  /* code of the selected salt correction */
    char s_outfile[FILE_MAX];     /* name of the file where to write the results */
};

/* Contains the result of the present analysis*/
struct thermodynamic {       
    double   d_total_enthalpy;  /* enthalpy of the helix-coil transition */
    double   d_total_entropy;   /* entropy of the helix-coil transition */
    double   d_tm;	        /* temperature of halt-denaturation */
    int      i_crick[NBNN];     /* number of each Crick's pair */
    int      i_mismatch[NBMM];  /* number of each mismach */
    int      i_inosine[NBIN];  /* number of each inosine mismach */
    int      i_dangends[NBDE]; /* number of each dangling end*/
};

#endif /* COMMON_H */
