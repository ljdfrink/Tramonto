/* This file was automatically generated.  Do not edit! */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAS_VALUES_H)
#include <values.h>
#include <unistd.h>
#include <string.h>
#endif
#include "mpi.h"
#include "az_aztec.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define NCOMP_MAX 5
extern double Betamu_att[NCOMP_MAX];
extern double Rho_b[NCOMP_MAX];
extern double Rho_b_RTF[NCOMP_MAX];
#define NMER_MAX     100
extern double Rho_seg_b[NMER_MAX];
extern double Rho_seg_RTF[NMER_MAX];
#define UNIFORM_INTERFACE  0
extern int Lsteady_state;
extern int Unk2Comp[NMER_MAX];
extern int Lseg_densities;
#define DENSITY        0
#define NEQ_TYPE       13 
extern int Phys2Unk_first[NEQ_TYPE];
double integrand_att_freen_bulk(int iunk,int inode_box,double **x);
#define THETA_PAIRPOT_RCUT    2
double int_stencil(double **x,int inode_box,int iunk,int sten_type);
double integrand_att_freen(int iunk,int inode_box,double **x);
