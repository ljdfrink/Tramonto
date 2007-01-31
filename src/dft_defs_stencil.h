/* This file was automatically generated.  Do not edit! */
int stencil_deltaLogical(int sten);
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAS_VALUES_H)
#include <values.h>
#include <unistd.h>
#endif
#include "mpi.h"
#include "az_aztec.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define THETA_FN_SIG          5
#define WTC          2
#define THETA_CR_RPM_MSA      3
#define DELTAC     1 
extern int Type_coul;
#define THETA_PAIRPOT_RCUT    2
extern int Type_attr;
#define THETA_FN_R            1
#define DELTA_FN_R            0
#define NONE       -1
#define NONE      -1
#define NONE -1
#define NONE        -1
extern int Type_func;
#define THETA_CR_DATA         4
#define TRUE  1
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(TRUE) && !defined(_CON_CONST_H_)
#define TRUE  1
#endif
#define DELTA_FN_BOND         6
#define CMS_SCFT     1
#define CMS          0
extern int Type_poly;
#define FALSE 0
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
#define NSTEN        7
extern int Sten_Type[NSTEN];
void setup_stencil_logicals();