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
#define PI    M_PI
double uCOULOMB_Integral(double r,int i,int j);
double uCOULOMB_ATT_CnoS(double r,int i,int j);
double uCOULOMB_DERIV1D(double r,double x,double z1,double z2);
double uCOULOMB_ATT_noCS(double r,int i,int j);
#define NCOMP_MAX 5
extern double Sigma_ff[NCOMP_MAX][NCOMP_MAX];
double uCOULOMB_ATT_CS(double r,int i,int j);
double uCOULOMB_CS_DERIV1D(double r,double x,double z1,double z2,double rcut);
void uCOULOMB_setparams(int context,int i,int j,double *param1,double *param2,double *param3);
#define NWALL_MAX_TYPE 50 
extern double Cut_ww[NWALL_MAX_TYPE][NWALL_MAX_TYPE];
#define WALL_WALL   2
#define NWALL_MAX 600 
extern int WallType[NWALL_MAX];
extern double Cut_wf[NCOMP_MAX][NWALL_MAX_TYPE];
extern double Elec_param_w[NWALL_MAX];
extern double Charge_f[NCOMP_MAX];
#define WALL_FLUID  1
extern double Cut_ff[NCOMP_MAX][NCOMP_MAX];
extern double *Charge;
#define FLUID_FLUID 0
void uCOULOMB_CS_setparams(int context,int i,int j,double *param1,double *param2,double *param3);
double uCOULOMB(double r,double z1,double z2);
extern double Temp_elec;
double uCOULOMB_CS(double r,double z1,double z2,double rcut);
