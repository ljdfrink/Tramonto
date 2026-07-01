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
#include "az_aztec_defs.h"
#include "rf_allo.h"
#include "dft_basic_lin_prob_mgr_wrapper.h"
#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "Tramonto_ConfigDefs.h"
#define TRUE  1
#define NDIM_MAX  3
#define NWALL_MAX_TYPE 20
#define NCOMP_MAX 6
extern double Esize_x[NDIM_MAX];
int offset_to_node_box(int *ijk_box,int *offset,int *reflect_flag);
void node_box_to_ijk_box(int node_box,int *ijk_box);
extern double Area_surf_el[3];
extern int ***Surf_normal;
extern int **Nelems_S;
#define NONE       -1
#define CONST_POTENTIAL  1
#define CONST_CHARGE     2
#define CHARGE_REGULATED 4
/* These are choices for reaction types at surfaces */
#define CHARGEREG1_SURFDENS 0
#define CHARGEREG1_POT      1
#define CHARGEREG2_SURFDENS 2  /* Chan Model */
#define NREACT_MAX 10
extern int ChargeRegTF;
extern int Nreact_perSurf[NWALL_MAX_TYPE];
extern double Density_rxnSites[NWALL_MAX_TYPE];
extern double **Charge_w_sum_els;
extern int Rxn_ID[NWALL_MAX_TYPE][NREACT_MAX];
extern int Rxn_type[NWALL_MAX_TYPE];
extern double Frac_reactSites[NWALL_MAX_TYPE][NREACT_MAX];
extern double K_react[NWALL_MAX_TYPE][NREACT_MAX];
extern double DeltaCharge_frwdRxn[NWALL_MAX_TYPE][NREACT_MAX];
extern int N_fluidCompReact[NWALL_MAX_TYPE][NREACT_MAX];
extern int N_fluidCompProd[NWALL_MAX_TYPE][NREACT_MAX];
extern int ***Fluid_reactComp;
extern int ***Fluid_prodComp;
extern double  Sigma_ff[NCOMP_MAX][NCOMP_MAX];
extern int     Type_func;

#define NWALL_MAX 600 
extern int WallType[NWALL_MAX];
#define NWALL_MAX_TYPE 20 
extern int Type_bc_elec[NWALL_MAX_TYPE];
extern int Nlists_HW;
extern int *B2L_node;
double integrand_surface_charge(int iunk,int inode_box,int iwall,double **x);
double integrand_surface_chargeOLD(int iunk,int inode_box,int iwall,double **x);
#define PI    3.141592653589793238462643383279502884197169399375
extern double Temp_elec;
extern int Nodes_x[NDIM_MAX];
extern int Ndim;
extern int *B2G_node;
void node_to_ijk(int node,int *ijk);
double integrand_maxwell_stress_freen(int iunk,int inode_box,double **x);
extern double *Deltac_b;
#define NCOMP_MAX 6
extern double Rho_b[NCOMP_MAX];
extern double Rho_b_RTF[NCOMP_MAX];
#define NMER_MAX     300
extern double Rho_seg_b[NMER_MAX];
extern double Rho_seg_RTF[NMER_MAX];
#define UNIFORM_INTERFACE  0
extern int Type_interface;
extern int Icomp_to_polID[NCOMP_MAX];
extern int Grafted[NCOMP_MAX];
#define WJDC3        5 
extern int Type_poly;
#define FALSE 0
#if !defined(_CON_CONST_H_)
#define _CON_CONST_H_
#endif
#if !defined(FALSE) && !defined(_CON_CONST_H_)
#define FALSE 0
#endif
extern int Grafted_Logical;
double integrand_elec_MSAcorr_freen_bulk(int iunk,int inode_box,double **x);
#define THETA_CR_GENERAL_MSA  7
#define DELTAC_GENERAL 2
#define THETA_CR_RPM_MSA      3
double int_stencil(double **x,int inode_box,int iunk,int sten_type);
#define DELTAC_RPM     1 
extern int Type_coul;
double integrand_elec_MSAcorr_freen(int iunk,int inode_box,double **x);
extern double Charge_f[NCOMP_MAX];
#define POISSON        1
#define DENSITY        0
#define NEQ_TYPE       12 
extern int Phys2Unk_first[NEQ_TYPE];
extern int Unk2Comp[NMER_MAX];
extern int Lseg_densities;
double integrand_elec_PB_freen(int iunk,int inode_box,double **x);
double calc_deriv_epot(int,int,int *,double **);
double calc_deriv_epot(int idim,int inode0,int *int_type,double **x);
