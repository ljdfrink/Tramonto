/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
 
/*
 *  FILE: dft_newton.c
 *
 *  This file solves the dft problem using Newton's method.
 */
#include "mpi.h" 


#ifdef HAVE_DFT_SCHUR_SOLVER
#include "dft_schur_solver.h"
#endif

#include "dft_globals_const.h" 
#include "rf_allo.h"
/* #include "dft_basic_lin_prob_mgr_wrapper.h"
   #include "dft_poly_lin_prob_mgr_wrapper.h" */

static void print_resid_norm(int iter);
void fill_test(double **x, int flag);

#ifdef NUMERICAL_JACOBIAN
static void do_numerical_jacobian(double **x);
#endif


/* NEWTON SOLVER using dft_Linprobmgr */
int solve_problem(double **x, double **x2)
/*
 * x and x2 in dft_main are allocated [Nunk_per_node][Nnodes_box]
 * x2 only relevent for Lbinodal calculations, so can mostly be ignored
 */

{
  int iter,iunk,i;
  double **xOwned;
  int geq[NMER_MAX], ginveq[NMER_MAX], cmseq[NCOMP_MAX], densityeq[NCOMP_MAX] ;
  int indnonlocaleq[NMER_MAX], depnonlocaleq[NMER_MAX];
/*  int gequ[] = {6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42};
  int ginvequ[] = {43, 41, 39, 37, 35, 33, 31, 29, 27, 25, 23, 21, 19, 17, 15, 13, 11, 9, 7};
  int cmsequ[] = {0, 1, 2};
  int densityequ[] = { 3, 4, 5};*/
  int count_density,count_cms_field,count_geqn,count_ginv_eqn,count_indnonlocal,count_depnonlocal,index_save;
  int one_particle_size;


  /* Construct dft_Linprobmgr with information on number of unknowns*/
  int is_poly = 0;
  int debug = 0;
 if (Type_poly != NONE) {
   
   count_density=count_cms_field=count_geqn=count_ginv_eqn=0;
   for (iunk=0;iunk<Nunk_per_node;iunk++){
     switch(Unk2Phys[iunk]){
     case DENSITY:
       densityeq[count_density++]=iunk; break; 
     case POISSON:
     case CMS_FIELD:                  
       cmseq[count_cms_field++]=iunk; break; 
     case CMS_G:                  
       if ((iunk-Geqn_start[0])%2 == 0)
	 geq[count_geqn++]=iunk; 
       else
	 ginveq[count_ginv_eqn++]=iunk; 
       break;
     } 
   }
   /* now invert the order of the ginverse equations ! */
   for (i=0;i<count_ginv_eqn/2;i++){
     index_save = ginveq[i];
     ginveq[i] = ginveq[count_ginv_eqn-1-i];
     ginveq[count_ginv_eqn-1-i]=index_save;
   }
   LinProbMgr_manager = dft_poly_lin_prob_mgr_create(Nunk_per_node, Aztec.options, Aztec.params, MPI_COMM_WORLD);
   dft_poly_lin_prob_mgr_setgequationids(LinProbMgr_manager, Ngeqn_tot/2, geq);
   dft_poly_lin_prob_mgr_setginvequationids(LinProbMgr_manager, Ngeqn_tot/2, ginveq);
   dft_poly_lin_prob_mgr_setcmsequationids(LinProbMgr_manager, Ncomp, cmseq);
   dft_poly_lin_prob_mgr_setdensityequationids(LinProbMgr_manager, Ncomp, densityeq);
   /*dft_poly_lin_prob_mgr_setfieldondensityislinear(LinProbMgr_manager,TRUE);*/
 }
 else if (Type_func > NONE) {

   count_density=count_indnonlocal=count_depnonlocal=0;
   one_particle_size=FALSE;
   if ((Lhard_surf && Nlists_HW == 2) || (!Lhard_surf && Nlists_HW == 1)) one_particle_size=TRUE;
   for (iunk=0;iunk<Nunk_per_node;iunk++){
     switch(Unk2Phys[iunk]){
     case POISSON:
     case DENSITY:
       densityeq[count_density++]=iunk; break; 
     case RHOBAR_ROSEN:                  
       if (one_particle_size && (
             (iunk > Phys2Unk_first[RHOBAR_ROSEN]+1 && 
              iunk < Phys2Unk_first[RHOBAR_ROSEN]+Nrho_bar_s)  ||
              iunk >= Phys2Unk_first[RHOBAR_ROSEN]+Nrho_bar_s+Ndim) )
	 depnonlocaleq[count_depnonlocal++]=iunk; 
       else
	 indnonlocaleq[count_indnonlocal++]=iunk; 
       break;
     case BOND_WTC:
     case CAVITY_WTC:
       indnonlocaleq[count_indnonlocal++]=iunk; break;
     } 
   }
   LinProbMgr_manager = dft_hardsphere_lin_prob_mgr_create(Nunk_per_node, Aztec.options, Aztec.params, MPI_COMM_WORLD);
   dft_hardsphere_lin_prob_mgr_setindnonlocalequationids(LinProbMgr_manager, count_indnonlocal, indnonlocaleq);
   dft_hardsphere_lin_prob_mgr_setdepnonlocalequationids(LinProbMgr_manager, count_depnonlocal, depnonlocaleq);
   dft_hardsphere_lin_prob_mgr_setdensityequationids(LinProbMgr_manager, count_density, densityeq);
   if (Type_attr != NONE || Type_poly_TC || Mesh_coarsening || Type_coul != NONE)
                dft_hardsphere_lin_prob_mgr_seta22blockisdiagonal(LinProbMgr_manager, FALSE);
   else         dft_hardsphere_lin_prob_mgr_seta22blockisdiagonal(LinProbMgr_manager, TRUE);
 }
 else{
   LinProbMgr_manager = dft_basic_lin_prob_mgr_create(Nunk_per_node, Aztec.options, Aztec.params, MPI_COMM_WORLD);
 }


  /* Give Nodal Row and Column maps */
  (void) dft_linprobmgr_setnodalrowmap(LinProbMgr_manager, Nnodes_per_proc, L2G_node);
  (void) dft_linprobmgr_setnodalcolmap(LinProbMgr_manager, Nnodes_box     , B2G_node);

  /* Linprobmgr can now set up its own numbering scheme, set up unknown-based Maps */
  (void) dft_linprobmgr_finalizeblockstructure(LinProbMgr_manager);

  /* Set initial guess on owned nodes and reconcile ghost nodes using importr2c */
  xOwned = (double **) array_alloc(2, Nunk_per_node, Nnodes_per_proc, sizeof(double));
  set_initial_guess(Iguess1, xOwned);

  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xOwned, x);

  /* If requested, write out initial guess */
   if (Iwrite == VERBOSE) print_profile_box(x,"rho_init.dat");

  /* Do same for second solution vector when Lbinodal is true */
  if (Lbinodal) {
    set_initial_guess(BINODAL_FLAG, xOwned);
    (void) dft_linprobmgr_importr2c(LinProbMgr_manager, xOwned, x2);
    if (Iwrite == VERBOSE) print_profile_box(x2,"rho_init2.dat");
  }
  safe_free((void **) &xOwned);

  if (Loca.method != -1)
    iter = solve_continuation(x, x2);
  else
     iter = newton_solver(x, NULL);

  /* Call the destructor for the dft_Linprobmgr */
  dft_linprobmgr_destruct(LinProbMgr_manager);

  return(iter);
}

/****************************************************************************/
/****************************************************************************/
void box2owned(double** xBox, double** xOwned) {
  /* Routine to restrict 2d box array to 2d node array */
  /* This is the opposite of importr2c, and strips ghost nodes */
  int iunk, inode, ibox;

  for (ibox=0; ibox<Nnodes_box; ibox++) {
    inode = B2L_node[ibox];
    if (inode!=-1) {
      for (iunk=0; iunk<Nunk_per_node; iunk++)
        xOwned[iunk][inode] = xBox[iunk][ibox];
    }
  }
}

/****************************************************************************/
/****************************************************************************/
int newton_solver(double** x, void* con_ptr) {

  int iter=0;
  int iunk, ibox;
  int converged=FALSE, converged2=TRUE;
  double** delta_x;
  delta_x = (double **) array_alloc(2, Nunk_per_node, Nnodes_box, sizeof(double));
  char filename[20]="matrix.dat";

  do {
    iter++;

    (void) dft_linprobmgr_initializeproblemvalues(LinProbMgr_manager);

    /* Call Matrix and Residual Fill routine, resid_only_flag=FALSE)*/
    fill_resid_and_matrix_control(x, iter,FALSE); 
    /*fill_test(x, FALSE);*/

    (void) dft_linprobmgr_finalizeproblemvalues(LinProbMgr_manager);
    if (Iwrite != NO_SCREEN) print_resid_norm(iter);
    (void) dft_linprobmgr_setupsolver(LinProbMgr_manager);
    (void) dft_linprobmgr_solve(LinProbMgr_manager);
    
    /* I am assuming getLhs returns box coordinates (e.g. Column Map)!! */
    (void) dft_linprobmgr_getlhs(LinProbMgr_manager, delta_x);
    /*
      for (iunk=0; iunk<Nunk_per_node; iunk++)
      for (ibox=0; ibox<Nnodes_box; ibox++) {
      printf("delta_x[%d][%d] = %g\n", iunk, ibox,delta_x[iunk][ibox]);
      }
    */
#ifdef NUMERICAL_JACOBIAN
    do_numerical_jacobian(x);
#endif
  //dft_linprobmgr_writeMatrix(LinProbMgr_manager,filename,NULL,NULL);
    
    if (con_ptr != NULL) converged2 =
      continuation_hook_conwrap(x, delta_x, con_ptr, Newton_rel_tol, Newton_abs_tol);

    /* Do: x += delta_x, and check for convergence */
    converged = update_solution(x, delta_x, iter);



  } while (iter < Max_Newton_iter && (!converged || !converged2));

  if (!converged || !converged2) {
    if (Proc==0 && Iwrite!=NO_SCREEN) printf("\tNewton Solver: Failed to converge in %d iterations\n",iter);
    iter = -iter;
  }
  else
    if (Proc==0 && Iwrite!=NO_SCREEN) printf("\tNewton Solver: Successful convergence in %d iterations\n",iter);

  safe_free((void **) &delta_x);
  return iter;
}

/****************************************************************************/
/****************************************************************************/
int update_solution(double** x, double** delta_x, int iter) {
/* Routine to update solution vector x using delta_x and 
 * to test for convergence of Newton's method.
 * 
 * Note: Modifications to Newton's method, including damping, have not yet
 *       been translated form previous update_solution.
 *       iter  value may be used in some damping methods
 */

  int iunk, ibox, inode;
  double updateNorm=0.0, temp,frac_min,frac;
  char *yo = "newupdate solution";

   /* Certain unknowns - specifically densities and Gs in CMS DFT cannot be less than 0.
      Here we locate problems, and scale the entire update vector to prevent this from 
      happening. */
  frac_min=1.0;
  for (ibox=0; ibox<Nnodes_box; ibox++) { /* find minimum update fraction in entire domain */
    for (iunk=0; iunk<Nunk_per_node; iunk++){
      if (Unk2Phys[iunk]==CMS_G || Unk2Phys[iunk]==DENSITY || Unk2Phys[iunk]==CMS_FIELD){
         if(x[iunk][ibox]+delta_x[iunk][ibox]<0.0){
             frac = min(1.0,x[iunk][ibox]/(-delta_x[iunk][ibox]));
             frac = max(frac,Min_update_frac);
         } 
         else frac=1.0;

         if (frac<frac_min) frac_min=frac;
      }
    }
  }

  frac_min=gmin_double(frac_min);
  if (Proc==0 && Iwrite != NO_SCREEN)
      printf("\tUPDATE FRAC = %g percent\n",frac_min*100);

  for (ibox=0; ibox<Nnodes_box; ibox++) {

    /* Increment updateNorm only for owned nodes (inode=-1 for ghosts) */
    inode = B2L_node[ibox];
    if (inode != -1) {
      for (iunk=0; iunk<Nunk_per_node; iunk++) {
        temp =(frac_min*delta_x[iunk][ibox])/(Newton_rel_tol*x[iunk][ibox] + Newton_abs_tol);
        updateNorm +=  temp*temp;
      }
    }

    /* Update all solution componenets */
    for (iunk=0; iunk<Nunk_per_node; iunk++){
      if ((Unk2Phys[iunk]==DENSITY || Unk2Phys[iunk]==CMS_G || Unk2Phys[iunk]==CMS_FIELD || Unk2Phys[iunk]==BOND_WTC || Unk2Phys[iunk]==CAVITY_WTC) && 
            x[iunk][ibox]+frac_min*delta_x[iunk][ibox] <1.e-15){
            x[iunk][ibox]=0.1*x[iunk][ibox];
      }
      else if ((iunk==Phys2Unk_first[RHOBAR_ROSEN] || iunk==(Phys2Unk_first[CAVITY_WTC]+1)) && 
              x[iunk][ibox]+frac_min*delta_x[iunk][ibox] > 1.0){
              x[iunk][ibox]+=0.5*(1.0-x[iunk][ibox]);
      }
      else{
         x[iunk][ibox] += frac_min*delta_x[iunk][ibox];
      }
    }
  }

  updateNorm = sqrt(gsum_double(updateNorm));

  if (Proc==0 && Iwrite != NO_SCREEN)
    printf("\n\t\t%s: Weighted norm of update vector =  %g\n", yo, updateNorm);

  if (updateNorm > 1.0) return(FALSE);
  else                  return(TRUE);

}
/*****************************************************************************************************/
#ifdef NUMERICAL_JACOBIAN
static void do_numerical_jacobian(double **x)
/* This routine compares the analytic and numerical jacobians for the     */
/* purpose of checking the accuracy of an analytic Jacobian. It is only   */
/* written for serial runs and will be slow, so only for small problems.  */
/* Three files are created: ja, jn, jd  containing the analytic jacobian, */
/* numerical jacobian, and the difference between them, printed as the    */
/* row, column, and value.                                                */
{
  double **full=NULL, del;
  int i, j, N=Nunk_per_node*Nnodes, count=0,iunk,junk,jnode,inode;
  char filename[20];
  FILE *ifp;
  double **resid, **resid_tmp;
  resid = (double **) array_alloc(2, Nunk_per_node, Nnodes_per_proc, sizeof(double));
  resid_tmp = (double **) array_alloc(2, Nunk_per_node, Nnodes_per_proc, sizeof(double));

  dft_linprobmgr_getrhs(LinProbMgr_manager, resid);

  /* print out analytic matrix in MSR format */
  sprintf(filename, "ja%0d",Proc);
  dft_linprobmgr_writeMatrix(LinProbMgr_manager,filename,NULL,NULL);

  /* compute numerical jacobian by repeated residual fill calls */

  full = (double **) array_alloc(2,N,N,sizeof(double));
  if (full==NULL){printf("Not enough memory for full numerical jacobian\n"); exit(-1);}

  for (iunk=0;iunk<Nunk_per_node;iunk++){
    for (inode=0; inode<Nnodes; inode++) {
      i=inode+Nnodes*iunk;
      del=1.e-6*fabs(x[iunk][inode])+1.e-12;
      x[iunk][inode] += del;

      for (junk=0; junk<Nunk_per_node; junk++) 
          for (jnode=0; jnode<Nnodes_per_proc; jnode++) resid_tmp[junk][jnode] = 0.0;

      (void) dft_linprobmgr_initializeproblemvalues(LinProbMgr_manager);
      if (Sten_Type[POLYMER_CR]) fill_resid_and_matrix_P(x,1,TRUE,NODAL_FLAG);
      else                       fill_resid_and_matrix(x,1,TRUE,NODAL_FLAG);

      dft_linprobmgr_getrhs(LinProbMgr_manager, resid_tmp);

      for (junk=0; junk<Nunk_per_node; junk++){ 
          for (jnode=0; jnode<Nnodes_per_proc; jnode++){ 
          j=jnode+Nnodes*junk;
          full[j][i] = (resid[junk][jnode] - resid_tmp[junk][jnode])/del;
      }}
      x[iunk][inode] -= del;
      if (count==100 || i==N-1){
        if (Iwrite != NO_SCREEN)printf("Proc: %d :: ith row=%d, out of %d rows.\n",Proc,i,N-1);
        count=0;
      } count++;
    }
  }
  if (Iwrite != NO_SCREEN)printf("Proc: %d finished computing numerical jacobian !!\n",Proc);

  /* print out nonzero entries of numerical jacobian */
  sprintf(filename, "jn%0d",Proc);
  ifp = fopen(filename,"w");
  for (i=0; i<Nnodes*Nunk_per_node; i++) {
    for (j=0; j<Nnodes*Nunk_per_node; j++) {
      if (fabs(full[i][j]) > 1.0e-10)
       fprintf(ifp,"%d  %d   %g\n",i,j,full[i][j]);
    }
  }
  fclose(ifp);
}
#endif

/*
  / compute and print the difference between the two jacobans /
  sprintf(filename, "jd%0d",Proc);
  ifp = fopen(filename,"w");
  for (i=0; i< N; i++) {
     full[i][i] -= Aztec.val[i];
     for (j= Aztec.bindx[i]; j<Aztec.bindx[i+1]; j++) {
       full[i][Aztec.bindx[j]] -= Aztec.val[j];
     }
   }

  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      if (fabs(full[i][j]) > 1.e-4)
       fprintf(ifp,"%d  %d   %g\n",i,j,full[i][j]);
    }
  }
  fclose(ifp);
  safe_free((void *) &full);
  if (Num_Proc>1) MPI_Barrier(MPI_COMM_WORLD);
    if (Iwrite != NO_SCREEN)printf("Proc: %d KILLING CODE AT END OF NUMERICAL JACOBIAN PRINT\n",Proc);
  exit(-1);
}
 endif*/
/****************************************************************************/


/*****************************************************************************************************/
static void print_resid_norm(int iter)
{
  int iunk, j;
  double norm=0.0;
  double **f;
  FILE *ifp;
  char filename[20]="Resid2.dat";

  if (Proc==0 && Iwrite==VERBOSE) ifp=fopen(filename,"w+");
  f = (double **) array_alloc(2, Nunk_per_node, Nnodes_per_proc, sizeof(double));
  dft_linprobmgr_getrhs(LinProbMgr_manager, f);

  for (j=0; j< Nnodes_per_proc; j++) {
    for (iunk=0; iunk<Nunk_per_node; iunk++) {
       norm += f[iunk][j] * f[iunk][j];
       if(Proc==0 && Iwrite==VERBOSE) fprintf(ifp," %d  %d  %14.11f  %14.11f\n",iunk,L2G_node[j],f[iunk][j],norm);
    }
  }
  if (Proc==0 && Iwrite==VERBOSE) fclose(ifp);

  safe_free((void **) &f);
  norm = gsum_double(norm);
  if (Proc==0) printf("\t\tResidual norm at iteration %d = %g\n",iter, sqrt(norm));
}
/*****************************************************************************************************/

void fill_test(double **x, int flag) 
/* Quick test problem for matrix loading, solve, and Newton's method */
/* For the first unknown at a node, x[0]=inode 
 * for subsequent unknowns at that node, x[iunk]^2=x[iunk-1]
 * So, when run at 3 unknowns per node, global node 16 should have
 * solutions 16, +-4, +-2, for the three unknowns at that node
 */
{
  int loc_inode, inode_box, iunk, inode;
  double f;

  for (loc_inode=0; loc_inode< Nnodes_per_proc; loc_inode++) {
    inode = L2G_node[loc_inode];
    inode_box = L2B_node[loc_inode];
    iunk = 0;

    f = x[iunk][inode_box] - inode - 1;
    dft_linprobmgr_insertrhsvalue(LinProbMgr_manager, iunk, loc_inode, -f);
    dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box,1.0);

    /* For iunk 1 and higher, the square of the unknow is the previous value */
    for (iunk=1; iunk<Nunk_per_node; iunk++) {
      f = x[iunk][inode_box]*x[iunk][inode_box] - x[iunk-1][inode_box];
      dft_linprobmgr_insertrhsvalue(LinProbMgr_manager, iunk, loc_inode, -f);
      dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk,inode_box, 2*x[iunk][inode_box]);
      dft_linprobmgr_insertonematrixvalue(LinProbMgr_manager,iunk,loc_inode,iunk-1,inode_box, -1.0);
    }
  } 
}

