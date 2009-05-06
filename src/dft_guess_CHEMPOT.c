/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/*
 *  FILE: dft_guess_CHEMPOT.c
 *
 *  This file contains a routine to set an initial guess for the
 *  chemical potential variable (diffusion problems).
 *
 */

#include "dft_guess_CHEMPOT.h"
 
/************************************************************/
/* setup_chem_pot: for cases with steady state profiles,
   set up regions of constant (electro)chemical potentials */
void setup_chem_pot(double **xOwned)
{
  int loc_inode,inode,ijk[3],icomp,iunk,i,nloop;
  double x_dist,x_tot;

  if (Lseg_densities) nloop=Nseg_tot;
  else                nloop=Ncomp;

  x_tot = Size_x[Grad_dim]-2.*X_const_mu;
  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
     inode     = L2G_node[loc_inode];
     node_to_ijk(inode,ijk); 
     x_dist = Esize_x[Grad_dim]*ijk[Grad_dim]-X_const_mu;

     for (i=0; i<nloop; i++){
        iunk = i+Phys2Unk_first[DIFFUSION];
        if (Lseg_densities) icomp=Unk2Comp[i];
        else                icomp=i;

        if (!Zero_density_TF[L2B_node[loc_inode]][icomp]){
           if (Ipot_ff_c == 1){
             xOwned[iunk][loc_inode] = log(xOwned[Phys2Unk_first[DENSITY]+i][loc_inode])
                           + Charge_f[icomp]*(xOwned[Phys2Unk_first[POISSON]][loc_inode]);

           }
           else{
               if (x_dist<0.)           xOwned[iunk][loc_inode]=Betamu_LBB[i];
               else if (x_dist > x_tot) xOwned[iunk][loc_inode]=Betamu_RTF[i];
               else  xOwned[iunk][loc_inode] = Betamu_LBB[i] + (Betamu_RTF[i]-Betamu_LBB[i])* x_dist/x_tot;
           }
       }
       else xOwned[iunk][loc_inode] = -VEXT_MAX;
     }
  }
  return;
}
/************************************************************/

