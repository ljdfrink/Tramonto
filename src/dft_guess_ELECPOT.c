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
 *  FILE: dft_ELECPOT.c
 *
 *  This file sets up an initial guess for an electrostatic potential
 *  field.
 *
 */

#include "dft_guess_ELECPOT.h"
 
/*********************************************************/
/*setup_elec_pot: set up the electrostatic potential initial guess*/
void setup_elec_pot(double **xInBox,int iguess)
{
  int loc_inode,iunk,inode_box,inode,ijk[3];
  double x_dist;

  for (loc_inode=0; loc_inode<Nnodes_per_proc; loc_inode++){
       inode_box = L2B_node[loc_inode];
       iunk = Phys2Unk_first[POISSON];
       if (Type_interface !=UNIFORM_INTERFACE && iguess==LINEAR){
           inode     = B2G_node[inode_box];
           node_to_ijk(inode,ijk); 
           x_dist = Esize_x[Grad_dim]*ijk[Grad_dim];

           xInBox[iunk][inode_box]  = Elec_pot_LBB + 
                                      (Elec_pot_RTF-Elec_pot_LBB)*
                                      x_dist/Size_x[Grad_dim];
       }
       else xInBox[iunk][inode_box] = 0.0;
  }
  return;
}
/************************************************************/
