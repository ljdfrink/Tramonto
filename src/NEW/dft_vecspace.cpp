/*@HEADER
// ***********************************************************************
// 
//                Tramonto: Molecular Theories Modeling Code
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
// Questions? Contact Laura J.D. Frink (ljfrink@sandia.gov)
// 
// ***********************************************************************
//@HEADER
*/
#include <cstdlib>
#include "dft_vecspace.h"


void create_vecspace(int vecspace_type, int N_loc, int * update, DFT_PMACHINE * machine,
		     DFT_VECSPACE ** vecspace){
 
  int * tmp_update = 0;
  if (N_loc>0) {
    tmp_update = new int[N_loc];
    for (int i=0; i<N_loc; i++) tmp_update[i] = update[i];
  }
  *vecspace = (DFT_VECSPACE *) new DFT_VECSPACE;
  (*vecspace)->vecspace_type = vecspace_type;
  (*vecspace)->N_loc = N_loc;
  (*vecspace)->update = tmp_update;
  (*vecspace)->machine = machine;
  // if (vecspace_type==DFT_AZTEC_VECSPACE) {
  //}
  //else {
  //  exit(1);
  //}
  return;
}

void destroy_vecspace(DFT_VECSPACE ** vecspace){

  delete [] (*vecspace)->update;
  delete (DFT_VECSPACE *) (*vecspace);
  return;
}

void dft_random_vector(DFT_VECSPACE * machine, double *x){

  return;
}