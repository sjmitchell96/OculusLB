/*  This file is part 
 *
 *  
 *  
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */

#ifndef KBC_LB_HELPERS_H
#define KBC_LB_HELPERS_H

#include "kbcLatticeDescriptors.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
struct kbcLbHelpers {
  static_assert(
    std::is_same<typename DESCRIPTOR::category_tag, descriptors::tag::KBC>::value,
    "DESCRIPTOR is tagged as KBC");

  static void computeFeq ( T fEq[DESCRIPTOR::q], T rho, const T u[DESCRIPTOR::d]); 

  static void computeFneq ( Cell<T,DESCRIPTOR> const& cell,
                            T fNeq[DESCRIPTOR::q], T rho, const T u[DESCRIPTOR::d]);

  static T kbcCollision(Cell<T, DESCRIPTOR>& cell, T rho, T u[DESCRIPTOR::d], const T& beta);
  //Only D3Q27 specialisation implemented currently - see kbcLbHelpersD3Q27.h    

  


  //general templates ... not optimised ... for loops etc. - delibaretely leave empty for now?? (use specialised only...)
  // Or maybe just one big specialised collision method here for d3q27 - don't even bother with seperate Eq, gamma functions ??? - decide tmrw
  // Probably go with above - avoid generalisation as don't need! - leave general collision function empty and go straight to inline specialisation??

  //computeEquilibrium (rho, u, uSqr)
  //static void computeEquilibrium(T fEq[DESCRIPTOR::q], T rho, T u, T uSqr) 
  //{
  

  //}  
  //computeGamma(rho, u)



  //Collision (cell, rho, u) - One big specialised function for d3q27 (in kbcLbHelpers3D.h), and leave this default empty for now?
    //compute feq (rho, u)
    //compute gamma (rho, u)
    //
  

};

} //namespace olb

#include "kbcLbHelpersD3Q27.h"

#endif