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

#include <cmath>
#include "kbcLatticeDescriptors.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
struct kbcLbHelpers {
  static_assert(
    std::is_same<typename DESCRIPTOR::category_tag, descriptors::tag::KBC>::value,
    "DESCRIPTOR is tagged as KBC");

  static T equilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], const T uSqr);

  static void computeFneq(Cell<T,DESCRIPTOR> const& cell, T fNeq[DESCRIPTOR::q], T rho, const T u[DESCRIPTOR::d]);

  static T kbcCollision(Cell<T, DESCRIPTOR>& cell, T rho, T u[DESCRIPTOR::d], const T& beta);
  //Only D3Q27 specialisation implemented currently - see kbcLbHelpersD3Q27.h    

};

} //namespace olb

#include "kbcLbHelpersD3Q27.h"

#endif