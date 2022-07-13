/*  
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
 * can be instantiated -- generic implementation.
 */
#ifndef KBC_MOMENTA_HH
#define KBC_MOMENTA_HH

#include "kbcMomenta.h"
#include "lbHelpers.h"
#include "kbcLbHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
T KBCBulkMomenta<T,DESCRIPTOR>::computeRho(Cell<T,DESCRIPTOR> const& cell) const
{
  return lbHelpers<T,DESCRIPTOR>::computeRho(cell);
}

template<typename T, typename DESCRIPTOR>
void KBCBulkMomenta<T,DESCRIPTOR>::computeU(Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d]) const
{
  T dummyRho;
  lbHelpers<T,DESCRIPTOR>::computeRhoU(cell, dummyRho, u);
}

template<typename T, typename DESCRIPTOR>
void KBCBulkMomenta<T,DESCRIPTOR>::computeJ(Cell<T,DESCRIPTOR> const& cell, T j[DESCRIPTOR::d]) const
{
  lbHelpers<T,DESCRIPTOR>::computeJ(cell, j);
}

template<typename T, typename DESCRIPTOR>
void KBCBulkMomenta<T,DESCRIPTOR>::computeStress (
  Cell<T,DESCRIPTOR> const& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  lbHelpers<T,DESCRIPTOR>::computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void KBCBulkMomenta<T,DESCRIPTOR>::computeRhoU (
  Cell<T,DESCRIPTOR> const& cell,
  T& rho, T u[DESCRIPTOR::d] ) const
{
  lbHelpers<T,DESCRIPTOR>::computeRhoU(cell, rho,u);
}

template<typename T, typename DESCRIPTOR>
void KBCBulkMomenta<T,DESCRIPTOR>::computeAllMomenta (
  Cell<T,DESCRIPTOR> const& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  lbHelpers<T,DESCRIPTOR>::computeAllMomenta(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void KBCBulkMomenta<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  T oldRho, u[DESCRIPTOR::d];
  computeRhoU(cell, oldRho, u);
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  T fNeq[DESCRIPTOR::q];
  kbcLbHelpers<T,DESCRIPTOR>::computeFneq(cell, fNeq, oldRho, u);
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = kbcLbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr) +
                 fNeq[iPop];
  }
}

template<typename T, typename DESCRIPTOR>
void KBCBulkMomenta<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{
  T rho, oldU[DESCRIPTOR::d];
  computeRhoU(cell, rho, oldU);
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  T fNeq[DESCRIPTOR::q];
  kbcLbHelpers<T,DESCRIPTOR>::computeFneq(cell, fNeq, rho, oldU);
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = kbcLbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);// +
                 //fNeq[iPop];
  }

}

template<typename T, typename DESCRIPTOR>
void KBCBulkMomenta<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{
  T oldRho, oldU[DESCRIPTOR::d];
  computeRhoU(cell, oldRho, oldU);
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  T fNeq[DESCRIPTOR::q];
  kbcLbHelpers<T,DESCRIPTOR>::computeFneq(cell, fNeq, oldRho, oldU);
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = kbcLbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);// + SM - REMOVE NEQ COMPONENT TEMPORARILY
                 //fNeq[iPop];
  }
}

//Not modified for KBC yet
template<typename T, typename DESCRIPTOR>
void KBCBulkMomenta<T,DESCRIPTOR>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = kbcLbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr) +
                 firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
  }
}

/////////////// Singletons //////////////////////////////////

namespace instances {

template<typename T, typename DESCRIPTOR>
KBCBulkMomenta<T,DESCRIPTOR>& getKBCBulkMomenta()
{
  static KBCBulkMomenta<T,DESCRIPTOR> bulkMomentaSingleton;
  return bulkMomentaSingleton;
}
}


}

#endif



