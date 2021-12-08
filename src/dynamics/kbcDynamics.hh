/*  This file is part of the
 *
 *  
 *  E-mail contact
 *  The most recent release of OpenLB can be downloaded at
 *  
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
 * A collection of entropic dynamics classes (e.g. Entropic, ForcedEntropic, Entropic, ForcedEntropic) with which a Cell object
 * can be instantiated -- generic implementation.
 */

#ifndef KBC_DYNAMICS_HH
#define KBC_DYNAMICS_HH

#include "lbHelpers.h"
#include "kbcDynamics.h"
#include "kbcLbHelpers.h"

//==============================================================================//
/////////////////////////// Class KBCdynamics ///////////////////////////////
//==============================================================================//

template<typename T, typename DESCRIPTOR>
KBCdynamics<T, DESCRIPTOR>::KBCdynamics(
  T beta_, Momenta<T, DESCRIPTOR>& momenta)
  : BasicDynamics<T, DESCRIPTOR>(momenta),
    _beta(beta)
{ }

template<typename T, typename DESCRIPTOR>
T KBCdynamics<T, DESCRIPTOR>::computeEquilibrium(
  int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return kbcLbHelpers<T, DESCRIPTOR>::equilibrium(iPop, rho, u);
}

template<typename T, typename DESCRIPTOR>
void KBCdynamics<T, DESCRIPTOR>::collide(
  Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics)
{
  typedef kbcLbHelpers<T, DESCRIPTOR> kbcLbH;

  //Compute rho, u
  T rho, u[DESCRIPTOR::d];
  *this -> _momenta.computeRhoU(cell, rho, u);
  T uSqr = kbcLbh::kbcCollision(cell, rho, u, _beta);

  //In helper instead ...
  //T uSqr = util::normSqr<T,L::d>(u);
  //Compute feq
  //T f[L::q], feq[L::q], fNeq[L::q];
  //for(int iPop = 0; iPop < L::q; ++iPop) {
  //    fEq[iPop] = kbcLbH::equilibrium(iPop, rho, u);
  //    fNeq[iPop] = cell[iPop] - fEq[iPop];
  //}
  //Compute delta si
  //Compute delta hi
  //Compute stabiliser gamma
  //Relax
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T KBCdynamics<T, DESCRIPTOR>::getBeta() const
{
  return _beta;
}

template<typename T, typename DESCRIPTOR>
void KBCdynamics<T, DESCRIPTOR>::setBeta(T beta) 
{
  _beta = beta;
}

#endif