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

namespace olb {

//==============================================================================//
/////////////////////////// Class KBCdynamics ///////////////////////////////
//==============================================================================//

template<typename T, typename DESCRIPTOR>
KBCdynamics<T, DESCRIPTOR>::KBCdynamics(
  T omega, Momenta<T, DESCRIPTOR>& momenta)
  : BasicDynamics<T, DESCRIPTOR>(momenta),
    _omega(omega), _beta(0.5 * omega)
{ }

template<typename T, typename DESCRIPTOR>
T KBCdynamics<T,DESCRIPTOR>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return kbcLbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
}

template<typename T, typename DESCRIPTOR>
void KBCdynamics<T, DESCRIPTOR>::collide(
  Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics)
{
  typedef kbcLbHelpers<T, DESCRIPTOR> kbcLbH;

  //Compute rho, u
  T rho, u[DESCRIPTOR::d];
  this -> _momenta.computeRhoU(cell, rho, u);
  T uSqr = kbcLbH::kbcCollision(cell, rho, u, _beta);

  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T KBCdynamics<T, DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
T KBCdynamics<T, DESCRIPTOR>::getBeta() const
{
  return _beta;
}

template<typename T, typename DESCRIPTOR>
void KBCdynamics<T, DESCRIPTOR>::setOmega(T omega) 
{
  _omega = omega;
}

//==============================================================================//
/////////////////////////// Class KBCGradDynamics ///////////////////////////////
//==============================================================================//

template<typename T, typename DESCRIPTOR>
KBCGradDynamics<T, DESCRIPTOR>::KBCGradDynamics(
  T omega, Momenta<T, DESCRIPTOR>& momenta)
  : KBCdynamics<T, DESCRIPTOR>(omega, momenta)
{ }

template<typename T, typename DESCRIPTOR>
void KBCGradDynamics<T, DESCRIPTOR>::collide(
  Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics)
{
  typedef kbcLbHelpers<T, DESCRIPTOR> kbcLbH;

  //Compute rho, u
  T rho, u[DESCRIPTOR::d];
  this -> _momenta.computeRhoU(cell, rho, u);

  //Store pre-collision velocities within velocity field
  //std::cout << "KBC collide " << u[0] << " " << u[1] << " " << u[2] << std::endl; 
  cell.template defineField<descriptors::VELOCITY>(u);
  //std::cout << "KBC collide " << u[0] << " " << u[1] << " " << u[2] << " " << cell.template getField<descriptors::VELOCITY>()[0] << " " << cell.template getField<descriptors::VELOCITY>()[1] << " " << cell.template getField<descriptors::VELOCITY>()[2] << std::endl; 
   

  T uSqr = kbcLbH::kbcCollision(cell, rho, u, KBCdynamics<T,DESCRIPTOR>::getBeta());

  statistics.incrementStats(rho, uSqr);
}


//==============================================================================//
/////////////////////////// Class KBCSpongeDynamics ///////////////////////////////
//==============================================================================//

template<typename T, typename DESCRIPTOR>
KBCSpongeDynamics<T, DESCRIPTOR>::KBCSpongeDynamics(
  T omega, Momenta<T, DESCRIPTOR>& momenta)
  : KBCdynamics<T, DESCRIPTOR>(omega, momenta)
{ }

template<typename T, typename DESCRIPTOR>
void KBCSpongeDynamics<T, DESCRIPTOR>::collide(
  Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics)
{
  typedef kbcLbHelpers<T, DESCRIPTOR> kbcLbH;

  //Compute rho, u
  T rho, u[DESCRIPTOR::d];
  this -> _momenta.computeRhoU(cell, rho, u);

  //Store pre-collision velocities within velocity field
  //std::cout << "KBC collide " << u[0] << " " << u[1] << " " << u[2] << std::endl; 
  //cell.template defineField<descriptors::TAU_EFF>(u);
  //std::cout << "KBC collide " << u[0] << " " << u[1] << " " << u[2] << " " << cell.template getField<descriptors::VELOCITY>()[0] << " " << cell.template getField<descriptors::VELOCITY>()[1] << " " << cell.template getField<descriptors::VELOCITY>()[2] << std::endl; 

  //COMPUTE EFFECTIVE RELAXATION (INCL SPONGE) AND PASS TO COLLISION 

  T tauEff = 1. / KBCdynamics<T,DESCRIPTOR>::getOmega() + cell.template getField<descriptors::TAU_EFF>();

  T betaEff = 0.5 * 1. / tauEff;

  //std::cout << "BETAEFF" << betaEff << std::endl;

  T uSqr = kbcLbH::kbcCollision(cell, rho, u, betaEff);

  statistics.incrementStats(rho, uSqr);
}

} //namespace


#endif