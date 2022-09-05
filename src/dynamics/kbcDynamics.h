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
 *
 * Entropic Multirelxation time - Karlin, Boesch, Chikatamarla (KBC) model:
 * Fabian Boesch, Shyam S. Chitakamarla, and Iliya V. Karlin
 * Entropic multirelaxation time lattice Boltzmann models for turbulent flows
 * Physical Review E 92, 043309 (2015)
 */

#ifndef KBC_DYNAMICS_H
#define KBC_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {

template <typename T, typename DESCRIPTOR > class Cell;

//Implementation of KBC collision
template <typename T, typename DESCRIPTOR> 
class KBCdynamics : public BasicDynamics<T, DESCRIPTOR> {
public:
  /// Constructor
  KBCdynamics(T omega, Momenta<T, DESCRIPTOR>& momenta);
  /// Compute equilibrium distribution 
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision
  void collide(Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// Get local relaxation parameter 
  T getOmega() const override;
  /// Get local relaxation parameter
  T getBeta() const;
  /// Set local relaxation parameter
  void setOmega(T omega) override;
private:
  //Some methods here for computing other collision parameters?

  T _omega; ///Relaxation parameter
  T _beta; //Relaxation parameter in EMRT form
};

//Implementation of KBC dynamics with compatibility for Grad's boundary
template <typename T, typename DESCRIPTOR> 
class KBCGradDynamics : public KBCdynamics<T, DESCRIPTOR> {
public:
  /// Constructor
  KBCGradDynamics(T omega, Momenta<T, DESCRIPTOR>& momenta);
  /// Collision
  void collide(Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
};

//Implementation of KBC dynamics with viscosity sponge layer at outlet
template <typename T, typename DESCRIPTOR> 
class KBCSpongeDynamics : public KBCdynamics<T, DESCRIPTOR> {
public:
  /// Constructor incl. sponge layer
  KBCSpongeDynamics(T omega, Momenta<T, DESCRIPTOR>& momenta);
  /// Collision
  void collide(Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
};

//Implementation of KBC dynamics with grad boundary and viscosity sponge layer at outlet
template <typename T, typename DESCRIPTOR> 
class KBCGradSpongeDynamics : public KBCdynamics<T, DESCRIPTOR> {
public:
  /// Constructor incl. sponge layer
  KBCGradSpongeDynamics(T omega, Momenta<T, DESCRIPTOR>& momenta);
  /// Collision
  void collide(Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
};

}
#endif