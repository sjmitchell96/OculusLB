/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias J. Krause
 *  E-mail contact: info@openlb.net
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
 * A helper for initialising 3D boundaries -- header file.
 */

#ifndef SUPER_VISCOSITY_SPONGE_3D_H
#define SUPER_VISCOSITY_SPONGE_3D_H

#include <vector>
#include "io/ostreamManager.h"
#include "utilities/functorPtr.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T, typename DESCRIPTOR> class ViscositySponge3D;
template<typename T, typename DESCRIPTOR> class SuperLattice3D;
template<typename T> class SuperGeometry3D;
template<typename T> class SuperIndicatorF3D;

/// A helper for initialising 3D boundaries for super lattices.
/** Here we have methods that initializes the local postprocessors and the
 * communicator (_commBC in SuperLattice) for boundary conditions
 * for a given global point or global range.
 *
 * This class is not intended to be derived from.
 */
template<typename T, typename DESCRIPTOR>
class sViscositySponge3D {
public:
  sViscositySponge3D(SuperLattice3D<T, DESCRIPTOR>& sLattice);
  ~sViscositySponge3D();

  //Sinusoid distribution sponge - DECIDE ON PARAMS - bounds and orientation
  
  void addSineSponge(FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator,
    IndicatorF3D<T>&& spongeIndicator, Vector<T, 3> orientation, T tauBase,
    T tauMax);
  void addSineSponge(SuperGeometry3D<T>& superGeometry,
                               IndicatorF3D<T>& spongeIndicator,
                               Vector<T, 3> orientation,
                               T tauBase, T tauMax,
                               std::vector<int> bulkMaterials = std::vector<int>(1,1));

  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  void addPoints2CommBC(FunctorPtr<SuperIndicatorF3D<T>>&& indicator);
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  void addPoints2CommBC(SuperGeometry3D<T>& superGeometry, int material);

  SuperLattice3D<T, DESCRIPTOR>& getSuperLattice();
  std::vector<ViscositySponge3D<T, DESCRIPTOR>*>& getBlockSponges();
  int getOverlap();
  void setOverlap(int overlap);

  void outputOn();
  void outputOff();

private:
  mutable OstreamManager clout;
  SuperLattice3D<T, DESCRIPTOR>& _sLattice;
  std::vector<ViscositySponge3D<T, DESCRIPTOR>*> _blockSponges;
  int _overlap;
  bool _output;
};

////////////////// Factory functions //////////////////////////////////

template<typename T, typename DESCRIPTOR>
void createViscositySponge3D(sViscositySponge3D<T, DESCRIPTOR>& sVS);

} // namespace olb

#endif