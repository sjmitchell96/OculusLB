/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
 * A helper for initialising 3D boundaries -- generic implementation.
 */
#ifndef VISCOSITY_SPONGE_3D_HH
#define VISCOSITY_SPONGE_3D_HH

#include "viscositySponge3D.h"
#include "momentaOnBoundaries3D.h"

namespace olb {

////////// Convenience wrappers for boundary functions ////////////////////////

template<typename T, typename DESCRIPTOR>
ViscositySponge3D<T,DESCRIPTOR>::ViscositySponge3D (
  BlockLatticeStructure3D<T,DESCRIPTOR>& block)
  : _block(block), _output(false)
{ }

template<typename T, typename DESCRIPTOR>
void ViscositySponge3D<T, DESCRIPTOR>::addSineSponge(
  BlockIndicatorF3D<T>& bulkIndicator,
  IndicatorF3D<T>&& spongeIndicator, Vector<T, 3> orientation,
  T tauBase, T tauMax)
{

  SineSpongeRegionGenerator3D<T, DESCRIPTOR>* spongeRegion =
    new SineSpongeRegionGenerator3D<T, DESCRIPTOR>(
      bulkIndicator.getBlockGeometryStructure(),
      std::forward<BlockIndicatorF3D<T>>(bulkIndicator),
      std::forward<IndicatorF3D<T>>(spongeIndicator),
      orientation, tauBase, tauMax);
  if (spongeRegion) {
    _block.addSpongeRegion(*spongeRegion);
  }
  
}

template<typename T, typename DESCRIPTOR>
void ViscositySponge3D<T,DESCRIPTOR>::outputOn()
{
  _output = true;
}

template<typename T, typename DESCRIPTOR>
void ViscositySponge3D<T,DESCRIPTOR>::outputOff()
{
  _output = false;
}
////////// Factory functions //////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
ViscositySponge3D<T,DESCRIPTOR>* createViscositySponge3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block)
{
  return new ViscositySponge3D <
         T, DESCRIPTOR > (block);
}

}  // namespace olb

#endif