/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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

#ifndef VISCOSITY_SPONGE_3D_H
#define VISCOSITY_SPONGE_3D_H

#include "dynamics/dynamics.h"
#include "core/unitConverter.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
class ViscositySponge3D {
public:
  ViscositySponge3D( BlockLatticeStructure3D<T,DESCRIPTOR>& block_ );
  ~ViscositySponge3D() { }

  //DECIDE ON PARAMS
  void addSineSponge(BlockIndicatorF3D<T>& bulkIndicator,
    IndicatorF3D<T>&& spongeIndicator, Vector<T, 3> orientation,
    T tauBase, T tauMax);

  void outputOn();
  void outputOff();

private:
  BlockLatticeStructure3D<T,DESCRIPTOR>& _block;
  bool _output;
};

////////// Factory functions //////////////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
ViscositySponge3D<T,DESCRIPTOR>* createViscositySponge3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block);
}

#endif