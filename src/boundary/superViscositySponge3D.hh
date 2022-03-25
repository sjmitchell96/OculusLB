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
 * A helper for initialising 3D boundaries -- generic implementation.
 */

#ifndef SUPER_VISCOSITY_SPONGE_3D_HH
#define SUPER_VISCOSITY_SPONGE_3D_HH

#include <vector>
#include "boundaryCondition3D.h"
#include "geometry/superGeometry3D.h"
#include "superViscositySponge3D.h"
#include "core/superLattice3D.h"
#include "functors/lattice/indicator/superIndicatorBaseF3D.h"

namespace olb {

///////// class sViscositySponge3D ///////////////////////////////

template<typename T, typename DESCRIPTOR>
sViscositySponge3D<T, DESCRIPTOR>::sViscositySponge3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice) :
  clout(std::cout,"sViscositySponge3D"),
  _sLattice(sLattice),
  _output(false)
{
}

template<typename T, typename DESCRIPTOR>
sViscositySponge3D<T, DESCRIPTOR>::~sViscositySponge3D()
{
  for (auto &iC : _blockSponges) {
    delete iC;
  }
}

template<typename T, typename DESCRIPTOR>
void sViscositySponge3D<T,DESCRIPTOR>::addSineSponge(
  FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator,
  IndicatorF3D<T>&& spongeIndicator,
  Vector<T, 3> orientation,
  T tauBase, T tauMax)
{
  if (_output) {
    clout.setMultiOutput(true);
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    if (_output) {
      clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
            << std::endl;
    }
    _blockSponges[iCloc]->addSineSponge(
    bulkIndicator->getExtendedBlockIndicatorF(iCloc),
    std::forward<IndicatorF3D<T>>(spongeIndicator), orientation, tauBase, tauMax);
  }
  if (_output) {
    clout.setMultiOutput(false);
  }
}

template<typename T, typename DESCRIPTOR>
void sViscositySponge3D<T,DESCRIPTOR>::addSineSponge(
  SuperGeometry3D<T>& superGeometry,
  IndicatorF3D<T>& spongeIndicator,
  Vector<T, 3> orientation,
  T tauBase, T tauMax,
  std::vector<int> bulkMaterials)
{
  addSineSponge(
    std::move(superGeometry.getMaterialIndicator(std::move(bulkMaterials))),
    std::move(spongeIndicator),
    orientation, tauBase, tauMax);
}

template<typename T, typename DESCRIPTOR>
void sViscositySponge3D<T, DESCRIPTOR>::addPoints2CommBC(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  if (_overlap == 0) {
    return;
  }

  SuperGeometry3D<T>& superGeometry = indicator->getSuperGeometry();
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    const int nX = superGeometry.getBlockGeometry(iCloc).getNx();
    const int nY = superGeometry.getBlockGeometry(iCloc).getNy();
    const int nZ = superGeometry.getBlockGeometry(iCloc).getNz();

    for (int iX = -_overlap; iX < nX+_overlap; ++iX) {
      for (int iY = -_overlap; iY < nY+_overlap; ++iY) {
        for (int iZ = -_overlap; iZ < nZ+_overlap; ++iZ) {
          if (iX < 0 || iX > nX - 1 ||
              iY < 0 || iY > nY - 1 ||
              iZ < 0 || iZ > nZ - 1 ) { // if within overlap
            if (superGeometry.getBlockGeometry(iCloc).getMaterial(iX,iY,iZ) != 0) {
              bool found = false;
              for (int iXo = -_overlap; iXo <= _overlap && !found; ++iXo) {
                for (int iYo = -_overlap; iYo <= _overlap && !found; ++iYo) {
                  for (int iZo = -_overlap; iZo <= _overlap && !found; ++iZo) {
                    const int nextX = iXo + iX;
                    const int nextY = iYo + iY;
                    const int nextZ = iZo + iZ;
                    if (indicator->getBlockIndicatorF(iCloc)(nextX, nextY, nextZ)
                        && nextX >= -_overlap && nextX < nX+_overlap
                        && nextY >= -_overlap && nextY < nY+_overlap
                        && nextZ >= -_overlap && nextZ < nZ+_overlap) {
                      _sLattice.get_commBC().add_cell(iCloc, iX, iY, iZ);
                      found = true;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void sViscositySponge3D<T, DESCRIPTOR>::addPoints2CommBC(SuperGeometry3D<T>& superGeometry, int material)
{
  addPoints2CommBC(superGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
SuperLattice3D<T, DESCRIPTOR>& sViscositySponge3D<T, DESCRIPTOR>::getSuperLattice()
{
  return _sLattice;
}

template<typename T, typename DESCRIPTOR>
std::vector<ViscositySponge3D<T, DESCRIPTOR>*>& sViscositySponge3D<T, DESCRIPTOR>::getBlockSponges()
{
  return _blockSponges;
}

template<typename T, typename DESCRIPTOR>
int sViscositySponge3D<T, DESCRIPTOR>::getOverlap()
{
  return _overlap;
}

template<typename T, typename DESCRIPTOR>
void sViscositySponge3D<T, DESCRIPTOR>::setOverlap(int overlap)
{
  _overlap = overlap;
}

//////////////// Output functions //////////////////////////////////
template<typename T, typename DESCRIPTOR>
void sViscositySponge3D<T, DESCRIPTOR>::outputOn()
{
  _output = true;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockSponges[iCloc]->outputOn();
  }
}

template<typename T, typename DESCRIPTOR>
void sViscositySponge3D<T, DESCRIPTOR>::outputOff()
{
  _output = false;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockSponges[iCloc]->outputOff();
  }
}


////////////////// Factory functions //////////////////////////////////

template<typename T, typename DESCRIPTOR>
void createViscositySponge3D(sViscositySponge3D<T, DESCRIPTOR>& sVS)
{
  int nC = sVS.getSuperLattice().getLoadBalancer().size();
  sVS.setOverlap(0);
  for (int iC = 0; iC < nC; iC++) {
    ViscositySponge3D<T, DESCRIPTOR>* blockSponge =
      createViscositySponge3D(
        sVS.getSuperLattice().getExtendedBlockLattice(iC));
    sVS.getBlockSponges().push_back(blockSponge);
  }
}

} // namespace olb

#endif