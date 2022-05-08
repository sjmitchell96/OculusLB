/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  Generic collision, which modifies the particle distribution
 *  functions, implemented by Orestis Malaspinas, 2007
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

#ifndef SPONGE_REGIONS_3D_HH
#define SPONGE_REGIONS_3D_HH

#include "core/postProcessing.h"
#include "core/blockLattice3D.h"
#include "spongeRegions3D.hh"

#include <cmath>

namespace olb {

template<typename T, typename DESCRIPTOR>
SpongeRegionGenerator3D<T,DESCRIPTOR>::SpongeRegionGenerator3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{  }

template<typename T, typename DESCRIPTOR>
void SpongeRegionGenerator3D<T,DESCRIPTOR>::shift(
  int deltaX, int deltaY, int deltaZ)
{
  x0 += deltaX;
  x1 += deltaX;
  y0 += deltaY;
  y1 += deltaY;
  z0 += deltaZ;
  z1 += deltaZ;
}

//Class SineSpongeRegion3D
template<typename T, typename DESCRIPTOR> 
SineSpongeRegion3D<T,DESCRIPTOR>::SineSpongeRegion3D(
    BlockGeometryStructure3D<T>& blockGeometryStructure_,
    BlockIndicatorF3D<T>&& bulkIndicator_,
    IndicatorF3D<T>&& spongeIndicator_, Vector<T, 3> orientation_,
    T tauBase_, T tauMax_ )
: blockGeometryStructure(blockGeometryStructure_), bulkIndicator(std::forward<BlockIndicatorF3D<T>>(bulkIndicator_)), 
  spongeIndicator(std::forward<IndicatorF3D<T>>(spongeIndicator_)), orientation(orientation_), tauBase(tauBase_), tauMax(tauMax_)
{
  //Decide on orientation (distanceStart) based on input
  if (orientation[0] == 1 && orientation[1] == 0 && orientation[2] == 0) {
    dir = 0;     
  //Compute start distance from sponge indicator min
    spStartPhys = spongeIndicator.getMin()[dir];
    spEndPhys = spongeIndicator.getMax()[dir];
  }
  else if (orientation[0] == 0 && orientation[1] == 1 && orientation[2] == 0) {
    dir = 1;     
    spStartPhys = spongeIndicator.getMin()[dir];
    spEndPhys = spongeIndicator.getMax()[dir];
  }
  else {
    std::cout << "SPONGE WARNING: INVALID ORIENTATION PLANE" << std::endl;
    dir = 0;
  }
  //Initialise quantities for tau_effective calc (sine formula)
  amp = 0.5 * (tauMax - tauBase);
  angFreq = M_PI / (spEndPhys - spStartPhys);
}

template<typename T, typename DESCRIPTOR> 
void SineSpongeRegion3D<T,DESCRIPTOR>::initialise(BlockLattice3D<T, DESCRIPTOR>& blockLattice)
{
  //Loop over all nodes
  int iX;
  for (iX = 0; iX <= blockLattice.getNx() - 1; ++iX) {
    for (int iY = 0; iY <= blockLattice.getNy() - 1; ++iY) {
      for (int iZ = 0; iZ <= blockLattice.getNz() - 1; ++iZ) {
        T physR[3];

        blockGeometryStructure.getPhysR(physR,iX,iY,iZ);
        int mat = blockGeometryStructure.getMaterial(iX, iY, iZ);
        bool isBulk = ( mat == 1 || mat == 2 || mat == 3 || mat ==4 );//bulkIndicator(iX, iY, iZ);
        bool isSponge = spongeIndicator(physR);

        if (isBulk && isSponge) {
          T d = physR[dir];   
          const T tauEffOld = blockLattice.get(iX,iY,iZ).template getField<descriptors::TAU_EFF>();
          //std::cout << "TAUEFFOLD" << tauEffOld << std::endl;
          T tauEff = tauEffOld + ( amp * sin(angFreq * (d - spStartPhys - 0.5 
            * (spEndPhys - spStartPhys))) + /*tauBase*/ + amp);
          if (tauBase + tauEff > tauMax) {
            tauEff = tauMax - tauBase;
          }
          //std::cout << iX << " " << iY << " " << iZ << std::endl;
          //std::cout << "TAU EFF " << tauEff << std::endl;
          blockLattice.get(iX,iY,iZ).template defineField<descriptors::TAU_EFF>(&tauEff);
          //std::cout << "TAUEFF" << tauEff + tauBase << std::endl;
        }
        //else { 
        //  const T tauEff = 0.;
          //std::cout << "TAU EFF " << tauEff << std::endl;
        //  blockLattice.get(iX,iY,iZ).template defineField<descriptors::TAU_EFF>(&tauEff);
        //}
      }
    }
  }
}

//Class SineSpongeRegionGenerator3D
template<typename T, typename DESCRIPTOR> 
SineSpongeRegionGenerator3D<T,DESCRIPTOR>::SineSpongeRegionGenerator3D(
    BlockGeometryStructure3D<T>& blockGeometryStructure_,
    BlockIndicatorF3D<T>&& bulkIndicator_,
    IndicatorF3D<T>&& spongeIndicator_, Vector<T, 3> orientation_,
    T tauBase_, T tauMax_ )
: SpongeRegionGenerator3D<T,DESCRIPTOR>(-1,-1,-1,-1,-1,-1),
  blockGeometryStructure(blockGeometryStructure_), bulkIndicator(std::forward<BlockIndicatorF3D<T>>(bulkIndicator_)), 
  spongeIndicator(std::forward<IndicatorF3D<T>>(spongeIndicator_)), orientation(orientation_), tauBase(tauBase_), tauMax(tauMax_)
{

}

template<typename T, typename DESCRIPTOR> 
SpongeRegion3D<T, DESCRIPTOR>* SineSpongeRegionGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new SineSpongeRegion3D<T,DESCRIPTOR>
    (this->blockGeometryStructure,
     std::forward<BlockIndicatorF3D<T>>(this->bulkIndicator),
     std::forward<IndicatorF3D<T>>(this->spongeIndicator),
     this->orientation, this->tauBase, this->tauMax);
}

template<typename T, typename DESCRIPTOR> 
SpongeRegionGenerator3D<T, DESCRIPTOR>* SineSpongeRegionGenerator3D<T,DESCRIPTOR>::clone() const
{
return new SineSpongeRegionGenerator3D<T,DESCRIPTOR>
    (this->blockGeometryStructure,
     std::forward<BlockIndicatorF3D<T>>(this->bulkIndicator),
     std::forward<IndicatorF3D<T>>(this->spongeIndicator),
     this->orientation, this->tauBase, this->tauMax);
}

}

#endif