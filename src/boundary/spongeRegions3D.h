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

#ifndef SPONGE_REGIONS_3D_H
#define SPONGE_REGIONS_3D_H

#include "core/postProcessing.h"
#include "momentaOnBoundaries.h"
#include "core/blockLattice3D.h"

namespace olb {

template<typename T> class BlockIndicatorF3D;
 
//Sponge region defines sponge properties and applies change in effective relaxation time 
template<typename T, typename DESCRIPTOR>
class SpongeRegion3D{
public:
  virtual ~SpongeRegion3D() { };
  virtual void initialise(BlockLattice3D<T,DESCRIPTOR>& blockLattice)=0; 
  virtual void initialiseSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)=0;   
};

//Class SpongeRegionGenerator3D
template<typename T, typename DESCRIPTOR>
class SpongeRegionGenerator3D{
public:
  SpongeRegionGenerator3D(int x0_, int x1_, int y0_, int y1_,
    int z0_, int z1_); 
  virtual ~SpongeRegionGenerator3D() { };
  void shift(int deltaX, int deltaY, int deltaZ);
  virtual SpongeRegion3D<T, DESCRIPTOR>* generate() const =0;
  virtual SpongeRegionGenerator3D<T, DESCRIPTOR>* clone() const =0;
protected:
  int x0, x1, y0, y1, z0, z1;
};

//Class sineSponge - derived from SpongeRegion
template<typename T, typename DESCRIPTOR>
class SineSpongeRegion3D : public SpongeRegion3D<T,DESCRIPTOR> {
public:
  SineSpongeRegion3D(BlockGeometryStructure3D<T>& blockGeometryStructure_,
                     BlockIndicatorF3D<T>&& bulkIndicator_,
                     IndicatorF3D<T>&& spongeIndicator_, Vector<T, 3> orientation_,
                     T tauBase_, T tauMax_);
  void initialise(BlockLattice3D<T, DESCRIPTOR>& blockLattice) override; 
  void initialiseSubDomain(BlockLattice3D<T, DESCRIPTOR>& blockLattice,
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override
    {
        this->initialise(blockLattice);
    }
private:
  BlockGeometryStructure3D<T>& blockGeometryStructure;
  BlockIndicatorF3D<T>&& bulkIndicator;
  IndicatorF3D<T>&& spongeIndicator;
  Vector<T, 3> orientation;
  T tauBase, tauMax;
  int dir;
  T spStartPhys, spEndPhys, amp, angFreq; 
};

//class sineSpongeGenerator - derived from SpongeRegionGenerator
template<typename T, typename DESCRIPTOR>
class SineSpongeRegionGenerator3D : public SpongeRegionGenerator3D<T,DESCRIPTOR>{
public:
  SineSpongeRegionGenerator3D(BlockGeometryStructure3D<T>& blockGeometryStructure_,
                              BlockIndicatorF3D<T>&& bulkIndicator_,
                              IndicatorF3D<T>&& spongeIndicator_, Vector<T, 3> orientation_,
                              T tauBase_, T tauMax_);
  SpongeRegion3D<T, DESCRIPTOR>* generate() const override;
  SpongeRegionGenerator3D<T, DESCRIPTOR>* clone() const override;
private:
  BlockGeometryStructure3D<T>& blockGeometryStructure;
  BlockIndicatorF3D<T>&& bulkIndicator;
  IndicatorF3D<T>&& spongeIndicator;
  Vector<T, 3> orientation;
  T tauBase, tauMax; 
};

}

#endif