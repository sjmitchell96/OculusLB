/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Kratzke, Mathias J. Krause
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

#ifndef OFF_BOUNDARY_POST_PROCESSORS_3D_H
#define OFF_BOUNDARY_POST_PROCESSORS_3D_H

#include "core/postProcessing.h"
#include "core/blockLattice3D.h"

namespace olb {

/**
* SM - This class computes the Quadratic Bouzidi BC
*/

template<typename T, typename DESCRIPTOR>
class ZeroVelocityBouzidiQuadraticPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  ZeroVelocityBouzidiQuadraticPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int x, y, z;
  int xN, yN, zN, xB, yB, zB;
  int iPop, opp, iPop2;
  T q, dist;
};

/**
* This class computes the Linear Bouzidi BC
*/

template<typename T, typename DESCRIPTOR>
class ZeroVelocityBouzidiLinearPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  ZeroVelocityBouzidiLinearPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int x, y, z;
  int xN, yN, zN, xB, yB, zB;
  int iPop, opp, iPop2;
  T q, dist;
};

template<typename T, typename DESCRIPTOR>
class ZeroVelocityBounceBackPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  ZeroVelocityBounceBackPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;

private:
  int x, y, z;
  int xN, yN, zN;
  int iPop, opp;
  T dist;
};

template<typename T, typename DESCRIPTOR>
class VelocityBouzidiLinearPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  VelocityBouzidiLinearPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int x, y, z;
  int xN, yN, zN, xB, yB, zB;
  int iPop, opp, iPop2;
  T q, dist;
  T ufrac;
};

template<typename T, typename DESCRIPTOR>
class VelocityBounceBackPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  VelocityBounceBackPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int x, y, z;
  int xN, yN, zN;
  int iPop, opp;
  T dist;
};

//Grad BC post processor - TODO
template<typename T, typename DESCRIPTOR>
class ZeroVelocityGradPostProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  ZeroVelocityGradPostProcessor3D(int x_, int y_, int z_,
    std::vector<T> distances_, std::vector<unsigned> iMissing_, BlockGeometryStructure3D<T>& blockGeometryStructure_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int x, y, z;
  std::vector<unsigned> iNotMissing; 
  std::vector<unsigned> iClean, iDirty;
  std::vector<unsigned> xFclean, yFclean, zFclean;
  std::vector<unsigned> xFdirty, yFdirty, zFdirty;
  std::vector<T> distancesClean, distancesDirty;
  std::vector<unsigned> pStencilX = std::vector<unsigned>(2);
  std::vector<unsigned> pStencilY = std::vector<unsigned>(2);
  std::vector<unsigned> pStencilZ = std::vector<unsigned>(2);
  BlockGeometryStructure3D<T>& blockGeometryStructure;
  int nMissing;
  
  std::vector<int> getNext(int x, int y, int z);
  std::vector<int> getPosition() override;
  // computePi
  // updateMissing
  // computeUf
  //
};

//New class - dirtyNodePostProcessor, cleanNodePostProcessor

/**
* SM - Quadratic Bouzidi BC Generator
*/

template<typename T, typename DESCRIPTOR>
class ZeroVelocityBouzidiQuadraticPostProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  ZeroVelocityBouzidiQuadraticPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int x, y, z;
  int iPop;
  T dist;
};

/**
* Linear Bouzidi BC Generator
*/

template<typename T, typename DESCRIPTOR>
class ZeroVelocityBouzidiLinearPostProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  ZeroVelocityBouzidiLinearPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int x, y, z;
  int iPop;
  T dist;
};

template<typename T, typename DESCRIPTOR>
class ZeroVelocityBounceBackPostProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  ZeroVelocityBounceBackPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int x, y, z;
  int iPop;
  T dist;
};

template<typename T, typename DESCRIPTOR>
class VelocityBouzidiLinearPostProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  VelocityBouzidiLinearPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int x, y, z;
  int iPop;
  T dist;
};

template<typename T, typename DESCRIPTOR>
class VelocityBounceBackPostProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  VelocityBounceBackPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int x, y, z;
  int iPop;
  T dist;
};

/**
* Grad BC Generator - TODO
*/

template<typename T, typename DESCRIPTOR>
class ZeroVelocityGradPostProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  ZeroVelocityGradPostProcessorGenerator3D(int x_, int y_, int z_, std::vector<T> distances_,
  std::vector<unsigned> iMissing_, BlockGeometryStructure3D<T>& blockGeometryStructure_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int x, y, z;
  std::vector<T> distances;
  std::vector<unsigned> iMissing;
  BlockGeometryStructure3D<T>& blockGeometryStructure; 
};

}

#endif
