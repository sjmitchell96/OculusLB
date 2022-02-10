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

/** \file
 * A helper for initialising 3D boundaries -- header file.
 */

#ifndef OFF_BOUNDARY_INSTANTIATOR_3D_H
#define OFF_BOUNDARY_INSTANTIATOR_3D_H

#include "offBoundaryCondition3D.h"
#include "geometry/blockGeometry3D.h"
#include "geometry/blockGeometryStatistics3D.h"
#include "core/cell.h"
#include "io/ostreamManager.h"
#include "core/util.h"
#include "io/stlReader.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"

namespace olb {

/**
* This class gets the needed processors from BoundaryManager and adds them
* to the Processor Vector of the DESCRIPTOR
*/

template<typename T, typename DESCRIPTOR, class BoundaryManager>
class OffBoundaryConditionInstantiator3D: public OffLatticeBoundaryCondition3D<T,
  DESCRIPTOR> {
public:
  OffBoundaryConditionInstantiator3D(BlockLatticeStructure3D<T, DESCRIPTOR>& block_, T epsFraction_ = 0.0001);
  ~OffBoundaryConditionInstantiator3D() override;

  void addOnePointZeroVelocityBoundary(int x, int y, int z, int iPop, T dist) override;
  void addTwoPointZeroVelocityBoundary(int x, int y, int z, int iPop, T dist) override;
  void addMultiPointZeroVelocityBoundary(int x, int y, int z, std::vector<T> distances,
                                         std::vector<unsigned> iMissing, BlockGeometryStructure3D<T>& blockGeometryStructure) override;
  void addOnePointVelocityBoundary(int x, int y, int z, int iPop, T dist) override;
  void addTwoPointVelocityBoundary(int x, int y, int z, int iPop, T dist) override;

  void addOffDynamics(int x, int y, int z, T location[DESCRIPTOR::d]) override;
  void addOffDynamics(int x, int y, int z, T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q]) override;
  void addOffDynamics(BlockIndicatorF3D<T>& indicator) override;

  void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist) override;
  void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, T distances[DESCRIPTOR::q]);
  void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ, IndicatorF3D<T>& geometryIndicator, BlockIndicatorF3D<T>& bulkIndicator);
  void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ, IndicatorF3D<T>& geometryIndicator, std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addZeroVelocityBoundary(BlockIndicatorF3D<T>& boundaryIndicator, BlockIndicatorF3D<T>& bulkIndicator, IndicatorF3D<T>& geometryIndicator) override;

  //void addZeroVelocityGradBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist, std::vector<T> distances, int nLinks, std::vector<int> iLinks, std::vector<int> iBulk) override;
  void addZeroVelocityGradBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int xB, int yB, int zB, std::vector<T> distances, std::vector<unsigned> iMissing);
  void addZeroVelocityGradBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ, IndicatorF3D<T>& geometryIndicator, BlockIndicatorF3D<T>& bulkIndicator);
  void addZeroVelocityGradBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ, IndicatorF3D<T>& geometryIndicator, std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addZeroVelocityGradBoundary(BlockIndicatorF3D<T>& boundaryIndicator, BlockIndicatorF3D<T>& bulkIndicator, IndicatorF3D<T>& geometryIndicator) override;

  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist);
  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, T distances[DESCRIPTOR::q]);
  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ, IndicatorF3D<T>& geometryIndicator, BlockIndicatorF3D<T>& bulkIndicator);
  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ, IndicatorF3D<T>& geometryIndicator, std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addVelocityBoundary(BlockIndicatorF3D<T>& boundaryIndicator, BlockIndicatorF3D<T>& bulkIndicator, IndicatorF3D<T>& geometryIndicator) override;

  void setBoundaryIntersection(int iX, int iY, int iZ, int iPop, T distance);
  bool getBoundaryIntersection(int iX, int iY, int iZ, int iPop, T point[DESCRIPTOR::d]);

  void defineU(int iX, int iY, int iZ, int iPop, const T u[DESCRIPTOR::d]) override;
  void defineU(BlockIndicatorF3D<T>& indicator, BlockIndicatorF3D<T>& bulkIndicator, AnalyticalF3D<T,T>& u) override;

  void outputOn() override;
  void outputOff() override;

  BlockLatticeStructure3D<T, DESCRIPTOR>& getBlock() override;
  BlockLatticeStructure3D<T, DESCRIPTOR> const& getBlock() const override;
private:
  BlockLatticeStructure3D<T, DESCRIPTOR>& block;
  //std::vector<Momenta<T, DESCRIPTOR>*> momentaVector;
  std::vector<Dynamics<T, DESCRIPTOR>*> dynamicsVector;
  bool _output;
  mutable OstreamManager clout;
};

///////// class OffBoundaryConditionInstantiator3D ////////////////////////

template<typename T, typename DESCRIPTOR, class BoundaryManager>
BlockLatticeStructure3D<T, DESCRIPTOR>& OffBoundaryConditionInstantiator3D<T, DESCRIPTOR,
                        BoundaryManager>::getBlock()
{
  return block;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
BlockLatticeStructure3D<T, DESCRIPTOR> const& OffBoundaryConditionInstantiator3D<T, DESCRIPTOR,
                        BoundaryManager>::getBlock() const
{
  return block;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::OffBoundaryConditionInstantiator3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& block_, T epsFraction_) :
  block(block_), _output(false), clout(std::cout,"BoundaryConditionInstantiator3D")
{
  this->_epsFraction = epsFraction_;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::~OffBoundaryConditionInstantiator3D()
{
  for (unsigned iDynamics = 0; iDynamics < dynamicsVector.size(); ++iDynamics) {
    delete dynamicsVector[iDynamics];
  }
  /*
  for (unsigned iMomenta = 0; iMomenta < dynamicsVector.size(); ++iMomenta) {
    delete momentaVector[iMomenta];
  }*/
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addOnePointZeroVelocityBoundary(
  int x, int y, int z, int iPop, T dist)
{
  PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::getOnePointZeroVelocityBoundaryProcessor
    (x, y, z, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addTwoPointZeroVelocityBoundary(
  int x, int y, int z, int iPop, T dist)
{
  PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::getTwoPointZeroVelocityBoundaryProcessor
    (x, y, z, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addMultiPointZeroVelocityBoundary(
  int x, int y, int z, std::vector<T> distances, std::vector<unsigned> iMissing, BlockGeometryStructure3D<T>& blockGeometryStructure)
{
  PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::getMultiPointZeroVelocityBoundaryProcessor
    (x, y, z, distances, iMissing, blockGeometryStructure);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addOnePointVelocityBoundary(
  int x, int y, int z, int iPop, T dist)
{
  PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::getOnePointVelocityBoundaryProcessor
    (x, y, z, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addTwoPointVelocityBoundary(
  int x, int y, int z, int iPop, T dist)
{
  PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::getTwoPointVelocityBoundaryProcessor
    (x, y, z, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addOffDynamics(
  int x, int y, int z, T location[DESCRIPTOR::d])
{
  Dynamics<T,DESCRIPTOR>* dynamics
    = BoundaryManager::getOffDynamics(location);
  this->getBlock().defineDynamics(x,x,y,y,z,z, dynamics);
  dynamicsVector.push_back(dynamics);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addOffDynamics(
  int x, int y, int z, T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q])
{
  Dynamics<T,DESCRIPTOR>* dynamics
    = BoundaryManager::getOffDynamics(location, distances);
  this->getBlock().defineDynamics(x,x,y,y,z,z, dynamics);
  dynamicsVector.push_back(dynamics);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addOffDynamics(
  BlockIndicatorF3D<T>& indicator)
{
  if ( !indicator.isEmpty() ) {
    const Vector<int,3> min = indicator.getMin();
    const Vector<int,3> max = indicator.getMax();

    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        for (int iZ = min[2]; iZ <= max[2]; ++iZ) {
          if (indicator(iX,iY,iZ)) {
            T location[DESCRIPTOR::d];
            indicator.getBlockGeometryStructure().getPhysR(location, iX,iY,iZ);
            addOffDynamics(iX, iY, iZ, location);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist)
{
  const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
  if (blockGeometryStructure.getMaterial(x-c[0], y-c[1], z-c[2]) != 1) {
    addOnePointZeroVelocityBoundary(x, y, z, iPop, dist);
  }
  else {
    addTwoPointZeroVelocityBoundary(x, y, z, iPop, dist);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::
addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, T distances[DESCRIPTOR::q])
{
  typedef DESCRIPTOR L;
  //T location[DESCRIPTOR::d];
  //location[0] = blockGeometryStructure.physCoordX(x);
  //location[1] = blockGeometryStructure.physCoordY(y);
  //location[2] = blockGeometryStructure.physCoordZ(z);
  //T distancesCopy[L::q];
  //T spacing = blockGeometryStructure.getDeltaR();
  //for (int iPop = 1; iPop < L::q ; ++iPop) {
  //  distancesCopy[iPop] = spacing*(1.-distances[iPop]);
  //  if (distances[iPop] == -1)
  //    distancesCopy[iPop] = -1;
  //}
  //addOffDynamics(x, y, z, location, distancesCopy);

  for (int iPop = 1; iPop < L::q ; ++iPop) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
      addZeroVelocityBoundary(blockGeometryStructure, x-c[0], y-c[1], z-c[2], iPop, distances[iPop]);
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ,
  IndicatorF3D<T>& geometryIndicator, BlockIndicatorF3D<T>& bulkIndicator)
{
  T distances[DESCRIPTOR::q];
  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    distances[iPop] = -1;
  }

  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
    const int iXn = iX + c[0];
    const int iYn = iY + c[1];
    const int iZn = iZ + c[2];
    if (blockGeometryStructure.isInside(iXn,iYn,iZn) && bulkIndicator(iXn,iYn,iZn)) {
      T dist = -1;
      T physR[3];
      blockGeometryStructure.getPhysR(physR,iXn,iYn,iZn);
      T voxelSize=blockGeometryStructure.getDeltaR();

      Vector<T,3> physC(physR);

      Vector<T,3> direction(-voxelSize*c[0],-voxelSize*c[1],-voxelSize*c[2]);
      T cPhysNorm = voxelSize*sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

      if (!geometryIndicator.distance(dist,physC,direction,blockGeometryStructure.getIcGlob() ) ) {
        T epsX = voxelSize*c[0]*this->_epsFraction;
        T epsY = voxelSize*c[1]*this->_epsFraction;
        T epsZ = voxelSize*c[2]*this->_epsFraction;

        Vector<T,3> physC2(physC);
        physC2[0] += epsX;
        physC2[1] += epsY;
        physC2[2] += epsZ;
        Vector<T,3> direction2(direction);
        direction2[0] -= 2.*epsX;
        direction2[1] -= 2.*epsY;
        direction2[2] -= 2.*epsZ;

        if ( !geometryIndicator.distance(dist,physC2,direction2,blockGeometryStructure.getIcGlob())) {
          clout << "ERROR: no boundary found at (" << iXn << "," << iYn << "," << iZn <<") ~ ("
                << physR[0] << "," << physR[1] << "," << physR[2] <<"), "
                << "in direction " << util::opposite<DESCRIPTOR >(iPop)
                << std::endl;
        }
        T distNew = (dist - sqrt(epsX*epsX+epsY*epsY+epsZ*epsZ))/cPhysNorm;
        if (distNew < 0.5) {
          dist = 0;
        }
        else {
          dist = 0.5 * cPhysNorm;
          clout << "WARNING: distance at (" << iXn << "," << iYn << "," << iZn <<") ~ ("
                << physR[0] << "," << physR[1] << "," << physR[2] <<"), "
                << "in direction " << util::opposite<DESCRIPTOR >(iPop) << ": "
                << distNew
                << " rounded to "
                << dist/cPhysNorm
                << std::endl;
        }
      }
      distances[util::opposite<DESCRIPTOR >(iPop)] = dist/cPhysNorm;
    } // bulkMaterials if
  } // iPop loop
  addZeroVelocityBoundary(blockGeometryStructure, iX, iY, iZ, distances);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ,
  IndicatorF3D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial3D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  addZeroVelocityBoundary(blockGeometryStructure, iX, iY, iZ,
                          geometryIndicator,
                          bulkIndicator);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityBoundary(
  BlockIndicatorF3D<T>& boundaryIndicator, BlockIndicatorF3D<T>& bulkIndicator, IndicatorF3D<T>& geometryIndicator)
{
  if ( !boundaryIndicator.isEmpty() ) {
    const Vector<int,3> min = boundaryIndicator.getMin();
    const Vector<int,3> max = boundaryIndicator.getMax();

    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        for (int iZ = min[2]; iZ <= max[2]; ++iZ) {
          if (boundaryIndicator(iX,iY,iZ)) {
            addZeroVelocityBoundary(boundaryIndicator.getBlockGeometryStructure(), iX, iY, iZ,
                                    geometryIndicator, bulkIndicator);
          }
        }
      }
    }
  }
}

/* OLD -------------------------------------------

//Zero velocity Grad boundary
template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityGradBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist, std::vector<T> distances, std::vector<int> iLinks) //CHANGE TO ACCEPT ENTIRE DISTANCES VECTOR AS WELL AS SINGLE LINK DISTANCE
{
  //const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
  //T pos[3];
  //blockGeometryStructure.getPhysR(pos,x - c[0],y - c[1],z - c[2]);
  //blockGeometryStructure.getPhysR(pos,x,y,z);
  //int material = blockGeometryStructure.getMaterial(x - c[0],y - c[1],z - c[2]);
  //if (blockGeometryStructure.getMaterial(x-c[0], y-c[1], z-c[2]) != 1) { 
  //  std::cout << "ADD BOUNCEBACK" << " " << x+c[0] << " " << y+c[1] << " " << z+c[2] << " " << iPop << " " << dist << " " << blockGeometryStructure.getMaterial(x-c[0], y-c[1], z-c[2]) << " " << pos[2] << " " << material << std::endl;
  //  addOnePointZeroVelocityBoundary(x, y, z, iPop, dist);
  //}
  //else {
    //std::cout << "ADD MULTI" << " " << x+c[0] << " " << y+c[1] << " " << z+c[2] << " " << iPop << " " << dist << " " << blockGeometryStructure.getMaterial(x-c[0], y-c[1], z-c[2]) << " " << pos[2] << " " << material << std::endl;
    //std::cout << "ADD MULTI" << std::endl;
    //addMultiPointZeroVelocityBoundary(x, y, z, iPop, dist, distances);
  //}

  //std::cout << "AddZVGB " << x << " " << y << " " << z << endl;

  bool isPeriodic = true;
  //std::vector<int> xF(nLinks);
  //std::vector<int> yF(nLinks);
  //std::vector<int> zF(nLinks);

  //std::vector<int> pStencilX(2);
  //std::vector<int> pStencilY(2);
  //std::vector<int> pStencilZ(2);
  
  int xs[3];
  int jPop;

  for (int i=0; i<iLinks.size(); ++i) {
    jPop = iLinks[i];
    Vector<int,3> c = descriptors::c<DESCRIPTOR>(jPop);
    xs[0] = x + c[0];
    xs[1] = y + c[1];
    xs[2] = z + c[2];

    //xs is ghost node
    if (isPeriodic && (
        xs[0] <= 0 || xs[0] >= blockGeometryStructure.getNx() - 1 ||
        xs[1] <= 0 || xs[1] >= blockGeometryStructure.getNy() - 1 ||
        xs[2] <= 0 || xs[2] >= blockGeometryStructure.getNz() - 1 )) 
      addOnePointZeroVelocityBoundary(x, y, z, iPop, dist);
    else {
      int nextX = x - descriptors::c<DESCRIPTOR>(jPop,0);
      int nextY = y - descriptors::c<DESCRIPTOR>(jPop,1);
      int nextZ = z - descriptors::c<DESCRIPTOR>(jPop,2);

      //std::cout << jPop << std::endl;
      //std::cout << blockGeometryStructure.getNx() - 1 << " " << blockGeometryStructure.getNy() - 1 << " " << blockGeometryStructure.getNz() - 1 << std::endl;
      //std::cout << descriptors::c<DESCRIPTOR>(jPop,0) << " " << descriptors::c<DESCRIPTOR>(jPop,1) << " " << descriptors::c<DESCRIPTOR>(jPop,2) << std::endl;
      //std::cout << xs[0] << " " << xs[1] << " " << xs[2] << std::endl;
      //std::cout << x << " " << y << " " << z << std::endl;

      //xf out of lattice domain
      if (isPeriodic && (nextX < 0 || nextX >= blockGeometryStructure.getNx() ||
          nextY < 0 || nextY >= blockGeometryStructure.getNy() ||
          nextZ < 0 || nextZ >= blockGeometryStructure.getNz())) {
        nextX = (nextX + blockGeometryStructure.getNx()) % blockGeometryStructure.getNx(); 
        nextY = (nextY + blockGeometryStructure.getNy()) % blockGeometryStructure.getNy(); 
        nextZ = (nextZ + blockGeometryStructure.getNz()) % blockGeometryStructure.getNz(); 
      }

      xF[i] = nextX;
      yF[i] = nextY;
      zF[i] = nextZ;
    
      //Determine stencil points for pressure tensor gradients
      if (blockGeometryStructure.getMaterial(x + 1, y, z) == 1) {
        pStencilX[0] = x + 1;
        pStencilX[1] = x;
      }
      else if (blockGeometryStructure.getMaterial(x - 1, y, z) == 1) {
        pStencilX[0] = x;
        pStencilX[1] = x - 1;
      }
      else 
        std::cout << "GradBC Warning: No suitable stencil point in x" << std::endl;

      if (blockGeometryStructure.getMaterial(x, y + 1, z) == 1) {
        pStencilY[0] = y + 1;
        pStencilY[1] = y;
      }
      else if (blockGeometryStructure.getMaterial(x, y - 1, z) == 1) {
        pStencilY[0] = y;
        pStencilY[1] = y - 1;
      }
      else 
        std::cout << "GradBC Warning: No suitable stencil point in y" << std::endl;

      if (blockGeometryStructure.getMaterial(x, y, z + 1) == 1) {
        pStencilZ[0] = z + 1;
        pStencilZ[1] = z;
      }
      else if (blockGeometryStructure.getMaterial(x, y, z - 1) == 1) {
        pStencilZ[0] = z;
        pStencilZ[1] = z - 1;
      }
      else 
        std::cout << "GradBC Warning: No suitable stencil point in z" << std::endl;

      
      //std::cout << "STENCIL" << pStencilX[0] << " "  << pStencilX[1] << " "  << pStencilY[0] << " "  << pStencilY[1] << " "  << pStencilZ[0] << " "  << pStencilZ[1] << std::endl; 
      addMultiPointZeroVelocityBoundary(x, y, z, jPop, dist, distances,
                                        xF, yF, zF,
                                        pStencilX, pStencilY, pStencilZ,
                                        nLinks, iLinks, iBulk);
    }
  }                                  
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::
addZeroVelocityGradBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, std::vector<T> distances)
{
  typedef DESCRIPTOR L;
  //T location[DESCRIPTOR::d];
  //location[0] = blockGeometryStructure.physCoordX(x);
  //location[1] = blockGeometryStructure.physCoordY(y);
  //location[2] = blockGeometryStructure.physCoordZ(z);
  //T distancesCopy[L::q];
  //T spacing = blockGeometryStructure.getDeltaR();
  //for (int iPop = 1; iPop < L::q ; ++iPop) {
  //  distancesCopy[iPop] = spacing*(1.-distances[iPop]);
  //  if (distances[iPop] == -1)
  //    distancesCopy[iPop] = -1;
  //}
  //addOffDynamics(x, y, z, location, distancesCopy);

  //Dummy run to record number of link velocities
  //int nLinks = 0;
  std::vector<int> iLinks;
  //std::vector<int> iBulk = {0};

  //construct std::vector of distances for safe passing by value
  std::vector<T> distancesVector;

  for (int iPop = 1; iPop < L::q ; ++iPop) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      nLinks++;
      iLinks.push_back(iPop);
      distancesVector.push_back(distances[iPop]);
    }
  //  else 
  //  iBulk.push_back(iPop);
  //}
  for (int iPop = 1; iPop < L::q ; ++iPop) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
      addZeroVelocityGradBoundary(blockGeometryStructure, x-c[0], y-c[1], z-c[2], iPop, distances[iPop], distancesVector, iLinks); 
    }
  }
}

*/


//CURRENTLY HERE - JUST FINISHED ADDZVGB FUNCTIONS, 
//NEXT - MAKE SURE COMPILES  UP TO THIS POINT??, CHECK NEXT FUNCTIONS (ADD MULTI POINT), UPDATE POST PROCESSOR, GET IT WORKING??

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityGradBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int xB, int yB, int zB, std::vector<T> distances, std::vector<unsigned> iMissing)
{
 int xs[3];
 int iPop;
 bool isPeriodic = true;
 //bool isFluid;
 Vector<int, 3> cf;

  //Determine whether all xF nodes are fluid
  //for (unsigned i=0; i<iMissing.size(); ++i) {
  //  cf = descriptors::c<DESCRIPTOR>(iMissing[i]);
  //  if (blockGeometryStructure.getMaterial(xB + cf[0], yB + cf[1], zB + cf[2]) != 1)
  //    isFluid = false;
  //}

  for (unsigned i=0; i<iMissing.size(); ++i) {
    iPop = iMissing[i];
    Vector<int,3> cs = descriptors::c<DESCRIPTOR>(util::opposite<DESCRIPTOR>(iPop));
    xs[0] = xB + cs[0];
    xs[1] = yB + cs[1];
    xs[2] = zB + cs[2];

    //xs is ghost node
    if (isPeriodic && (
      xs[0] <= 0 || xs[0] >= blockGeometryStructure.getNx() - 1 ||
      xs[1] <= 0 || xs[1] >= blockGeometryStructure.getNy() - 1 ||
      xs[2] <= 0 || xs[2] >= blockGeometryStructure.getNz() - 1 )) 
        addOnePointZeroVelocityBoundary(xB, yB, zB, util::opposite<DESCRIPTOR>(iPop), distances[i]);
    //else if (!isFluid){
    //  std::cout << "NOT FLUID " << zB << std::endl; 
    //    addOnePointZeroVelocityBoundary(xB, yB, zB, util::opposite<DESCRIPTOR>(iPop), distances[i]);
    //}
    else {
       //std::cout << xB << " " << yB << " " << zB << std::endl;
       //std::cout << xs[0] << " " << xs[1] << " " << xs[2] << std::endl;
       addMultiPointZeroVelocityBoundary(xB, yB, zB, distances, iMissing, blockGeometryStructure);     
    }
  }
} 

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityGradBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ,
  IndicatorF3D<T>& geometryIndicator, BlockIndicatorF3D<T>& bulkIndicator)
{  
  
  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    const Vector<int,3> cb = descriptors::c<DESCRIPTOR>(iPop);
    const int iXb = iX + cb[0]; //Boundary node
    const int iYb = iY + cb[1];
    const int iZb = iZ + cb[2];

    if (blockGeometryStructure.isInside(iXb,iYb,iZb) && bulkIndicator(iXb, iYb, iZb)) {

      T physR[3];
      blockGeometryStructure.getPhysR(physR,iXb,iYb,iZb);
      T voxelSize=blockGeometryStructure.getDeltaR();
      Vector<T,3> physC(physR);
      std::vector<T> distances;
      std::vector<unsigned> iMissing;

      for (int jPop = 1; jPop < DESCRIPTOR::q; ++jPop) {
        const Vector<int,3> cs = descriptors::c<DESCRIPTOR>(jPop);
        const int iXs = iXb + cs[0]; //Solid node
        const int iYs = iYb + cs[1];
        const int iZs = iZb + cs[2]; 
        T dist = -1;

        if (blockGeometryStructure.isInside(iXs,iYs,iZs) && !bulkIndicator(iXs, iYs, iZs)) {
         
          Vector<T,3> direction(voxelSize*cs[0],voxelSize*cs[1],voxelSize*cs[2]);
          T cPhysNorm = voxelSize*sqrt(cs[0]*cs[0]+cs[1]*cs[1]+cs[2]*cs[2]);

          if (!geometryIndicator.distance(dist,physC,direction,blockGeometryStructure.getIcGlob())) {
            clout << "Distance not found initially" << std::endl;
            T epsX = voxelSize*cs[0]*this->_epsFraction;
            T epsY = voxelSize*cs[1]*this->_epsFraction;
            T epsZ = voxelSize*cs[2]*this->_epsFraction;

            Vector<T,3> physC2(physC);
            physC2[0] -= epsX;
            physC2[1] -= epsY;
            physC2[2] -= epsZ;
            Vector<T,3> direction2(direction);
            direction2[0] += 2.*epsX;
            direction2[1] += 2.*epsY;
            direction2[2] += 2.*epsZ;
                      
            if ( !geometryIndicator.distance(dist,physC2,direction2,blockGeometryStructure.getIcGlob())) {
              clout << "ERROR: no boundary found at (" << iXb << "," << iYb << "," << iZb <<") ~ ("
                    << physR[0] << "," << physR[1] << "," << physR[2] <<"), "
                    << "in direction " << jPop
                    << std::endl;
            }
            T distNew = (dist - sqrt(epsX*epsX+epsY*epsY+epsZ*epsZ))/cPhysNorm;
            if (distNew < 0.5) {
              dist = 0;
              clout << "DISTANCE WARNING" << std::endl;
            }
            else {
              dist = 0.5 * cPhysNorm;
              clout << "WARNING: distance at (" << iXb << "," << iYb << "," << iZb <<") ~ ("
                    << physR[0] << "," << physR[1] << "," << physR[2] <<"), "
                    << "in direction " << jPop << ": "
                    << distNew
                    << " rounded to "
                    << dist/cPhysNorm
                    << std::endl;
            }
          }
          distances.push_back(dist/cPhysNorm);
          iMissing.push_back(util::opposite<DESCRIPTOR >(jPop));
        } // bulkMaterials if j
      } //jPop loop
      addZeroVelocityGradBoundary(blockGeometryStructure, iXb, iYb, iZb, distances, iMissing);
    }   // bulkMaterials if i
  } // iPop loop
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityGradBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ,
  IndicatorF3D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial3D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  addZeroVelocityGradBoundary(blockGeometryStructure, iX, iY, iZ,
                          geometryIndicator,
                          bulkIndicator);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityGradBoundary(
  BlockIndicatorF3D<T>& boundaryIndicator, BlockIndicatorF3D<T>& bulkIndicator, IndicatorF3D<T>& geometryIndicator)
{
  if ( !boundaryIndicator.isEmpty() ) {
    const Vector<int,3> min = boundaryIndicator.getMin();
    const Vector<int,3> max = boundaryIndicator.getMax();

    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        for (int iZ = min[2]; iZ <= max[2]; ++iZ) {
          if (boundaryIndicator(iX,iY,iZ)) {
            addZeroVelocityGradBoundary(boundaryIndicator.getBlockGeometryStructure(), iX, iY, iZ,
                                    geometryIndicator, bulkIndicator);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist)
{
  const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
  if (blockGeometryStructure.getMaterial(x-c[0], y-c[1], z-c[2]) != 1) {
    addOnePointVelocityBoundary(x, y, z, iPop, dist);
  }
  else {
    addTwoPointVelocityBoundary(x, y, z, iPop, dist);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::
addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, T distances[DESCRIPTOR::q])
{
  T location[DESCRIPTOR::d];
  blockGeometryStructure.getPhysR(location, x,y,z);
  T distancesCopy[DESCRIPTOR::q];
  T spacing = blockGeometryStructure.getDeltaR();
  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    distancesCopy[iPop] = spacing*(1.-distances[iPop]);
    if ( util::nearZero(distances[iPop]+1) ) {
      distancesCopy[iPop] = -1;
    }
  }
  addOffDynamics(x, y, z, location, distancesCopy);

  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    if (!util::nearZero(distances[iPop]+1)) {
      const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
      addVelocityBoundary(blockGeometryStructure, x-c[0], y-c[1], z-c[2], iPop, distances[iPop]);
    }
  }
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ,
  IndicatorF3D<T>& geometryIndicator, BlockIndicatorF3D<T>& bulkIndicator)
{
  T distances[DESCRIPTOR::q];
  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    distances[iPop] = -1;
  }

  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
    const int iXn = iX + c[0];
    const int iYn = iY + c[1];
    const int iZn = iZ + c[2];
    if (blockGeometryStructure.isInside(iXn,iYn,iZn) && bulkIndicator(iXn,iYn,iZn)) {
      T dist = -1;
      T physR[3];
      blockGeometryStructure.getPhysR(physR,iXn,iYn,iZn);
      T voxelSize=blockGeometryStructure.getDeltaR();

      Vector<T,3> physC(physR);
      Vector<T,3> direction(-voxelSize*c[0],-voxelSize*c[1],-voxelSize*c[2]);
      T cPhysNorm = voxelSize*sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

      if (!geometryIndicator.distance(dist,physC,direction,blockGeometryStructure.getIcGlob() ) ) {
        T epsX = voxelSize*c[0]*this->_epsFraction;
        T epsY = voxelSize*c[1]*this->_epsFraction;
        T epsZ = voxelSize*c[2]*this->_epsFraction;

        Vector<T,3> physC2(physC);
        physC2[0] += epsX;
        physC2[1] += epsY;
        physC2[2] += epsZ;
        Vector<T,3> direction2(direction);
        direction2[0] -= 2.*epsX;
        direction2[1] -= 2.*epsY;
        direction2[2] -= 2.*epsZ;

        if ( !geometryIndicator.distance(dist,physC2,direction2,blockGeometryStructure.getIcGlob())) {
          clout << "ERROR: no boundary found at (" << iXn << "," << iYn << "," << iZn <<") ~ ("
                << physR[0] << "," << physR[1] << "," << physR[2] <<"), "
                << "in direction " << util::opposite<DESCRIPTOR >(iPop)
                << std::endl;
        }
        T distNew = (dist - sqrt(epsX*epsX+epsY*epsY+epsZ*epsZ))/cPhysNorm;
        if (distNew < 0.5) {
          dist = 0;
        }
        else {
          dist = 0.5 * cPhysNorm;
          clout << "WARNING: distance at (" << iXn << "," << iYn << "," << iZn <<") ~ ("
                << physR[0] << "," << physR[1] << "," << physR[2] <<"), "
                << "in direction " << util::opposite<DESCRIPTOR >(iPop) << ": "
                << distNew
                << " rounded to "
                << dist/cPhysNorm
                << std::endl;
        }
      }
      distances[util::opposite<DESCRIPTOR >(iPop)] = dist/cPhysNorm;
    } // bulk indicator if
  } // iPop loop
  addVelocityBoundary(blockGeometryStructure, iX, iY, iZ, distances);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ,
  IndicatorF3D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  BlockIndicatorMaterial3D<T> bulkIndicator(blockGeometryStructure, bulkMaterials);
  addVelocityBoundary(blockGeometryStructure, iX, iY, iZ,
                      geometryIndicator,
                      bulkIndicator);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary(
  BlockIndicatorF3D<T>& boundaryIndicator,
  BlockIndicatorF3D<T>& bulkIndicator,
  IndicatorF3D<T>&      geometryIndicator)
{
  if ( !boundaryIndicator.isEmpty() ) {
    const Vector<int,3> min = boundaryIndicator.getMin();
    const Vector<int,3> max = boundaryIndicator.getMax();

    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        for (int iZ = min[2]; iZ <= max[2]; ++iZ) {
          if (boundaryIndicator(iX,iY,iZ)) {
            addVelocityBoundary(bulkIndicator.getBlockGeometryStructure(), iX, iY, iZ,
                                geometryIndicator, bulkIndicator);
          }
        }
      }
    }
  }
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::
setBoundaryIntersection(int iX, int iY, int iZ, int iPop, T distance)
{
  this->getBlock().getDynamics(iX, iY, iZ)->setBoundaryIntersection(iPop, distance);
  if (_output) {
    clout << "setBoundaryIntersection(" << iX << ", " << iY << ", " << iZ << " )" << std::endl;
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
bool OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::
getBoundaryIntersection(int iX, int iY, int iZ, int iPop, T point[DESCRIPTOR::d])
{
  return this->getBlock().getDynamics(iX, iY, iZ)->getBoundaryIntersection(iPop, point);
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::
defineU(int iX, int iY, int iZ, int iPop, const T u[DESCRIPTOR::d])
{
  this->getBlock().getDynamics(iX, iY, iZ)->defineU(iPop, u);
  if (_output) {
    clout << "defineU(" << iX << ", " << iY << ", " << iZ << " )" << std::endl;
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::defineU(
  BlockIndicatorF3D<T>& indicator, BlockIndicatorF3D<T>& bulkIndicator, AnalyticalF3D<T,T>& u)
{
  if ( !indicator.isEmpty() ) {
    const Vector<int,3> min = indicator.getMin();
    const Vector<int,3> max = indicator.getMax();

    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        for (int iZ = min[2]; iZ <= max[2]; ++iZ) {
          if (indicator(iX, iY, iZ)) {
            for (int q = 1; q < DESCRIPTOR::q ; ++q) {
              // Get direction
              const int iXn = iX + descriptors::c<DESCRIPTOR>(q,0);
              const int iYn = iY + descriptors::c<DESCRIPTOR>(q,1);
              const int iZn = iZ + descriptors::c<DESCRIPTOR>(q,2);
              if (bulkIndicator.getBlockGeometryStructure().isInside(iXn,iYn,iZn) &&
                  bulkIndicator(iXn,iYn,iZn)) {
                T intersection[3] = { };
                const int opp = util::opposite<DESCRIPTOR>(q);
                if (this->getBoundaryIntersection(iX, iY, iZ, opp, intersection)) {
                  T vel[3]= { };
                  u(vel, intersection);
                  this->defineU(iX, iY, iZ, opp, vel);
                }
              }
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::outputOn()
{
  _output = true;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::outputOff()
{
  _output = false;
}

}

#endif
