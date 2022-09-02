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

#ifndef OFF_BOUNDARY_POST_PROCESSORS_3D_HH
#define OFF_BOUNDARY_POST_PROCESSORS_3D_HH

#include "offBoundaryPostProcessors3D.h"
#include "core/blockLattice3D.h"
#include "core/util.h"
#include "core/cell.h"

namespace olb {

/////////// SM - QuadraticBouzidiPostProcessor3D /////////////////////////////////////

/* Bouzidi Interpolation scheme of second order
 *
 * fluid nodes               wall  solid node
 * --o------<-o->------<-o->-----<-o->--|----x----
 *            xB2        xB        x         xN
 * directions: --> iPop
 *             <-- opp
 *
*/

template<typename T, typename DESCRIPTOR>
ZeroVelocityBouzidiQuadraticPostProcessor3D<T,DESCRIPTOR>::
ZeroVelocityBouzidiQuadraticPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_)
  : x(x_), y(y_), z(z_), iPop(iPop_)
{
#ifndef QUIET
  if (dist_ < 0 || dist_ > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
              << dist_ << std::endl;
#endif
  typedef DESCRIPTOR L;
  const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
  opp = util::opposite<L>(iPop);
  
  xN = x + c[0]; 
  yN = y + c[1];
  zN = z + c[2];

  if (dist_ >= 0.5) {
    xB1 = x - c[0];
    yB1 = y - c[1];
    zB1 = z - c[2];

    xB2 = x - 2 * c[0];
    yB2 = y - 2 * c[1];
    zB2 = z - 2 * c[2];

    a1 = 1. / (dist_ * (2. * dist_ + 1)); 
    a2 = (2. * dist_ - 1) / dist_;
    a3 = (1. - 2. * dist_) / (1. + 2. * dist_);

    iPop2 = opp;
  } else {
    xB1 = x;
    yB1 = y;
    zB1 = z;

    xB2 = x - c[0]; 
    yB2 = y - c[1];
    zB2 = z - c[2]; 

    a1 = dist_ * (2. * dist_ + 1.);
    a2 = (1. + 2. * dist_) * (1. - 2. * dist_);
    a3 = - dist_ * (1. - 2.* dist_);

    iPop2 = iPop;
  }
  
    /*std::cout << "ZeroVelocityQuadratic (" << x << "," << y << "," << z <<
      "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," << zN <<
      "), opp: " << opp << ", bP: (" << xB << "," << yB << "," << zB <<
      "), dist: " << dist << ", q: " << q << std::endl;*/
  
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBouzidiQuadraticPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBouzidiQuadraticPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  blockLattice.get(x, y, z)[opp] = a1 * blockLattice.get(xN, yN, zN)[iPop] +
                                   a2 * blockLattice.get(xB1, yB1, zB1)[iPop2] +
                                   a3 * blockLattice.get(xB2, yB2, zB2)[iPop2];
                                   
}

/////////// LinearBouzidiPostProcessor3D /////////////////////////////////////

/* Bouzidi Interpolation scheme of first order
 *
 * fluid nodes               wall  solid node
 * --o-------<-o->-----<-o->--|----x----
 *            xB         x        xN
 * directions: --> iPop
 *             <-- opp
 *
*/

template<typename T, typename DESCRIPTOR>
ZeroVelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
ZeroVelocityBouzidiLinearPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_)
  : x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{
#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
  opp = util::opposite<L>(iPop);
  
  xN = x + c[0]; 
  yN = y + c[1];
  zN = z + c[2];

  if (dist >= 0.5) {
    xB = x - c[0];
    yB = y - c[1];
    zB = z - c[2];
    q = 1/(2*dist);
    iPop2 = opp;
  } else {
    xB = x;
    yB = y;
    zB = z;
    q = 2*dist;
    iPop2 = iPop;
  }
  
    /*std::cout << "ZeroVelocityLinear (" << x << "," << y << "," << z <<
      "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," << zN <<
      "), opp: " << opp << ", bP: (" << xB << "," << yB << "," << zB <<
      "), dist: " << dist << ", q: " << q << std::endl;*/
  
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  blockLattice.get(x, y, z)[opp] = q*blockLattice.get(xN, yN, zN)[iPop] +
                                   (1-q)*blockLattice.get(xB, yB, zB)[iPop2];
  //BOUNCEBACK blockLattice.get(x, y, z)[opp] = blockLattice.get(xN, yN, zN)[iPop];
}

template<typename T, typename DESCRIPTOR>
VelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
VelocityBouzidiLinearPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_)
  : x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{
#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];
  zN = z + c[2];

  if (dist >= 0.5) {
    xB = x - c[0];
    yB = y - c[1];
    zB = z - c[2];
    q = 1/(2*dist);
    ufrac = q;
    iPop2 = opp;
  } else {
    xB = x;
    yB = y;
    zB = z;
    q = 2*dist;
    iPop2 = iPop;
    ufrac = 1;
  }
  /*
    std::cout << "VelocityLinear (" << x << "," << y << "," << z <<
      "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," << zN <<
      "), opp: " << opp << ", bP: (" << xB << "," << yB << "," << zB <<
      "), dist: " << dist << ", q: " << q << std::endl;
  */
}

template<typename T, typename DESCRIPTOR>
void VelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void VelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(xN, yN, zN);
  T u = ufrac*dynamics->getVelocityCoefficient(iPop);
  dynamics->defineRho( blockLattice.get(xN, yN, zN), blockLattice.get(x, y, z).computeRho() );
  T j = u;// * blockLattice.get(x, y, z).computeRho();
  blockLattice.get(x, y, z)[opp] = q*blockLattice.get(xN, yN, zN)[iPop] +
                                   (1-q)*blockLattice.get(xB, yB, zB)[iPop2] + j;
}


//////// CornerBouzidiPostProcessor3D ///////////////////

template<typename T, typename DESCRIPTOR>
ZeroVelocityBounceBackPostProcessor3D<T,DESCRIPTOR>::
ZeroVelocityBounceBackPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_)
  : x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{
#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  const Vector<int,3> c = descriptors::c<L>(iPop);
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];
  zN = z + c[2];
  /*
    std::cout << "Corner (" << x << "," << y << "," << z <<
        "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," << zN <<
        "), dist: " << dist << std::endl;
  */
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBounceBackPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityBounceBackPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  blockLattice.get(x, y, z)[opp] = blockLattice.get(xN, yN, zN)[iPop];
}

template<typename T, typename DESCRIPTOR>
VelocityBounceBackPostProcessor3D<T,DESCRIPTOR>::
VelocityBounceBackPostProcessor3D(int x_, int y_, int z_, int iPop_, T dist_)
  : x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{
#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
              << dist << std::endl;
#endif
  typedef DESCRIPTOR L;
  const Vector<int,3> c = descriptors::c<L>(iPop);
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];
  zN = z + c[2];

  /*
    std::cout << "Corner (" << x << "," << y << "," << z <<
        "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," << zN <<
        "), dist: " << dist << std::endl;
  */
}

template<typename T, typename DESCRIPTOR>
void VelocityBounceBackPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void VelocityBounceBackPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(xN, yN, zN);
  T u = dynamics->getVelocityCoefficient(iPop);
  dynamics->defineRho( blockLattice.get(xN, yN, zN), blockLattice.get(x, y, z).computeRho() );
  T j = u;//*blockLattice.get(x, y, z).computeRho();
  blockLattice.get(x, y, z)[opp] = blockLattice.get(xN, yN, zN)[iPop] + j;
}

//Zero velocity Grad post processor
template<typename T, typename DESCRIPTOR>
ZeroVelocityGradPostProcessor3D<T,DESCRIPTOR>::
ZeroVelocityGradPostProcessor3D(int x_, int y_, int z_,
                                std::vector<T> distances_,
                                std::vector<unsigned> iMissing_, BlockGeometryStructure3D<T>& blockGeometryStructure_)
  : x(x_), y(y_), z(z_), blockGeometryStructure(blockGeometryStructure_)
{

  //Option 
  int mode =1;

  //Check bogus distances  
  #ifndef QUIET
  for (unsigned i = 0; i < iMissing_.size(); i++ ) {
    if (distances_[i] < 0 || distances_[i] > 1)
      std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
                << distances_[i] << std::endl;
  }
  #endif
 
  //Generate iNotMissing  
  for (unsigned i = 0; i < DESCRIPTOR::q; i++) {
    bool missing = false;
    for (unsigned j = 0; j < iMissing_.size(); ++j) {
      if (i == iMissing_[j]) {
        missing = true;
        break;
      }
    }
    if (!missing)
      iNotMissing.push_back(i);
  }

  //Sort boundary nodes based on adjacent fluid node properties
  for (unsigned i = 0; i < iMissing_.size(); ++i){
    Vector<int, 3> cf = descriptors::c<DESCRIPTOR>(iMissing_[i]);
    //std::vector<int> xFtemp = getNext(x + cf[0], y + cf[1], z + cf[2]);
    std::vector<int> xFtemp = {x + cf[0], y + cf[1], z + cf[2]};

    bool isClean = true;
    /*for (int j = 0; j < DESCRIPTOR::q; ++j) {
      Vector<int, 3> cf2 = descriptors::c<DESCRIPTOR>(j);
      //std::vector<int> xFtemp2 = getNext(xFtemp[0] + cf2[0], xFtemp[1] + cf2[1], xFtemp[2] + cf2[2]);
      std::vector<int> xFtemp2 = {xFtemp[0] + cf2[0], xFtemp[1] + cf2[1], xFtemp[2] + cf2[2]};
      if (blockGeometryStructure.getMaterial(xFtemp2[0], xFtemp2[1], xFtemp2[2]) != 1) {
        isClean = false;
        break;
      }
    }*/
    
    if (blockGeometryStructure.getMaterial(xFtemp[0], xFtemp[1], xFtemp[2]) != 1) {
      isClean = false;
    }

    if (isClean) {
      iClean.push_back(iMissing_[i]);
      xFclean.push_back(xFtemp[0]);
      yFclean.push_back(xFtemp[1]);
      zFclean.push_back(xFtemp[2]);
      distancesClean.push_back(distances_[i]);
    }
    else if (mode == 1) { //Include dirty nodes into bc
      iDirty.push_back(iMissing_[i]);
      xFdirty.push_back(xFtemp[0]);
      yFdirty.push_back(xFtemp[1]);
      zFdirty.push_back(xFtemp[2]);
      distancesDirty.push_back(distances_[i]);
    }
  }

  //Xb material
  int xBmat = 6;

  //Determine stencil points for pressure tensor gradients
  if ((blockGeometryStructure.getMaterial(x + 1, y, z) == 1) ||
       (blockGeometryStructure.getMaterial(x + 1, y, z) == xBmat)) {
    pStencilX[0] = x + 1;
    pStencilX[1] = x;
  }
  else if ((blockGeometryStructure.getMaterial(x - 1, y, z) == 1) ||
  (blockGeometryStructure.getMaterial(x - 1, y, z) == xBmat)) {
    pStencilX[0] = x;
    pStencilX[1] = x - 1;
  }
  else 
    std::cout << "GradBC Warning: No suitable stencil point in x" << std::endl;
  if ((blockGeometryStructure.getMaterial(x, y + 1, z) == 1) ||
     (blockGeometryStructure.getMaterial(x, y + 1, z) == xBmat)) {
    pStencilY[0] = y + 1;
    pStencilY[1] = y;
  }
  else if ((blockGeometryStructure.getMaterial(x, y - 1, z) == 1) ||
           (blockGeometryStructure.getMaterial(x, y - 1, z) == xBmat)) {
    pStencilY[0] = y;
    pStencilY[1] = y - 1;
  }
  else 
    std::cout << "GradBC Warning: No suitable stencil point in y" << std::endl;
  if ((blockGeometryStructure.getMaterial(x, y, z + 1) == 1) ||
      (blockGeometryStructure.getMaterial(x, y, z + 1) == xBmat)) {
    pStencilZ[0] = z + 1;
    pStencilZ[1] = z;
  }
  else if ((blockGeometryStructure.getMaterial(x, y, z - 1) == 1) ||
           (blockGeometryStructure.getMaterial(x, y, z - 1) == xBmat)) {
    pStencilZ[0] = z;
    pStencilZ[1] = z - 1;
  }
  else 
    std::cout << "GradBC Warning: No suitable stencil point in z" << std::endl;

  //Set nMissing based on mode
  if (mode == 1)
    nMissing = iClean.size();
  else
    nMissing = iClean.size() + iDirty.size();
}

template<typename T, typename DESCRIPTOR>
std::vector<int> ZeroVelocityGradPostProcessor3D<T, DESCRIPTOR>::
getNext(int x, int y, int z)
{
      bool isPeriodic = true;
      if (isPeriodic && (x < 0 || x >= blockGeometryStructure.getNx() ||
          y < 0 || y >= blockGeometryStructure.getNy() ||
          z < 0 || z >= blockGeometryStructure.getNz())) {
        x = (x + blockGeometryStructure.getNx()) % blockGeometryStructure.getNx(); 
        y = (y + blockGeometryStructure.getNy()) % blockGeometryStructure.getNy(); 
        z = (z + blockGeometryStructure.getNz()) % blockGeometryStructure.getNz();
      }

      return std::vector<int> {x, y, z};
}

template<typename T, typename DESCRIPTOR>
std::vector<int> ZeroVelocityGradPostProcessor3D<T, DESCRIPTOR>::
getPosition()
{
      return std::vector<int> {x, y, z};
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityGradPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

template<typename T, typename DESCRIPTOR>
void ZeroVelocityGradPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
using namespace olb::util::tensorIndices3D;

//std::cout << "PROCESS1" << std::endl;

//Compute target density and velocity
T rhoTarget = 1.; //Initialised as sum of weights
T uTarget[3] = {0.,0.,0.};
T invDist;
T ui[3];
T invCs2 = descriptors::invCs2<T, DESCRIPTOR>();
T cs2 = 1. / invCs2;

//std::cout << cs2 << std::endl;

//Target density
for (unsigned i = 0; i < iClean.size(); ++i) {
  rhoTarget += blockLattice.get(x - descriptors::c<DESCRIPTOR>(iClean[i],0),
    y - descriptors::c<DESCRIPTOR>(iClean[i],1), 
    z - descriptors::c<DESCRIPTOR>(iClean[i],2))
    [util::opposite<DESCRIPTOR>(iClean[i])];
    //std::cout << "CLEAN " << iClean[i] << std::endl;
}

for (unsigned i = 0; i < iDirty.size(); ++i) {
  rhoTarget += blockLattice.get(x - descriptors::c<DESCRIPTOR>(iDirty[i],0),
    y - descriptors::c<DESCRIPTOR>(iDirty[i],1), 
    z - descriptors::c<DESCRIPTOR>(iDirty[i],2))
    [util::opposite<DESCRIPTOR>(iDirty[i])];
    //std::cout << "DIRTY " << iDirty[i] << std::endl;
}

//std::cout << "PROCESS2" << std::endl;
for (unsigned i = 0; i < iNotMissing.size(); ++i) { 
  //std::cout << "NOT " << iNotMissing[i] << std::endl;
  rhoTarget += blockLattice.get(x, y, z)[iNotMissing[i]];
}

//Target velocity - only option 1 works for now...
for (unsigned i = 0; i < iClean.size(); ++i) {
    //std::cout << xFclean[i] << " " << yFclean[i] << " " << zFclean[i] << std::endl;
    invDist = distancesClean[i] / (distancesClean[i] + 1);
    blockLattice.get(xFclean[i], yFclean[i], zFclean[i]).computeU(ui);
    //std::cout << ui[0] << " " << ui[1] << " " << ui[2] << std::endl;
    uTarget[0] += invDist * ui[0];
    uTarget[1] += invDist * ui[1];
    uTarget[2] += invDist * ui[2];
}

//std::cout << iDirty.size() << std::endl;
//for (unsigned i = 0; i < iDirty.size(); ++i) {
//    invDist = distancesDirty[i] / (distancesDirty[i] + 1);
//    //std::cout << invDist << std::endl;
//
//    uTarget[0] += invDist * blockLattice.get(xFdirty[i], yFdirty[i], zFdirty[i]).template getField<descriptors::VELOCITY>()[0];
//    uTarget[1] += invDist * blockLattice.get(xFdirty[i], yFdirty[i], zFdirty[i]).template getField<descriptors::VELOCITY>()[1];
//    uTarget[2] += invDist * blockLattice.get(xFdirty[i], yFdirty[i], zFdirty[i]).template getField<descriptors::VELOCITY>()[2];
//}

//std::cout << "PROCESS5" << std::endl;
T invNmissing = 1. / nMissing;
uTarget[0] = uTarget[0] * invNmissing;
uTarget[1] = uTarget[1] * invNmissing;
uTarget[2] = uTarget[2] * invNmissing;

//std::cout << uTarget[0] << " " << uTarget[1] << " " << uTarget[2] << std::endl;
//std::cout << "PROCESS6" << std::endl;
//Compute pressure tensor components - use previous time step velocities 
T pi[util::TensorVal<DESCRIPTOR>::n];
T cs2Beta = cs2 / blockLattice.getDynamics(x, y, z)->getOmega();   

T dx_ux = blockLattice.get(pStencilX[0], y, z).template getField<descriptors::VELOCITY>()[0] - 
          blockLattice.get(pStencilX[1], y, z).template getField<descriptors::VELOCITY>()[0]; 
T dy_ux = blockLattice.get(x, pStencilY[0], z).template getField<descriptors::VELOCITY>()[0] - 
          blockLattice.get(x, pStencilY[1], z).template getField<descriptors::VELOCITY>()[0]; 
T dz_ux = blockLattice.get(x, y, pStencilZ[0]).template getField<descriptors::VELOCITY>()[0] - 
          blockLattice.get(x, y, pStencilZ[1]).template getField<descriptors::VELOCITY>()[0]; 
T dx_uy = blockLattice.get(pStencilX[0], y, z).template getField<descriptors::VELOCITY>()[1] -
          blockLattice.get(pStencilX[1], y, z).template getField<descriptors::VELOCITY>()[1]; 
T dy_uy = blockLattice.get(x, pStencilY[0], z).template getField<descriptors::VELOCITY>()[1] -
          blockLattice.get(x, pStencilY[1], z).template getField<descriptors::VELOCITY>()[1]; 
T dz_uy = blockLattice.get(x, y, pStencilZ[0]).template getField<descriptors::VELOCITY>()[1] -
          blockLattice.get(x, y, pStencilZ[1]).template getField<descriptors::VELOCITY>()[1]; 
T dx_uz = blockLattice.get(pStencilX[0], y, z).template getField<descriptors::VELOCITY>()[2] -
          blockLattice.get(pStencilX[1], y, z).template getField<descriptors::VELOCITY>()[2];  
T dy_uz = blockLattice.get(x, pStencilY[0], z).template getField<descriptors::VELOCITY>()[2] -
          blockLattice.get(x, pStencilY[1], z).template getField<descriptors::VELOCITY>()[2]; 
T dz_uz = blockLattice.get(x, y, pStencilZ[0]).template getField<descriptors::VELOCITY>()[2] -
          blockLattice.get(x, y, pStencilZ[1]).template getField<descriptors::VELOCITY>()[2]; 

pi[xx] = uTarget[0] * uTarget[0] - cs2Beta * 2. * dx_ux; 
pi[yy] = uTarget[1] * uTarget[1] - cs2Beta * 2. * dy_uy; 
pi[zz] = uTarget[2] * uTarget[2] - cs2Beta * 2. * dz_uz; 
pi[xy] = uTarget[0] * uTarget[1] - cs2Beta * (dx_uy + dy_ux);
pi[xz] = uTarget[0] * uTarget[2] - cs2Beta * (dx_uz + dz_ux);
pi[yz] = uTarget[1] * uTarget[2] - cs2Beta * (dy_uz + dz_uy);


//Replace all missing
  for (unsigned i = 0; i < iClean.size(); ++i) {
    Vector<int, 3> ci = descriptors::c<DESCRIPTOR>(iClean[i]);

  //  std::cout << "grad " << pi[xx] << " " << pi[yy] << " " << pi[zz] << " " << pi[xy] << " " << pi[xz] << " " << pi[yz] << " "
  //<< uTarget[0] << " " << uTarget[1] << " " << uTarget[2] << " " << ci[0] << " " << ci[1] << " " << ci[2] << " " <<
  //rhoTarget << " " << cs2 << " " << invCs2 <<  " " << descriptors::t<T,DESCRIPTOR>(iMissing[i]) << std::endl;

    blockLattice.get(x, y, z)[iClean[i]] = descriptors::t<T,DESCRIPTOR>(iClean[i]) * ( rhoTarget * 
      (1. + invCs2 * (uTarget[0] * ci[0] + uTarget[1] * ci[1] + uTarget[2] * ci[2]) +
      0.5 * invCs2 * invCs2 * 
      (pi[xx] * (ci[0] * ci[0] - cs2) +
       pi[yy] * (ci[1] * ci[1] - cs2) +
       pi[zz] * (ci[2] * ci[2] - cs2) + 2. * (
       pi[xy] * ci[0] * ci[1] +
       pi[xz] * ci[0] * ci[2] + 
       pi[yz] * ci[1] * ci[2]))) - 1.); //(fi - ti) needs to be stored
  }

for (unsigned i = 0; i < iDirty.size(); ++i) {
    Vector<int, 3> ci = descriptors::c<DESCRIPTOR>(iDirty[i]);

  //  std::cout << "grad " << pi[xx] << " " << pi[yy] << " " << pi[zz] << " " << pi[xy] << " " << pi[xz] << " " << pi[yz] << " "
  //<< uTarget[0] << " " << uTarget[1] << " " << uTarget[2] << " " << ci[0] << " " << ci[1] << " " << ci[2] << " " <<
  //rhoTarget << " " << cs2 << " " << invCs2 <<  " " << descriptors::t<T,DESCRIPTOR>(iMissing[i]) << std::endl;

    blockLattice.get(x, y, z)[iDirty[i]] = descriptors::t<T,DESCRIPTOR>(iDirty[i]) * ( rhoTarget * 
      (1. + invCs2 * (uTarget[0] * ci[0] + uTarget[1] * ci[1] + uTarget[2] * ci[2]) +
      0.5 * invCs2 * invCs2 * 
      (pi[xx] * (ci[0] * ci[0] - cs2) +
       pi[yy] * (ci[1] * ci[1] - cs2) +
       pi[zz] * (ci[2] * ci[2] - cs2) + 2. * (
       pi[xy] * ci[0] * ci[1] +
       pi[xz] * ci[0] * ci[2] + 
       pi[yz] * ci[1] * ci[2]))) - 1.); //(fi - ti) needs to be stored
  }

}

//////// SM - QuadraticBouzidiBoundaryPostProcessorGenerator ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ZeroVelocityBouzidiQuadraticPostProcessorGenerator3D<T,DESCRIPTOR>::
ZeroVelocityBouzidiQuadraticPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
ZeroVelocityBouzidiQuadraticPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new ZeroVelocityBouzidiQuadraticPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ZeroVelocityBouzidiQuadraticPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ZeroVelocityBouzidiQuadraticPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this->z, this->iPop, this->dist);
}

////////  LinearBouzidiBoundaryPostProcessorGenerator ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ZeroVelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::
ZeroVelocityBouzidiLinearPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
ZeroVelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new ZeroVelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ZeroVelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ZeroVelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
VelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::
VelocityBouzidiLinearPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
VelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new VelocityBouzidiLinearPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
VelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new VelocityBouzidiLinearPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this->z, this->iPop, this->dist);
}

/////////// CornerBouzidiBoundaryPostProcessorGenerator /////////////////////////////////////

template<typename T, typename DESCRIPTOR>
ZeroVelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::
ZeroVelocityBounceBackPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
ZeroVelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new ZeroVelocityBounceBackPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ZeroVelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ZeroVelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
VelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::
VelocityBounceBackPostProcessorGenerator3D(int x_, int y_, int z_, int iPop_, T dist_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), iPop(iPop_), dist(dist_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
VelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new VelocityBounceBackPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->iPop, this->dist);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
VelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new VelocityBounceBackPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this->z, this->iPop, this->dist);
}

//Grad post processor generator - TODO
template<typename T, typename DESCRIPTOR>
ZeroVelocityGradPostProcessorGenerator3D<T,DESCRIPTOR>::
ZeroVelocityGradPostProcessorGenerator3D(int x_, int y_, int z_, std::vector<T> distances_,
  std::vector<unsigned> iMissing_, BlockGeometryStructure3D<T>& blockGeometryStructure_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), distances(distances_), iMissing(iMissing_), blockGeometryStructure(blockGeometryStructure_) 
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
ZeroVelocityGradPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new ZeroVelocityGradPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->distances, this->iMissing, this->blockGeometryStructure);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ZeroVelocityGradPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ZeroVelocityGradPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this-> z, this-> distances, this-> iMissing, this->blockGeometryStructure);
}

}  // namespace olb

#endif
