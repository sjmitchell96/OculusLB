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
                                std::vector<int> iMissing_)
  : x(x_), y(y_), z(z_), distances(distances_), iMissing(iMissing_)
{

  /*
#ifndef QUIET
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << z << "): "
              << dist << std::endl;
#endif

  //Take only distances as input

  //Generate iMissing and iNotMissing
  


  //Cycle through distances and get next fluid nodes (remember periodic) - seperate method

  //Check if nodes are clean or dirty - seperate method

  //If clean, add to clean vector

  //If dirty, snf if option2, add to dirty vector
  //If option 1, just ignore cell (don't add to uF vector)

  //



  //Convert link/bulk arrays to missing/notMissing
  for (unsigned i = 0; i < iLinks_.size(); ++i) 
    iMissing.push_back(util::opposite<DESCRIPTOR>(iLinks_[i]));
  for (unsigned i = 0; i < iBulk_.size(); ++i)
    iNotMissing.push_back(util::opposite<DESCRIPTOR>(iBulk_[i])); 
  opp = util::opposite<DESCRIPTOR>(iPop);

  */
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

//Process dirty nodes vector - seperate object - update f_missing for these ones so that uF can be calculated here

//

/*
//Compute target density and velocity
T rhoTarget = 0.;
T uTarget[3] = {0.,0.,0.};
T invDist;
//T ui[3];
T invCs2 = descriptors::invCs2<T, DESCRIPTOR>();
T cs2 = 1. / invCs2;

for (unsigned i = 0; i < iMissing.size(); ++i) {
  rhoTarget += blockLattice.get(x - descriptors::c<DESCRIPTOR>(iMissing[i],0),
    y - descriptors::c<DESCRIPTOR>(iMissing[i],1), 
    z - descriptors::c<DESCRIPTOR>(iMissing[i],2))
    [util::opposite<DESCRIPTOR>(iMissing[i])]; 

    invDist = distances[i] / (distances[i] + 1);
    uTarget[0] += invDist * blockLattice.get(xF[i], yF[i], zF[i]).template 
      getField<descriptors::VELOCITY>()[0];
    uTarget[1] += invDist * blockLattice.get(xF[i], yF[i], zF[i]).template
     getField<descriptors::VELOCITY>()[1];
    uTarget[2] += invDist * blockLattice.get(xF[i], yF[i], zF[i]).template
     getField<descriptors::VELOCITY>()[2];
    //blockLattice.get(xF[i], yF[i], zF[i]).computeU(ui);
    //uTarget[0] += invDist * ui[0];
    //uTarget[1] += invDist * ui[1];
    //uTarget[2] += invDist * ui[2];
}
for (unsigned i = 0; i < iNotMissing.size(); ++i) 
  rhoTarget += blockLattice.get(x, y, z)[iNotMissing[i]];

T invNmissing = 1. / iMissing.size();
uTarget[0] = uTarget[0] * invNmissing;
uTarget[1] = uTarget[1] * invNmissing;
uTarget[2] = uTarget[2] * invNmissing;

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

Vector<int, 3> cOpp = descriptors::c<DESCRIPTOR>(opp);

blockLattice.get(x, y, z)[opp] = descriptors::t<T,DESCRIPTOR>(opp) * rhoTarget * 
  (1. + invCs2 * (uTarget[0] * cOpp[0] + uTarget[1] * cOpp[1] + uTarget[2] * cOpp[2]) +
  0.5 * invCs2 * invCs2 * 
  (pi[xx] * (cOpp[0] * cOpp[0] - cs2) +
   pi[yy] * (cOpp[1] * cOpp[1] - cs2) +
   pi[zz] * (cOpp[2] * cOpp[2] - cs2) + 2. * (
   pi[xy] * cOpp[0] * cOpp[1] +
   pi[xz] * cOpp[0] * cOpp[2] + 
   pi[yz] * cOpp[1] * cOpp[2])));
*/
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
ZeroVelocityGradPostProcessorGenerator3D(int x_, int y_, int z_, std::vector<T> distances_, std::vector<int> iMissing_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_, x_, y_, y_, z_, z_),
    x(x_), y(y_), z(z_), distances(distances_), iMissing(iMissing_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
ZeroVelocityGradPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new ZeroVelocityGradPostProcessor3D<T,DESCRIPTOR>
         ( this->x, this->y, this->z, this->distances, this->iMissing);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ZeroVelocityGradPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ZeroVelocityGradPostProcessorGenerator3D<T,DESCRIPTOR>
         (this->x, this->y, this-> z, this-> distances, this-> iMissing);
}

}  // namespace olb

#endif
