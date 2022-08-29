/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2013 Mathias J. Krause, Thomas Henn, Tim Dornieden
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

/* cylinder3d.cpp:
 * This example examines a steady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
 * condition set by p = 0.
 * Inspired by "Benchmark Computations of Laminar Flow Around
 * a Cylinder" by M.Schäfer and S.Turek. For high resolution, low
 * latticeU, and enough time to converge, the results for pressure drop, drag
 * and lift lie within the estimated intervals for the exact results.
 * An unsteady flow with Karman vortex street can be created by changing the
 * Reynolds number to Re=100.
 * It also shows the usage of the STL-reader and explains how
 * to set boundary conditions automatically.
 */


#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code
#endif
//Include 2D code for 2D vtk writer!
#include "olb2D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb2D.hh"   // include full template code
#endif

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <type_traits>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;

//Collision choice
//#define WALE
//#define Smagorinsky
#define KBC
#define sponge 

//Boundary condition choice
#define Bouzidi 
//#define Grad

#ifdef WALE
#define DESCRIPTOR WALED3Q19Descriptor
#elif defined (Smagorinsky)
#define DESCRIPTOR D3Q19<>
#elif defined (KBC)
#if defined (sponge)
#define DESCRIPTOR D3Q27descriptorKBCSponge
#else
#define DESCRIPTOR D3Q27descriptorKBC
#endif
#endif

#ifdef Grad
#define DESCRIPTOR D3Q27descriptorKBCGrad
#endif

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( Grid3D<T,DESCRIPTOR>& grid,
		      Vector<T,3> const& origin,
                      Vector<T,3> const& extend,
                      IndicatorCylinder3D<T>& indicatorCylinder) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  auto& converter  = grid.getConverter();
  auto& sGeometry  = grid.getSuperGeometry();
  const T deltaX   = converter.getPhysDeltaX();
  const T diameter    = 2. * indicatorCylinder.getRadius();

  const Vector<T,3> cylinderOrigin = indicatorCylinder.getCenter1();
  const Vector<T,3> cylinderExtend = indicatorCylinder.getCenter2() - cylinderOrigin;
  const Vector<T,3> gridOrigin  = 
	  grid.getSuperGeometry().getStatistics().getMinPhysR(0);
  const Vector<T,3> gridExtend = 
          grid.getSuperGeometry().getStatistics().getPhysExtend(0);

  sGeometry.rename(0,1);
  sGeometry.rename(1, 5, indicatorCylinder);

  //Material number for section of blade to read pressures
  const Vector<T,3> pressureSection1Origin {cylinderOrigin[0] - diameter,
	                                          cylinderOrigin[1] - diameter,
                                            cylinderOrigin[2] + 0.5 * 
                                            cylinderExtend[2] - 0.5 * deltaX};
  const Vector<T,3> pressureSection1Extend {2. * diameter, 2. * diameter, deltaX};
  IndicatorCuboid3D<T> pressureSection1(pressureSection1Extend,
                                        pressureSection1Origin);
  sGeometry.rename(5, 7, pressureSection1);

  //Front face
  {
    const Vector<T,3> indiOrigin {origin[0] - deltaX / 2.,
                                  origin[1] - deltaX / 2.,
                                  origin[2] - 50 * deltaX / 2.};
    const Vector<T,3> indiExtend {deltaX,
                                  extend[1] +deltaX,
                                  extend[2] + 100 * deltaX};
    IndicatorCuboid3D<T> ff(indiExtend, indiOrigin);
    sGeometry.rename(1, 3, ff);
  }

  //Upper face
  {
    const Vector<T,3> indiOrigin {origin[0] + deltaX / 2,
                                  origin[1] + extend[1] - deltaX / 2.,
                                  origin[2] - 50 * deltaX / 2.};
    const Vector<T,3> indiExtend {extend[0] + deltaX,
                                  deltaX,
                                  extend[2] + 100 * deltaX};
    IndicatorCuboid3D<T> uf(indiExtend, indiOrigin);
    sGeometry.rename(1, 3, uf);
  }

  //Lower face
  {
    const Vector<T,3> indiOrigin {origin[0] + deltaX / 2,
                                  origin[1] - deltaX / 2.,
                                  origin[2] - 50 * deltaX / 2.};
    const Vector<T,3> indiExtend {extend[0] + deltaX,
                                  deltaX,
                                  extend[2] + 100 * deltaX};
    IndicatorCuboid3D<T> lf(indiExtend, indiOrigin);
    sGeometry.rename(1, 3, lf);
  }

  //Rear face
  {
    const Vector<T,3> indiOrigin {origin[0] + extend[0] - deltaX / 2.,
                                  origin[1] + deltaX / 2.,
                                  origin[2] - 50 * deltaX / 2.};
    const Vector<T,3> indiExtend {deltaX,
                                  extend[1] - deltaX,
                                  extend[2] + 100 * deltaX};
    IndicatorCuboid3D<T> rf(indiExtend, indiOrigin);
    sGeometry.rename(1, 4, rf);
  }

  // Removes all not needed boundary voxels outside the surface
  sGeometry.clean();

  #ifdef Grad
  IndicatorLayer3D<T> cylinderLayer(indicatorCylinder, deltaX);
  sGeometry.rename(1,6,cylinderLayer);
  #endif

  sGeometry.checkForErrors();
  sGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void setupRefinement(Grid3D<T,DESCRIPTOR>& coarseGrid,
		     Vector<T,3> const& domainOrigin,
		     Vector<T,3> const& domainExtend,
         IndicatorCylinder3D<T>& indicatorCylinder,
		     const int n) {

  T radius = indicatorCylinder.getRadius();
  T diameter = 2. * radius;
  Vector<T,3> cylinderOrigin = indicatorCylinder.getCenter1();
  OstreamManager clout(std::cout, "setupRefinement");
  clout << "Setup Refinement ..." << std::endl;

  //Origin of sphere bounding box
  Vector<T,3> cylinderBoxOrigin = {cylinderOrigin[0] - diameter / 2.,
	                                 cylinderOrigin[1] - diameter / 2., 
				                           cylinderOrigin[2]};

  //Heights around wing box for each refinement level
  //x,y heights in negative direction //Innermost
  const Vector<T,2> hn4 = {0.1 * diameter, 0.1 * diameter}; 
  const Vector<T,2> hp4 = {0.5 * diameter, 0.1 * diameter}; // '' positive

  const Vector<T,2> hn3 = {0.3 * diameter, 0.3 * diameter};
  const Vector<T,2> hp3 = {6.0 * diameter, 0.3 * diameter};

  const Vector<T,2> hn2 = {0.7 * diameter, 0.7 * diameter};
  const Vector<T,2> hp2 = {12.0 * diameter, 0.7 * diameter};

  const Vector<T,2> hn1 = {1.5 * diameter, 1.5 * diameter}; //Outermost
  const Vector<T,2> hp1 = {16.0 * diameter, 1.5 * diameter};

  if(n >= 1) {
    // Refinement around the wing box - level 1
    T coarseDeltaX = coarseGrid.getConverter().getPhysDeltaX();
    Vector<T,3> fineOrigin = 
      {cylinderBoxOrigin[0] - hn1[0], cylinderBoxOrigin[1] - hn1[1], domainOrigin[2]};
    Vector<T,3> fineExtend 
      = {diameter+hp1[0] + hn1[0], diameter + hp1[1] + hn1[1], domainExtend[2]};

    auto& fineGrid = coarseGrid.refine(fineOrigin, fineExtend, false, false,
		                       true, false);
    prepareGeometry(fineGrid, domainOrigin, domainExtend, indicatorCylinder);

    Vector<T,3> origin = fineGrid.getOrigin() +
      Vector<T,3> {0., 0., 0.5 * coarseDeltaX};
    Vector<T,3> extend = fineGrid.getExtend() -
      Vector<T,3> {0., 0., 0.5 * coarseDeltaX};

    Vector<T,3> extendXZ = {extend[0], 0., extend[2]};
    Vector<T,3> extendYZ = {0., extend[1], extend[2]};
    coarseGrid.addFineCoupling(fineGrid, origin, extendXZ);
    coarseGrid.addFineCoupling(fineGrid, origin, extendYZ);

    Vector<T,3> extendX	= {extend[0], 0., 0.};
    Vector<T,3> extendY	= {0., extend[1], 0.};
    coarseGrid.addFineCoupling(fineGrid, origin + extendX, extendYZ);
    coarseGrid.addFineCoupling(fineGrid, origin + extendY, extendXZ);

    Vector<T,3> innerOrigin =
      origin + Vector<T,3> {coarseDeltaX, coarseDeltaX, 0.};
    Vector<T,3> innerExtendXZ = {extend[0] - 2. * coarseDeltaX, 0., extend[2]};
    Vector<T,3> innerExtendYZ = {0., extend[1] - 2. * coarseDeltaX, extend[2]};
    coarseGrid.addCoarseCoupling(fineGrid, innerOrigin, innerExtendXZ);
    coarseGrid.addCoarseCoupling(fineGrid, innerOrigin, innerExtendYZ);

    Vector<T,3> innerExtendX = {extend[0] - 2. * coarseDeltaX, 0., 0.};
    Vector<T,3> innerExtendY = {0, extend[1] - 2. * coarseDeltaX, 0.};
    coarseGrid.addCoarseCoupling(fineGrid, innerOrigin + innerExtendX,
		                 innerExtendYZ);
    coarseGrid.addCoarseCoupling(fineGrid, innerOrigin + innerExtendY,
				 innerExtendXZ);

    Vector<T,3> refinedOrigin = origin + Vector<T,3> {2. * coarseDeltaX,
	                                              2. * coarseDeltaX,
	                                              - 2. * coarseDeltaX};
    Vector<T,3> refinedExtend = extend - Vector<T,3> {4. * coarseDeltaX,
	                                              4. * coarseDeltaX,
	                                              - 4. * coarseDeltaX};
    IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
    coarseGrid.getSuperGeometry().reset(refined);

    if(n >= 2) {
      // Refinement around the wing box - level 2
      const T deltaX0 = fineGrid.getConverter().getPhysDeltaX();
      coarseDeltaX = deltaX0;	
      Vector<T,3> fineOrigin2 = {cylinderBoxOrigin[0] - hn2[0],
	                         cylinderBoxOrigin[1] - hn2[1],
	                         domainOrigin[2] - deltaX0};
      Vector<T,3> fineExtend2 = {diameter + hp2[0] + hn2[0],
	                         diameter+hp2[1] + hn2[1],
	                         domainExtend[2] + deltaX0};
      auto& fineGrid2 = fineGrid.refine(fineOrigin2, fineExtend2, false, false,
		                        true, false);
      prepareGeometry(fineGrid2, domainOrigin, domainExtend, indicatorCylinder);
      origin = fineGrid2.getOrigin() + Vector<T,3> {0., 0., 0.5 * coarseDeltaX};
      extend = fineGrid2.getExtend() - Vector<T,3> {0., 0., 0.5 * coarseDeltaX};
      extendXZ = {extend[0], 0., extend[2]};
      extendYZ = {0., extend[1], extend[2]};
      fineGrid.addFineCoupling(fineGrid2, origin, extendXZ);
      fineGrid.addFineCoupling(fineGrid2, origin, extendYZ);
      extendX = {extend[0], 0., 0.};
      extendY = {0., extend[1], 0.};
      fineGrid.addFineCoupling(fineGrid2, origin + extendX, extendYZ);
      fineGrid.addFineCoupling(fineGrid2, origin + extendY, extendXZ);

      innerOrigin = origin + Vector<T,3> {coarseDeltaX, coarseDeltaX, 0.};
      innerExtendXZ = {extend[0] - 2. * coarseDeltaX, 0, extend[2]};
      innerExtendYZ = {0., extend[1] - 2. * coarseDeltaX, extend[2]};
      fineGrid.addCoarseCoupling(fineGrid2, innerOrigin, innerExtendXZ);
      fineGrid.addCoarseCoupling(fineGrid2, innerOrigin, innerExtendYZ);

      innerExtendX = {extend[0] - 2. * coarseDeltaX, 0., 0.};
      innerExtendY = {0., extend[1] - 2. * coarseDeltaX, 0.};
      fineGrid.addCoarseCoupling(fineGrid2, innerOrigin + innerExtendX,
				 innerExtendYZ);
      fineGrid.addCoarseCoupling(fineGrid2, innerOrigin + innerExtendY,
				 innerExtendXZ);

      refinedOrigin = origin + 
        Vector<T,3> {2. * coarseDeltaX, 1.9 * coarseDeltaX, - 2. * coarseDeltaX};
      refinedExtend = extend - 
        Vector<T,3> {4. * coarseDeltaX, 3.9 * coarseDeltaX, - 4. * coarseDeltaX};
      IndicatorCuboid3D<T> refined2(refinedExtend, refinedOrigin);
      fineGrid.getSuperGeometry().reset(refined2);

      if(n >= 3) {
        // Refinement around the wing box - level 3
	      const T deltaX1 = fineGrid2.getConverter().getPhysDeltaX();
	      coarseDeltaX = deltaX1;
	      Vector<T,3> fineOrigin3 = {cylinderBoxOrigin[0] - hn3[0], 
		                   cylinderBoxOrigin[1] - hn3[1],
				               domainOrigin[2] - deltaX0 - deltaX1};
	      Vector<T,3> fineExtend3 = {diameter + hp3[0] + hn3[0],
		                   diameter + hp3[1] + hn3[1],
				               domainExtend[2] + deltaX0 + deltaX1};

	      auto& fineGrid3 = fineGrid2.refine(fineOrigin3, fineExtend3, false,
		                           false, true, false);
	      prepareGeometry(fineGrid3, domainOrigin, domainExtend, indicatorCylinder);

	      origin = fineGrid3.getOrigin() + Vector<T,3>{0., 0., 0.5 * coarseDeltaX};
	      extend = fineGrid3.getExtend() - Vector<T,3>{0., 0., 0.5 * coarseDeltaX};
        extendXZ = {extend[0], 0., extend[2]};
	      extendYZ = {0., extend[1], extend[2]};
	      fineGrid2.addFineCoupling(fineGrid3, origin, extendXZ);
	      fineGrid2.addFineCoupling(fineGrid3, origin, extendYZ);
                   extendX	= {extend[0], 0., 0.};
	      extendY	= {0., extend[1], 0.};
        fineGrid2.addFineCoupling(fineGrid3, origin + extendX, extendYZ);
	      fineGrid2.addFineCoupling(fineGrid3, origin + extendY, extendXZ);

       	innerOrigin = origin + Vector<T,3> {coarseDeltaX, coarseDeltaX, 0.};
        innerExtendXZ = {extend[0] - 2. * coarseDeltaX, 0., extend[2]};
	      innerExtendYZ = {0., extend[1] - 2. * coarseDeltaX, extend[2]};
	      fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin, innerExtendXZ);
	      fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin, innerExtendYZ);

       	innerExtendX = {extend[0] - 2. * coarseDeltaX, 0., 0.};
	      innerExtendY = {0., extend[1] - 2. * coarseDeltaX, 0.};
	      fineGrid2.addCoarseCoupling(fineGrid3,
		                    innerOrigin + innerExtendX,
				    innerExtendYZ);
      	fineGrid2.addCoarseCoupling(fineGrid3,
		            	    innerOrigin + innerExtendY,
          				    innerExtendXZ);

       	refinedOrigin = origin + Vector<T,3> {2. * coarseDeltaX, 
		                              2. * coarseDeltaX,
					      - 2. * coarseDeltaX};
	      refinedExtend = extend - Vector<T,3> {4. * coarseDeltaX,
	                                      4. * coarseDeltaX,
					      - 4. * coarseDeltaX};
	      IndicatorCuboid3D<T> refined3(refinedExtend, refinedOrigin);
	      fineGrid2.getSuperGeometry().reset(refined3);

	      if(n >= 4) {
	      // Refinement around the wing box - level 4 (current innermost)
          const T deltaX2 = fineGrid3.getConverter().getPhysDeltaX();
	          coarseDeltaX = deltaX2;
			
	        Vector<T,3> fineOrigin4 = {cylinderBoxOrigin[0] - hn4[0],
		                     cylinderBoxOrigin[1] - hn4[1],
				                 domainOrigin[2]-deltaX0-deltaX1-deltaX2};
	        Vector<T,3> fineExtend4 = 
	        {diameter + hp4[0] + hn4[0], diameter + hp4[1] + hn4[1],
               domainExtend[2] + deltaX0 + deltaX1 + deltaX2};

      	  auto& fineGrid4 = fineGrid3.refine(fineOrigin4, fineExtend4,
				             false, false, true, false);
	        prepareGeometry(fineGrid4, domainOrigin, domainExtend,
			        indicatorCylinder);
	        origin = fineGrid4.getOrigin() + 
               Vector<T,3>{0., 0., 0.5 * coarseDeltaX};
      	  extend = fineGrid4.getExtend() - 
            Vector<T,3>{0., 0., 0.5 * coarseDeltaX};
	        extendXZ = {extend[0], 0., extend[2]};
	        extendYZ = {0., extend[1], extend[2]};
            fineGrid3.addFineCoupling(fineGrid4, origin, extendXZ);
	        fineGrid3.addFineCoupling(fineGrid4, origin, extendYZ);
          extendX = {extend[0], 0., 0.};
	        extendY = {0., extend[1], 0.};
            fineGrid3.addFineCoupling(fineGrid4, origin + extendX, extendYZ);
	        fineGrid3.addFineCoupling(fineGrid4, origin + extendY, extendXZ);

      	  innerOrigin = origin + Vector<T,3> {coarseDeltaX, coarseDeltaX, 0.};
          innerExtendXZ = {extend[0] - 2. * coarseDeltaX , 0., extend[2]};
	        innerExtendYZ = {0., extend[1] - 2. * coarseDeltaX, extend[2]};
	        fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin, innerExtendXZ);
          fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin, innerExtendYZ);

      	  innerExtendX = {extend[0] - 2. * coarseDeltaX, 0., 0.};
	        innerExtendY = {0., extend[1] - 2. * coarseDeltaX, 0.};
	        fineGrid3.addCoarseCoupling(fineGrid4,
			              innerOrigin + innerExtendX,
				      innerExtendYZ);
	        fineGrid3.addCoarseCoupling(fineGrid4,
			              innerOrigin + innerExtendY,
				      innerExtendXZ);

      	  refinedOrigin = origin + Vector<T,3> {2. * coarseDeltaX,
		                                2. * coarseDeltaX,
					       	- 2. * coarseDeltaX};
	        refinedExtend = extend - Vector<T,3> {4. * coarseDeltaX,
		                                4. * coarseDeltaX,
			                      	- 4. * coarseDeltaX};
	        IndicatorCuboid3D<T> refined4(refinedExtend, refinedOrigin);
	        fineGrid3.getSuperGeometry().reset(refined4);
  	    }
      }
    }
  }
  clout << "Setup Refinement ... OK" << std::endl;
}

// Create lattice structures
void prepareLattice(Grid3D<T,DESCRIPTOR>& grid,
		    IndicatorCylinder3D<T>& indicatorCylinder,
		    const bool& bouzidiOn) {
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare lattice ..." << std::endl;

  auto& converter = grid.getConverter();
  auto& sGeometry = grid.getSuperGeometry();
  auto& sLattice  = grid.getSuperLattice();
  const T omega	= converter.getLatticeRelaxationFrequency();

  // Initialize dynamics
  #if defined(WALE)
   Dynamics<T,DESCRIPTOR>& bulkDynamics = 
    grid.addDynamics(std::unique_ptr<Dynamics<T,DESCRIPTOR>>(
      new WALEBGKdynamics<T,DESCRIPTOR>(
        omega, instances::getBulkMomenta<T,DESCRIPTOR>(),0.5))); 
  #elif defined(Smagorinsky)
    Dynamics<T,DESCRIPTOR>& bulkDynamics = 
      grid.addDynamics(std::unique_ptr<Dynamics<T,DESCRIPTOR>>(
        new SmagorinskyBGKdynamics<T,DESCRIPTOR>(
          omega, instances::getBulkMomenta<T,DESCRIPTOR>(), 0.1)));
  #elif defined(KBC)
    #if defined(Grad)
      Dynamics<T,DESCRIPTOR>& bulkDynamics = 
        grid.addDynamics(std::unique_ptr<Dynamics<T,DESCRIPTOR>>(
          new KBCGradDynamics<T,DESCRIPTOR>(
            omega, instances::getKBCBulkMomenta<T,DESCRIPTOR>())));
    #elif defined(sponge)
      Dynamics<T,DESCRIPTOR>& bulkDynamics = 
        grid.addDynamics(std::unique_ptr<Dynamics<T,DESCRIPTOR>>(
          new KBCSpongeDynamics<T,DESCRIPTOR>(
          //new KBCdynamics<T,DESCRIPTOR>(
            omega, instances::getKBCBulkMomenta<T,DESCRIPTOR>())));
    #else
     Dynamics<T,DESCRIPTOR>& bulkDynamics = 
        grid.addDynamics(std::unique_ptr<Dynamics<T,DESCRIPTOR>>(
          new KBCdynamics<T,DESCRIPTOR>(
            omega, instances::getKBCBulkMomenta<T,DESCRIPTOR>())));
    #endif
  #endif

  // Initialize boundary condition types
  //Interp
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc = 
    grid.getOnLatticeBoundaryCondition();
  createInterpBoundaryCondition3D<T,DESCRIPTOR>(bc);

  //Local
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& onbc =
    grid.getOnLatticeBoundaryCondition();
  createLocalBoundaryCondition3D<T,DESCRIPTOR>(onbc);

  sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc =
      grid.getOffLatticeBoundaryCondition();

  #if defined(Bouzidi)
    createBouzidiBoundaryCondition3D<T,DESCRIPTOR>(offBc);
  #elif defined(Grad)
    createGradBoundaryCondition3D<T,DESCRIPTOR>(offBc);
  #endif

  // Define dynamics
  sLattice.defineDynamics(sGeometry, 0,
		       	  &instances::getNoDynamics<T,DESCRIPTOR>());
//  sLattice.defineDynamics(sGeometry, 2,
//		       	  &instances::getNoDynamics<T,DESCRIPTOR>());

  auto bulkIndicator = sGeometry.getMaterialIndicator({1, 2, 3, 4, 6});
  sLattice.defineDynamics(bulkIndicator, &bulkDynamics);

  // Define boundary conditions
  bc.addVelocityBoundary(sGeometry, 3, omega);
  bc.addPressureBoundary(sGeometry, 4, omega); 

  #if defined(Bouzidi)
    // material=5, 7 --> no dynamics + bouzidi zero velocity
    sLattice.defineDynamics( sGeometry,5,&instances::getNoDynamics<T,DESCRIPTOR>() );
    sLattice.defineDynamics( sGeometry,7,&instances::getNoDynamics<T,DESCRIPTOR>() );
    offBc.addSecondOrderZeroVelocityBoundary( sGeometry,5,indicatorCylinder );
    offBc.addSecondOrderZeroVelocityBoundary( sGeometry,7,indicatorCylinder );
  #elif defined(Grad)
    sLattice.defineDynamics( sGeometry,5,&instances::getNoDynamics<T,DESCRIPTOR>() );
    sLattice.defineDynamics( sGeometry,7,&instances::getNoDynamics<T,DESCRIPTOR>() );
    offBc.addZeroVelocityGradBoundary( sGeometry,5,indicatorCylinder,std::vector<int>{1,6} );
    offBc.addZeroVelocityGradBoundary( sGeometry,7,indicatorCylinder,std::vector<int>{1,6} );
  #else
    //material=5,7 --> fullway bounceBack dynamics
    sLattice.defineDynamics( sGeometry, 5, &instances::getBounceBack<T, DESCRIPTOR>() );
    sLattice.defineDynamics( sGeometry, 7, &instances::getBounceBack<T, DESCRIPTOR>() );
  #endif


  //Define and initialise viscosity sponge zones
  //Sponge indicator
  const T physChord = 1.00;
  const T deltaX = converter.getPhysDeltaX();
  const Vector<T,3> spongeOrigin = {18. * physChord - deltaX /2000, - deltaX / 2,
    - 4. * deltaX};
  const Vector<T,3> spongeExtend = {2. * physChord + deltaX / 1000,
    10. * physChord + deltaX,
    3. * physChord + 8. * deltaX};
  IndicatorCuboid3D<T> spongeRegion(spongeExtend, spongeOrigin);
  //Orientation
  const Vector<T,3> spongeOrientation = {1., 0., 0.};
  //Min and max tau limits
  const T tauSpongeBase = 1. / omega;
  const T tauSpongeMax = 1.;
  std::vector<int> spongeMaterials = {1,2,3,4,6};

  sViscositySponge3D<T,DESCRIPTOR>& outletSponge =  grid.getViscositySponge();
  createViscositySponge3D(outletSponge);

  outletSponge.addSineSponge(sGeometry, spongeRegion, spongeOrientation,
    tauSpongeBase, tauSpongeMax, spongeMaterials);

/*
  //Sponge indicator 2 - x-z 
  const Vector<T,3> spongeOrigin2 = {0. * physChord - deltaX /2, 12. * physChord - deltaX / 2000,
    - 4 * deltaX};
  const Vector<T,3> spongeExtend2 = {32. * physChord + deltaX,
    4. * physChord + deltaX / 1000,
    8 * physChord + 8 * deltaX};
  IndicatorCuboid3D<T> spongeRegion2(spongeExtend2, spongeOrigin2);
  //Orientation
  const Vector<T,3> spongeOrientation2 = {0., 1., 0.};

  //sViscositySponge3D<T,DESCRIPTOR>& outletSponge2 =  grid.getViscositySponge();
  //createViscositySponge3D(outletSponge2);

  //outletSponge2.addSineSponge(sGeometry, spongeRegion2, spongeOrientation2,
  //  tauSpongeBase, tauSpongeMax, spongeMaterials);

  outletSponge.addSineSponge(sGeometry, spongeRegion2, spongeOrientation2,
    tauSpongeBase, tauSpongeMax, spongeMaterials);
*/
  sLattice.initialiseSponges();

  // Initial conditions - characteristic physical velocity and density for inflow
  AnalyticalConst3D<T,T> rhoF {1.};
  Vector<T,3> velocityV {converter.getCharLatticeVelocity(), 0., 0.};
  AnalyticalConst3D<T,T> uF(velocityV);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Set boundary values for start-scale inlet velocity
void setBoundaryValues(Grid3D<T,DESCRIPTOR>& grid, int iT) {

	OstreamManager clout(std::cout, "setBoundaryValues");

	auto& converter	= grid.getConverter();
	auto& sGeometry	= grid.getSuperGeometry();
	auto& sLattice	= grid.getSuperLattice();

    Vector<T,3> inVel {converter.getCharLatticeVelocity(), 0., 0.};
    //inVel[0] = converter.getCharLatticeVelocity();
    if (iT > 200 && iT < 400) {
      inVel[1] = converter.getCharLatticeVelocity();
    }
    T inRho = 1.0;
    
    AnalyticalConst3D<T,T> inRhoConst(inRho);
    AnalyticalConst3D<T,T> inVelConst(inVel);

    //sLattice.defineU(sGeometry, 2, inVelConst);
    //sLattice.defineU(sGeometry, 3, inVelConst);
    //sLattice.defineRhoU(sGeometry, 3, inRhoConst, inVelConst);
    sLattice.iniEquilibrium(sGeometry, 2, inRhoConst, inVelConst);
    sLattice.iniEquilibrium(sGeometry, 3, inRhoConst, inVelConst);
}

// Output results to vtk files
void getVTK(Grid3D<T,DESCRIPTOR>& grid, const std::string& prefix, int iT,
	    SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>& wssp,
	    SuperLatticeTimeAveragedF3D<T>& sAveragedWSSP) {

  OstreamManager clout( std::cout,"getVTK" );

  auto& converter = grid.getConverter();
  auto& sLattice  = grid.getSuperLattice();
  auto& sGeometry = grid.getSuperGeometry();

  SuperVTMwriter3D<T> vtmWriter(prefix);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure(sLattice, converter);
  SuperLatticeGeometry3D<T,DESCRIPTOR> geometry(sLattice, sGeometry);
  vtmWriter.addFunctor(geometry);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);
  vtmWriter.addFunctor(wssp);
  vtmWriter.addFunctor(sAveragedWSSP);

  if (iT==0) {
     vtmWriter.createMasterFile();
  }

  vtmWriter.write(iT);
}

void getStats( Grid3D<T,DESCRIPTOR>& grid, int iT,
               Timer<T> timer) {
  auto& sLattice  = grid.getSuperLattice();
  auto& converter = grid.getConverter();
  // Writes output on the console
  // Timer console output
  timer.update( iT );
  timer.printStep();
  // Lattice statistics console output
  sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
}

void getCylinderForce(Grid3D<T,DESCRIPTOR>& grid,
		   const std::string& filePath,
		   IndicatorCylinder3D<T>& cylinder) {
  auto& sLattice = grid.getSuperLattice();
  auto& superGeometry = grid.getSuperGeometry();
  auto& converter = grid.getConverter();
  T span = (cylinder.getCenter2() - cylinder.getCenter1())[2] * 0.5; 
  T diameter = cylinder.getRadius() * 2.;

  SuperLatticePhysDragBlade3D<T,DESCRIPTOR> drag(sLattice, superGeometry,
		                                 std::vector<int>{5,7},
						 converter,span*diameter);
  T cylinderForce[3];
  int input1[0];
  drag(cylinderForce,input1);
  //std::cout << "Lift = " << bladeForce[1] << endl;
  //std::cout << "Drag = " << bladeForce[0] << endl;

  ofstream myfile;
  std::string filename {filePath+".csv"};
  myfile.open(filename,fstream::app);
  myfile << "0.0  " << cylinderForce[1] << "  " << cylinderForce[0] << std::endl;
  myfile.close();
}

// Capture the pressure around the cylinder --- middle z
void capturePressure(Grid3D<T,DESCRIPTOR>& grid,
                     IndicatorCylinder3D<T>& indicatorCylinder,
                      int iT) {
	auto& sLattice = grid.getSuperLattice();
	auto& converter = grid.getConverter();
	const T dR = converter.getPhysDeltaX();
	SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
	AnalyticalFfromSuperF3D<T> interpolatePressure(pressure, true);
	const T radiusCylinder = indicatorCylinder.getRadius();
  const T cylinderCenterX = indicatorCylinder.getCenter1()[0];
  const T cylinderCenterY = indicatorCylinder.getCenter1()[1];
  const T lz = indicatorCylinder.getCenter2()[2] - indicatorCylinder.getCenter1()[2];

	ofstream myfile;
	std::string filename {"tmp/pressure_" + to_string(iT) + ".csv"};
	myfile.open(filename, fstream::app);

	for(int i = 0; i < 180; ++i) {
		const T degree = i * 2. * 3.14159 / 180.0;
		const T point[3] {cylinderCenterX
						  + (radiusCylinder + 1.5 * dR) * std::cos(degree),
						  cylinderCenterY
						  + (radiusCylinder + 1.5 * dR) * std::sin(degree),
						  lz/2.};
		T pressureAtPoint {};
		interpolatePressure(&pressureAtPoint, point);

		myfile << degree << " " << pressureAtPoint << std::endl;
	}

   myfile.close();
}

int main( int argc, char* argv[] ) {

  //Cylinder parameters
  const T diameter = 1.00;
  const T span = 3. * diameter;
  const Vector<T,3> cylinderOrigin = {5.0 * diameter + 0.00045,
                                      5.0 * diameter + 0.00045,
                                      -0.5 * span}; 
  const Vector<T,3> cylinderExtend = {0.0 * diameter,
                                      0.0 * diameter,
                                      2. * span}; 


  //Domain and simulation parameters
  const int N = 22; //14        // resolution of the model (coarse cells per chord)
  const int nRefinement = 4;	//Number of refinement levels (current max = 4)
  const T lDomainPhysx = 20.*diameter; //Length of domain in physical units (m)
  const T lDomainPhysy = 10.*diameter;
  const T lDomainPhysz = span; //
  const T maxPhysT = 100; // max. simulation time in s, SI unit
  const T physL = diameter; //Physical reference length (m)

  //Flow conditions
  const T Re = 3900.;       // Reynolds number
  const T Mach = 0.12;
  const T uC = Mach * 1./std::pow(3,0.5); //Lattice characteristic velocity
  const T physuC = 1.; //Physical characteristic velocity
  const T rho = 1;	//Density
  const T physNu = physuC * physL / Re;//m2/s

  //Options for blade surface boundary condition
  const bool bouzidiOn = true; //true = bouzidi, false = fullway bb

  //Time-loop options
  const int vtkIter   	   = 500; //Every 10% of max physical time
  const int statIter  	   = 10;
  const int checkIter 	   = 500;
  const int cylinderForceIter = 1;
  const int timeAvgIter    = 1000;
  const std::string checkpoint = "odd"; //load even or odd checkpoint

  //Names of output files
  std::string cylinderForceFile = "tmp/cylinderForces";

  //Characteristics needed by Grid3D
  const Characteristics<T> PhysCharacteristics(
    physL,
	  Re*physNu/physL,       //Reference velocity (m/s)
	  physNu,  //Kinematic viscosity (m2/s)
	  rho);     //Density (kg/m3)

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  // An overall indicator for all domain parts
  Vector<T,3> origin = {0.,0.,0.};
  Vector<T,3> extend = {lDomainPhysx, lDomainPhysy, lDomainPhysz};
  IndicatorCuboid3D<T> coarseDomain(extend, origin);

  // Indicator for blade
  IndicatorCylinder3D<T> cylinder(cylinderOrigin,
                                  cylinderOrigin + cylinderExtend,
                                  diameter * 0.5);

  // Construct a background coarse grid
  Grid3D<T,DESCRIPTOR> coarseGrid(
    coarseDomain,
    LatticeVelocity<T>(uC),
    N,
    PhysCharacteristics,
    false,false,true);

  //Overall domain dimensions
  const Vector<T,3> domainOrigin =
    	  coarseGrid.getSuperGeometry().getStatistics().getMinPhysR(0);
  const Vector<T,3> domainExtend =
    	  coarseGrid.getSuperGeometry().getStatistics().getPhysExtend(0);

  // === 2nd Step: Prepare Geometry ===
  prepareGeometry(coarseGrid, domainOrigin, domainExtend, cylinder);
  setupRefinement(coarseGrid, domainOrigin, domainExtend, cylinder,
                  nRefinement);

  // === 3rd Step: Prepare Lattice ===
  coarseGrid.forEachGrid(std::function<void(Grid3D<T,DESCRIPTOR>&,
    IndicatorCylinder3D<T>&, const bool&)>(prepareLattice),cylinder,bouzidiOn); 
  clout << "Total number of active cells: " << coarseGrid.getActiveVoxelN() 
	      << std::endl;

  //Reference to finest grid containing wing
  Grid3D<T,DESCRIPTOR>& cylinderGrid = coarseGrid.locate(
    Vector<T,3>(cylinderOrigin[0],cylinderOrigin[1],cylinderOrigin[2]+span/2.));

	//Functor vectors for 3D VTK
	std::vector<std::unique_ptr<SuperLatticePhysVelocity3D<T,DESCRIPTOR>>> sVel;
	std::vector<std::unique_ptr<SuperLatticePhysPressure3D<T,DESCRIPTOR>>> sP;
	std::vector<std::unique_ptr<SuperLatticePhysWallShearStressAndPressure3D<
    T,DESCRIPTOR>>> wssp;
	std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>> sAveragedWSSP;

  //Helper lambdas for initialisation and management of data
  auto initialiseVTK = [](Grid3D<T,DESCRIPTOR>& grid, auto& sVel, auto& sP,
	   auto& wssp, auto& sAveragedWSSP,
	   IndicatorCylinder3D<T>& indicatorCylinder) {
			auto& sGeometry = grid.getSuperGeometry();
			auto& sLattice = grid.getSuperLattice();
			auto& converter = grid.getConverter();

      typedef typename std::remove_reference<decltype(sVel)>::
        type::value_type sVelType; 
      typedef typename std::remove_reference<decltype(sP)>::
        type::value_type sPtype;
      typedef typename std::remove_reference<decltype(wssp)>::
        type::value_type wsspType;
      typedef typename std::remove_reference<decltype(sAveragedWSSP)>::
        type::value_type taType;

			sVel.push_back(sVelType(new typename sVelType::element_type(
        sLattice, converter)));
			sP.push_back(sPtype(new typename sPtype::element_type(
        sLattice, converter)));
			wssp.push_back(wsspType(new typename wsspType::element_type(
        sLattice,converter, sGeometry,7,indicatorCylinder)));
			sAveragedWSSP.push_back(taType(new typename taType::element_type(
        *wssp.back())));
	};

  auto loadCheckpoint = [](Grid3D<T,DESCRIPTOR>& grid, std::string&& id,
    const std::string& checkpoint, int& i_grid) {
      grid.getSuperLattice().load(id+".checkpoint");
			std::cout << checkpoint + " checkpoint loaded." << std::endl;
      i_grid++;
      id = "cylinder3d_"+std::to_string(i_grid)+"_"+checkpoint;
  };

  auto addTaEnsemble = [](Grid3D<T,DESCRIPTOR>& grid, 
				  auto& sAveragedWSSP, int& i_grid) {
					auto& sLattice = grid.getSuperLattice();
		      sLattice.communicate();
					sAveragedWSSP[i_grid]->addEnsemble();
					i_grid++;
	};
 
  auto writeTaVTK = [](Grid3D<T,DESCRIPTOR>& grid, std::string&& id, int& iT,
	  auto& wsspVector,
		auto& sAveragedWSSPVector, int& i_grid){
			auto& wssp= *wsspVector[i_grid];
			auto& sAveragedWSSP = *sAveragedWSSPVector[i_grid];
			getVTK(grid, id, iT, wssp,
        sAveragedWSSP);
			i_grid++;
      id = "cylinder3d_"+std::to_string(i_grid);
      std::cout << "Get results vtk" << endl;
	};

  auto saveCheckpoint = [](Grid3D<T,DESCRIPTOR>& grid, std::string&& id,
    const std::string& checkpoint, int& i_grid) {
      grid.getSuperLattice().save(id+".checkpoint");
			std::cout << checkpoint + " checkpoint saved." << std::endl;
      i_grid++;
      id = "cylinder3d_"+std::to_string(i_grid)+"_"+checkpoint;
  };

  //counter for multi-grid operations (TODO update grid structure to avoid this)
  int i_grid = 0;

	//Pass functor vectors and create new averaged functors for each grid 
  coarseGrid.forEachGrid(initialiseVTK, sVel, sP,
    wssp, sAveragedWSSP, cylinder);

	// === 4th Step: Main Loop with Timer ===
	clout << "starting simulation..." << endl;
	Timer<T> timer(
		coarseGrid.getConverter().getLatticeTime(maxPhysT),
		coarseGrid.getSuperGeometry().getStatistics().getNvoxel());
  timer.start();

  for ( int iT = 0; iT < 
    coarseGrid.getConverter().getLatticeTime( maxPhysT ); ++iT ) {

		// Load last checkpoint at startup
		if (iT == 0) {
      i_grid = 0;
			coarseGrid.forEachGrid(loadCheckpoint,"cylinder3d_"+
        std::to_string(i_grid)+"_"+checkpoint, checkpoint, i_grid);
		}

    // === 6th Step: Collide and Stream Execution ===
    coarseGrid.collideAndStream();

    //setBoundaryValues(coarseGrid,iT);

    // === 7th Step: Computation and Output of the Results ===
		//Add ensemble to time-averaged functors
		if ( iT % timeAvgIter == 0) {
			i_grid = 0;
			coarseGrid.forEachGrid(addTaEnsemble, sAveragedWSSP, i_grid);
		}

		//Add time-averaged functor vector into getVTK functions
    if ( iT % vtkIter == 0 ) {
			i_grid = 0;	
      coarseGrid.forEachGrid(writeTaVTK,"cylinder3d_"+std::to_string(i_grid),
        iT, wssp, sAveragedWSSP, i_grid);
		}

		// Save checkpoint
		if ( (iT % checkIter == 0) && (iT != 0)) {
			if (iT % (2 * checkIter) == 0) {
        i_grid = 0;
				coarseGrid.forEachGrid(saveCheckpoint,
          "cylinder3d_"+std::to_string(i_grid)+"_even", "even", i_grid);
			}
			else {
        i_grid = 0;
				coarseGrid.forEachGrid(saveCheckpoint,
          "cylinder3d_"+std::to_string(i_grid)+"_odd", "odd", i_grid);
			}
		}

		if ( iT % statIter == 0 ){ 
			getStats(coarseGrid, iT, timer);
			clout << "Get results stats" << endl;
   		}

		if (iT % cylinderForceIter == 0) {
			getCylinderForce(cylinderGrid,cylinderForceFile,cylinder);
      capturePressure(cylinderGrid, cylinder, iT);
		}

  }
	timer.stop();
	timer.printSummary();
}