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

//Turbulence model choice
//#define WALE
//#define Smagorinsky

#ifdef WALE
#define DESCRIPTOR WALED3Q19Descriptor
#elif defined (Smagorinsky)
#define DESCRIPTOR D3Q19<>
#else
#define DESCRIPTOR D3Q27descriptorKBCGrad
//#define DESCRIPTOR D3Q27descriptorKBC
#endif

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( Grid3D<T,DESCRIPTOR>& grid,
		      Vector<T,3> const& origin,
                      Vector<T,3> const& extend,
                      IndicatorBladeDca3D<T>& indicatorBlade) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  auto& converter  = grid.getConverter();
  auto& sGeometry  = grid.getSuperGeometry();
  const T deltaX   = converter.getPhysDeltaX();
  const T chord    = indicatorBlade.getChord();

  const Vector<T,3> bladeOrigin = indicatorBlade.getOrigin();
  const Vector<T,3> gridOrigin  = 
	  grid.getSuperGeometry().getStatistics().getMinPhysR(0);
  const Vector<T,3> gridExtend = 
          grid.getSuperGeometry().getStatistics().getPhysExtend(0);

  sGeometry.rename(0,1);
  sGeometry.rename(1, 5, indicatorBlade);

  //Material number for section of blade to read pressures
  const Vector<T,3> pressureSection1Origin {bladeOrigin[0] - chord,
	                                    bladeOrigin[1] - chord,
                                            bladeOrigin[2] + 0.02 + 
                                              0.5 * 0.2 * chord - 0.5 * deltaX};
  const Vector<T,3> pressureSection1Extend {2. * chord, 2. * chord, deltaX};
  IndicatorCuboid3D<T> pressureSection1(pressureSection1Extend,
                                        pressureSection1Origin);
  //sGeometry.rename(5, 6, pressureSection1);

  //Set material number for inflow
  {
    const Vector<T,3> wallOrigin {origin[0] - deltaX / 2.,
                                  origin[1] + deltaX / 2.,
                                  origin[2] - deltaX / 2.};
    const Vector<T,3> wallExtend {deltaX,
                                  extend[1] - deltaX,
                                  extend[2] + deltaX};
    IndicatorCuboid3D<T> inflow(wallExtend, wallOrigin);
    sGeometry.rename(1, 3, inflow);
  }
  //Outflow
  {
    const Vector<T,3> wallOrigin {origin[0] + extend[0] - deltaX / 2.,
                                  origin[1] + deltaX / 2.,
                                  origin[2] - deltaX / 2.};
    const Vector<T,3> wallExtend {deltaX,
                                  extend[1] - deltaX,
                                  extend[2] + deltaX};
    IndicatorCuboid3D<T> outflow(wallExtend, wallOrigin);
    sGeometry.rename(1, 4, outflow);
  }
  //Top wall
  {
    const Vector<T,3> wallOrigin {origin[0] - deltaX / 2.,
                                  origin[1] + extend[1] - deltaX / 2.,
                                  origin[2] - deltaX / 2.};
    const Vector<T,3> wallExtend {origin[0] + extend[0] + deltaX,
                                  origin[1] + deltaX,
                                  origin[2] + extend[2] + deltaX};
    IndicatorCuboid3D<T> topWall(wallExtend, wallOrigin);
    sGeometry.rename(1, 2, topWall);
  }
  //Bottom wall
  {
    const Vector<T,3> wallOrigin {origin[0] - deltaX / 2.,
                                  origin[1] - deltaX / 2.,
                                  origin[2] - deltaX / 2.};
    const Vector<T,3> wallExtend {extend[0] + deltaX,
                                  deltaX,
                                  extend[2] + deltaX};
    IndicatorCuboid3D<T> bottomWall(wallExtend, wallOrigin);
    sGeometry.rename(1, 2, bottomWall);
  }
  // Removes all not needed boundary voxels outside the surface
  sGeometry.clean();
  sGeometry.checkForErrors();
  sGeometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void setupRefinement(Grid3D<T,DESCRIPTOR>& coarseGrid,
		     Vector<T,3> const& domainOrigin,
		     Vector<T,3> const& domainExtend,
                     IndicatorBladeDca3D<T>& indicatorBlade,
		     const int n) {

  T chord = indicatorBlade.getChord();
  T thickness = indicatorBlade.getThickness();
  Vector<T,3> bladeOrigin = indicatorBlade.getOrigin();
  OstreamManager clout(std::cout, "setupRefinement");
  clout << "Setup Refinement ..." << std::endl;

  //Origin of sphere bounding box
  Vector<T,3> bladeBoxOrigin = {bladeOrigin[0] - chord / 2,
	                        bladeOrigin[1] - thickness / 2, 
				bladeOrigin[2]};

  //Heights around wing box for each refinement level
  //x,y heights in negative direction //Innermost
  const Vector<T,2> hn4 = {0.25 * chord, 0.3 * chord}; 
  const Vector<T,2> hp4 = {0.25 * chord, 0.3 * chord}; // '' positive

  const Vector<T,2> hn3 = {1.1 * chord, 0.5 * chord};
  const Vector<T,2> hp3 = {1.1 * chord, 0.5 * chord};

  const Vector<T,2> hn2 = {2.1 * chord, 0.9 * chord};
  const Vector<T,2> hp2 = {2.1 * chord, 0.9 * chord};

  const Vector<T,2> hn1 = {3.1 * chord, 1.3 * chord}; //Outermost
  const Vector<T,2> hp1 = {3.1 * chord, 1.3 * chord};

  if(n >= 1) {
    // Refinement around the wing box - level 1
    T coarseDeltaX = coarseGrid.getConverter().getPhysDeltaX();
    Vector<T,3> fineOrigin = 
      {bladeBoxOrigin[0] - hn1[0], bladeBoxOrigin[1] - hn1[1], domainOrigin[2]};
    Vector<T,3> fineExtend 
      = {chord+hp1[0] + hn1[0], thickness + hp1[1] + hn1[1], domainExtend[2]};

    auto& fineGrid = coarseGrid.refine(fineOrigin, fineExtend, false, false,
		                       true, false);
    prepareGeometry(fineGrid, domainOrigin, domainExtend, indicatorBlade);

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
      Vector<T,3> fineOrigin2 = {bladeBoxOrigin[0] - hn2[0],
	                         bladeBoxOrigin[1] - hn2[1],
	                         domainOrigin[2] - deltaX0};
      Vector<T,3> fineExtend2 = {chord + hp2[0] + hn2[0],
	                         thickness+hp2[1] + hn2[1],
	                         domainExtend[2] + deltaX0};
      auto& fineGrid2 = fineGrid.refine(fineOrigin2, fineExtend2, false, false,
		                        true, false);
      prepareGeometry(fineGrid2, domainOrigin, domainExtend, indicatorBlade);
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
        Vector<T,3> {2. * coarseDeltaX, 2. * coarseDeltaX, - 2. * coarseDeltaX};
      refinedExtend = extend - 
        Vector<T,3> {4. * coarseDeltaX, 4. * coarseDeltaX, - 4. * coarseDeltaX};
      IndicatorCuboid3D<T> refined2(refinedExtend, refinedOrigin);
      fineGrid.getSuperGeometry().reset(refined2);

      if(n >= 3) {
        // Refinement around the wing box - level 3
	const T deltaX1 = fineGrid2.getConverter().getPhysDeltaX();
	coarseDeltaX = deltaX1;
	Vector<T,3> fineOrigin3 = {bladeBoxOrigin[0] - hn3[0], 
		                   bladeBoxOrigin[1] - hn3[1],
				   domainOrigin[2] - deltaX0 - deltaX1};
	Vector<T,3> fineExtend3 = {chord + hp3[0] + hn3[0],
		                   thickness + hp3[1] + hn3[1],
				   domainExtend[2] + deltaX0 + deltaX1};

	auto& fineGrid3 = fineGrid2.refine(fineOrigin3, fineExtend3, false,
		                           false, true, false);
	prepareGeometry(fineGrid3, domainOrigin, domainExtend, indicatorBlade);

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
			
	  Vector<T,3> fineOrigin4 = {bladeBoxOrigin[0] - hn4[0],
		                     bladeBoxOrigin[1] - hn4[1],
				     domainOrigin[2]-deltaX0-deltaX1-deltaX2};
	  Vector<T,3> fineExtend4 = 
	    {chord + hp4[0] + hn4[0], thickness + hp4[1] + hn4[1],
             domainExtend[2] + deltaX0 + deltaX1 + deltaX2};

	  auto& fineGrid4 = fineGrid3.refine(fineOrigin4, fineExtend4,
				             false, false, true, false);
	  prepareGeometry(fineGrid4, domainOrigin, domainExtend,
			  indicatorBlade);
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
		    IndicatorBladeDca3D<T>& indicatorBlade,
		    const bool& bouzidiOn) {
  std::cout << "PREPAREL" << std::endl;
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
  #else 
    Dynamics<T,DESCRIPTOR>& bulkDynamics = 
      grid.addDynamics(std::unique_ptr<Dynamics<T,DESCRIPTOR>>(
        new KBCGradDynamics<T,DESCRIPTOR>(
          omega, instances::getBulkMomenta<T,DESCRIPTOR>())));
    //Dynamics<T,DESCRIPTOR>& bulkDynamics = 
    //  grid.addDynamics(std::unique_ptr<Dynamics<T,DESCRIPTOR>>(
    //    new KBCdynamics<T,DESCRIPTOR>(
     //     omega, instances::getBulkMomenta<T,DESCRIPTOR>())));
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

  //Bouzidi
  //sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc =
  //  grid.getOffLatticeBoundaryCondition();
  //createBouzidiBoundaryCondition3D<T,DESCRIPTOR>(offBc);

  //Grad
  sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc =
    grid.getOffLatticeBoundaryCondition();
  createGradBoundaryCondition3D<T,DESCRIPTOR>(offBc);
  //createBouzidiBoundaryCondition3D<T,DESCRIPTOR>(offBc);


  // Define dynamics
  sLattice.defineDynamics(sGeometry, 0,
		       	  &instances::getNoDynamics<T,DESCRIPTOR>());
  auto bulkIndicator = sGeometry.getMaterialIndicator({1, 2, 3, 4});
  sLattice.defineDynamics(bulkIndicator, &bulkDynamics);

  // Define boundary conditions
  onbc.addVelocityBoundary(sGeometry, 2, omega);
  onbc.addVelocityBoundary(sGeometry, 3, omega);
  bc.addPressureBoundary(sGeometry, 4, omega);

  //Options for blade surface boundary
  if ( bouzidiOn ) {
    // material=5,6 --> no dynamics + bouzidi zero velocity
    sLattice.defineDynamics( sGeometry,5,&instances::getNoDynamics<T,DESCRIPTOR>() );
    offBc.addZeroVelocityGradBoundary( sGeometry,5,indicatorBlade );
    //offBc.addZeroVelocityBoundary( sGeometry,5,indicatorBlade );
    sLattice.defineDynamics( sGeometry,6,&instances::getNoDynamics<T,DESCRIPTOR>() );
    offBc.addZeroVelocityGradBoundary( sGeometry,6,indicatorBlade );
    //offBc.addZeroVelocityBoundary( sGeometry,6,indicatorBlade );
  }
  else {
  // material=5 --> fullway bounceBack dynamics
  sLattice.defineDynamics( sGeometry, 5, &instances::getBounceBack<T, DESCRIPTOR>() );
  sLattice.defineDynamics( sGeometry, 6, &instances::getBounceBack<T, DESCRIPTOR>() );
  }
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

// Output results to vtk files
void getVTK(Grid3D<T,DESCRIPTOR>& grid, const std::string& prefix, int iT,
	    SuperLatticeTimeAveragedF3D<T>& sAveragedVel,
	    SuperLatticeTimeAveragedF3D<T>& sAveragedP,
	    SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>& wssp,
	    SuperLatticeTimeAveragedF3D<T>& sAveragedWSSP,
	    SuperLatticeYplus3D<T,DESCRIPTOR>& yPlus) {

  OstreamManager clout( std::cout,"getVTK" );

  auto& converter = grid.getConverter();
  auto& sLattice  = grid.getSuperLattice();
  auto& sGeometry = grid.getSuperGeometry();

  SuperVTMwriter3D<T> vtmWriter(prefix);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure(sLattice, converter);
  SuperLatticeGeometry3D<T,DESCRIPTOR> geometry(sLattice, sGeometry);
  SuperLatticeKnudsen3D<T,DESCRIPTOR> knudsen(sLattice);
  SuperLatticeRefinementMetricKnudsen3D<T,DESCRIPTOR> quality(sLattice,
		                                              converter);
  vtmWriter.addFunctor(geometry);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);
  vtmWriter.addFunctor(sAveragedVel);
  vtmWriter.addFunctor(sAveragedP);
  vtmWriter.addFunctor(wssp);
  vtmWriter.addFunctor(sAveragedWSSP);

  if (iT==0) {
     vtmWriter.createMasterFile();
  }

  vtmWriter.write(iT);

}

//2D VTK WRITE
void getVTK2D(Grid3D<T,DESCRIPTOR>& grid, const std::string& prefix, int iT,
	      SuperLatticeTimeAveragedF3D<T>& sAveragedVel,
	      SuperLatticeTimeAveragedF3D<T>& sAveragedP,
	      SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>& wss,
	      SuperLatticeTimeAveragedF3D<T>& sAveragedWSSP,
	      SuperLatticeYplus3D<T,DESCRIPTOR>& yPlus) {

  OstreamManager clout( std::cout,"getVTK2D" );

  auto& converter = grid.getConverter();
  auto& sLattice  = grid.getSuperLattice();
  auto& sGeometry = grid.getSuperGeometry();

  //1st 2D VTK - x-y along centre span
  //Hyperplane
  auto plane1 = 
    Hyperplane3D<T>().originAt({0.271,0.136,0.158}).normalTo({0,0,1});

  //3D functors
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure(sLattice, converter);
  SuperLatticeGeometry3D<T,DESCRIPTOR> geometry(sLattice, sGeometry);

  //Reduce 3D functors to 2D in plane1
  BlockReduction3D2D<T> plane1Velocity(velocity,plane1,
		                       BlockDataSyncMode::ReduceAndBcast,
				       BlockDataReductionMode::Discrete);
  BlockReduction3D2D<T> plane1Pressure(pressure,plane1,
		                       BlockDataSyncMode::ReduceAndBcast,
				       BlockDataReductionMode::Discrete);
  BlockReduction3D2D<T> plane1Geometry(geometry,plane1,
		                       BlockDataSyncMode::ReduceAndBcast,
				       BlockDataReductionMode::Discrete);
  BlockReduction3D2D<T> plane1sAveragedVel(sAveragedVel,plane1,
		                           BlockDataSyncMode::ReduceAndBcast,
					   BlockDataReductionMode::Discrete);
  BlockReduction3D2D<T> plane1sAveragedP(sAveragedP,plane1,
		                         BlockDataSyncMode::ReduceAndBcast,
					 BlockDataReductionMode::Discrete);
        
  //Might be worth outputting just images instead 
  //2D vtm writer
  BlockVTKwriter2D<T> vtkWriter(prefix);
  vtkWriter.addFunctor(plane1Velocity);
  vtkWriter.addFunctor(plane1Pressure);
  vtkWriter.addFunctor(plane1Geometry);
  vtkWriter.addFunctor(plane1sAveragedVel);
  vtkWriter.addFunctor(plane1sAveragedP);
  vtkWriter.write(iT);
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

void getBladeForce(Grid3D<T,DESCRIPTOR>& grid,
		   const std::string& filePath,
		   IndicatorBladeDca3D<T>& blade) {
  auto& sLattice = grid.getSuperLattice();
  auto& superGeometry = grid.getSuperGeometry();
  auto& converter = grid.getConverter();
  T theta = blade.getTheta();
  T span = blade.getSpan();
  T chord = blade.getChord();

  SuperLatticePhysDragBlade3D<T,DESCRIPTOR> drag(sLattice, superGeometry,
		                                 std::vector<int>{5,6},
						 converter,span*chord);
  T bladeForce[3];
  int input1[0];
  drag(bladeForce,input1);
  std::cout << "Lift = " << bladeForce[1] << endl;
  std::cout << "Drag = " << bladeForce[0] << endl;

  ofstream myfile;
  std::string filename {filePath+".csv"};
  myfile.open(filename,fstream::app);
  myfile << theta << "	" << bladeForce[1] << "	" << bladeForce[0] << std::endl;
  myfile.close();
}	

int main( int argc, char* argv[] ) {

  //Blade parameters
  const Vector<T,3> bladeOrigin = {0.2295,0.204,-0.02}; //Origin of blade
  const T chord = 0.051;
  const T thickness = 0.00382;
  const T span = 0.04;
  const T r1 = 0.1836;
  const T r2 = 0.00015;
  const T xp = 0.02538;
  const T theta = 0.00; //Pitch (+ve = anticlockwise)

  //Domain and simulation parameters
  const int N = 14;        // resolution of the model (coarse cells per chord)
  const int nRefinement = 4;	//Number of refinement levels (current max = 4)
  const T lDomainPhysx = 16.*chord; //Length of domain in physical units (m)
  const T lDomainPhysy = 8.*chord;
  const T lDomainPhysz = 0.2*chord; //
  const T maxPhysT = 100; // max. simulation time in s, SI unit
  const T physL = chord; //Physical reference length (m)

  //Flow conditions
  const T Re = 100000.;       // Reynolds number
  const T Mach = 0.1;
  const T uC = Mach * 1./std::pow(3,0.5); //Lattice characteristic velocity
  const T physNu = 1.468*std::pow(10,-5); //Kinematic viscosity
  const T rho = 1.2;	//Density

  //Options for blade surface boundary condition
  const bool bouzidiOn = true; //true = bouzidi, false = fullway bb

  //Time-loop options
  const int vtkIter   	   = 10; //Every 10% of max physical time
  //const int vtk2DIter      = 20;
  const int statIter  	   = 10;
  const int checkIter 	   = 1000;
  const int bladeForceIter = 1;
  const int timeAvgIter    = 1;
  const std::string checkpoint = "odd"; //load even or odd checkpoint

  //Names of output files
  std::string bladeForceFile = "tmp/bladeForces";

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
  IndicatorBladeDca3D<T> blade(bladeOrigin,chord, thickness, span, r1, r2, xp,
    theta);

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
  prepareGeometry(coarseGrid, domainOrigin, domainExtend, blade);
  setupRefinement(coarseGrid, domainOrigin, domainExtend, blade, nRefinement);

  // === 3rd Step: Prepare Lattice ===
  coarseGrid.forEachGrid(std::function<void(Grid3D<T,DESCRIPTOR>&,
    IndicatorBladeDca3D<T>&, const bool&)>(prepareLattice),blade,bouzidiOn); 
  clout << "Total number of active cells: " << coarseGrid.getActiveVoxelN() 
	      << std::endl;

  //Reference to finest grid containing wing
  Grid3D<T,DESCRIPTOR>& wingGrid = coarseGrid.locate(
    Vector<T,3>(bladeOrigin[0],bladeOrigin[1],bladeOrigin[2]+span/2.));

	//Instantiate vector of functors for time-averaged quantities
	std::vector<std::unique_ptr<SuperLatticePhysVelocity3D<T,DESCRIPTOR>>> sVel;
	std::vector<std::unique_ptr<SuperLatticePhysPressure3D<T,DESCRIPTOR>>> sP;
	std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>> sAveragedVel;
	std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>> sAveragedP;
	std::vector<std::unique_ptr<SuperLatticePhysWallShearStressAndPressure3D<
    T,DESCRIPTOR>>> wssp;
	std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>> sAveragedWSSP;
	std::vector<std::unique_ptr<SuperLatticeYplus3D<T,DESCRIPTOR>>> yPlus;

  //Helper lambdas for initialisation and management of data
  auto initialiseVTK = [](Grid3D<T,DESCRIPTOR>& grid, auto& sVel, auto& sP,
	  auto& sAveragedVel, auto& sAveragedP, auto& wssp, auto& sAveragedWSSP,
	  auto& yPlus, IndicatorBladeDca3D<T>& indicatorBlade) {
			auto& sGeometry = grid.getSuperGeometry();
			auto& sLattice = grid.getSuperLattice();
			auto& converter = grid.getConverter();

      typedef typename std::remove_reference<decltype(sVel)>::
        type::value_type sVelType; 
      typedef typename std::remove_reference<decltype(sP)>::
        type::value_type sPtype;
      typedef typename std::remove_reference<decltype(wssp)>::
        type::value_type wsspType;
      typedef typename std::remove_reference<decltype(yPlus)>::
        type::value_type yPlusType;
      typedef typename std::remove_reference<decltype(sAveragedP)>::
        type::value_type taType;

			sVel.push_back(sVelType(new typename sVelType::element_type(
        sLattice, converter)));
			sAveragedVel.push_back(taType(new typename taType::element_type(
        *sVel.back())));
			sP.push_back(sPtype(new typename sPtype::element_type(
        sLattice, converter)));
			sAveragedP.push_back(taType(new typename taType::element_type(
        *sP.back())));
			wssp.push_back(wsspType(new typename wsspType::element_type(
        sLattice,converter, sGeometry,6,indicatorBlade)));
			yPlus.push_back(yPlusType(new typename yPlusType::element_type(
        sLattice,converter,sGeometry,indicatorBlade,6)));
			sAveragedWSSP.push_back(taType(new typename taType::element_type(
        *wssp.back())));
	};

  auto loadCheckpoint = [](Grid3D<T,DESCRIPTOR>& grid, std::string&& id,
    const std::string& checkpoint, int& i_grid) {
      grid.getSuperLattice().load(id+".checkpoint");
			std::cout << checkpoint + " checkpoint loaded." << std::endl;
      i_grid++;
      id = "dcaBlade3d_"+std::to_string(i_grid)+"_"+checkpoint;
  };

  auto addTaEnsemble = [](Grid3D<T,DESCRIPTOR>& grid, auto& sAveragedVel,
				auto& sAveragedP, auto& sAveragedWSSP, int& i_grid) {
					auto& sLattice = grid.getSuperLattice();
		      sLattice.communicate();
					sAveragedVel[i_grid]->addEnsemble();
					sAveragedP[i_grid]->addEnsemble();
					sAveragedWSSP[i_grid]->addEnsemble();
					i_grid++;
	};
 
  auto writeTaVTK = [](Grid3D<T,DESCRIPTOR>& grid, std::string&& id, int& iT,
	  auto& sAveragedVelVector, auto& sAveragedPVector, auto& wsspVector,
		auto& sAveragedWSSPVector, auto& yPlusVector, int& i_grid){
		  auto& sAveragedVel = *sAveragedVelVector[i_grid]; 
			auto& sAveragedP = *sAveragedPVector[i_grid]; 
			auto& wssp= *wsspVector[i_grid];
			auto& sAveragedWSSP = *sAveragedWSSPVector[i_grid];
			auto& yPlus = *yPlusVector[i_grid];
			getVTK(grid, id, iT, sAveragedVel, sAveragedP, wssp,
        sAveragedWSSP, yPlus);
			i_grid++;
      id = "dcaBlade3d_"+std::to_string(i_grid);
      std::cout << "Get results vtk" << endl;
	};

  auto saveCheckpoint = [](Grid3D<T,DESCRIPTOR>& grid, std::string&& id,
    const std::string& checkpoint, int& i_grid) {
      grid.getSuperLattice().save(id+".checkpoint");
			std::cout << checkpoint + " checkpoint saved." << std::endl;
      i_grid++;
      id = "dcaBlade3d_"+std::to_string(i_grid)+"_"+checkpoint;
  };

  //counter for multi-grid operations (TODO update grid structure to avoid this)
  int i_grid = 0;

	//Pass functor vectors and create new averaged functors for each grid 
  coarseGrid.forEachGrid(initialiseVTK, sVel, sP, sAveragedVel, sAveragedP,
    wssp, sAveragedWSSP, yPlus, blade);

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
			coarseGrid.forEachGrid(loadCheckpoint,"dcaBlade3d_"+
        std::to_string(i_grid)+"_"+checkpoint, checkpoint, i_grid);
		}

    // === 6th Step: Collide and Stream Execution ===
    coarseGrid.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
		//Add ensemble to time-averaged functors
		if ( iT % timeAvgIter == 0) {
			i_grid = 0;
			coarseGrid.forEachGrid(addTaEnsemble, sAveragedVel, sAveragedP,
       sAveragedWSSP, i_grid);
		}

		//Add time-averaged functor vector into getVTK functions
    if ( iT % vtkIter == 0 ) {
			i_grid = 0;	
      coarseGrid.forEachGrid(writeTaVTK,"dcaBlade3d_"+std::to_string(i_grid),
        iT, sAveragedVel, sAveragedP, wssp, sAveragedWSSP, yPlus, i_grid);
		}

		// Save checkpoint
		if ( (iT % checkIter == 0) && (iT != 0)) {
			if (iT % (2 * checkIter) == 0) {
        i_grid = 0;
				coarseGrid.forEachGrid(saveCheckpoint,
          "dcaBlade3d_"+std::to_string(i_grid)+"_even", "even", i_grid);
			}
			else {
        i_grid = 0;
				coarseGrid.forEachGrid(saveCheckpoint,
          "dcaBlade3d_"+std::to_string(i_grid)+"_odd", "odd", i_grid);
			}
		}

		if ( iT % statIter == 0 ){ 
			getStats(coarseGrid, iT, timer);
			clout << "Get results stats" << endl;
   		}

		if (iT % bladeForceIter == 0) {
			getBladeForce(wingGrid,bladeForceFile,blade);
		}

  }
	timer.stop();
	timer.printSummary();
}