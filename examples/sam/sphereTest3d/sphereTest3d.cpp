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

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19<>


// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( Grid3D<T,DESCRIPTOR>& grid,
		      Vector<T,3> const& origin,
                      Vector<T,3> const& extend,
                      IndicatorSphere3D<T>& indicatorSphere) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  auto& converter		= grid.getConverter();
  auto& sGeometry		= grid.getSuperGeometry();
  const T deltaX		= converter.getPhysDeltaX();

  const Vector<T,3> gridOrigin = grid.getSuperGeometry().getStatistics().getMinPhysR(0);
  const Vector<T,3> gridExtend = grid.getSuperGeometry().getStatistics().getPhysExtend(0);

  sGeometry.rename(0,1);

  //Set material number for sphere
  sGeometry.rename(1, 5, indicatorSphere);

//Set material number for inflow
 {
    const Vector<T,3> wallOrigin {origin[0] - deltaX/2.,
                    origin[1] + deltaX/2.,
                    origin[2] - deltaX/2.};
    const Vector<T,3> wallExtend {deltaX,
                    extend[1] - deltaX,
                    extend[2] + deltaX};

    IndicatorCuboid3D<T> inflow(wallExtend, wallOrigin);
    sGeometry.rename(1, 3, inflow);
  }
  //Outflow
  {
    const Vector<T,3> wallOrigin {origin[0] + extend[0] - deltaX/2.,
                    origin[1] + deltaX/2.,
                    origin[2] - deltaX/2.};
    const Vector<T,3> wallExtend {deltaX,
                    extend[1] - deltaX,
                    extend[2] + deltaX};

    IndicatorCuboid3D<T> outflow(wallExtend, wallOrigin);
    sGeometry.rename(1, 4, outflow);
  }
  //Top wall
  {
    const Vector<T,3> wallOrigin {origin[0] - deltaX/2.,
                    origin[1] + extend[1] - deltaX/2.,
                    origin[2] - deltaX/2.};
    const Vector<T,3> wallExtend {origin[0] + extend[0] + deltaX,
                    origin[1] + deltaX,
                    origin[2] + extend[2] + deltaX};
 
    IndicatorCuboid3D<T> topWall(wallExtend, wallOrigin);
    sGeometry.rename(1, 2, topWall);
  }
  //Bottom wall
  {
    const Vector<T,3> wallOrigin {origin[0] - deltaX/2.,
                    origin[1] - deltaX/2.,
                    origin[2] - deltaX/2.};
    const Vector<T,3> wallExtend {extend[0] + deltaX,
                    deltaX,
                    extend[2] + deltaX};

    IndicatorCuboid3D<T> bottomWall(wallExtend, wallOrigin);
    sGeometry.rename(1, 2, bottomWall);
  }

  //Set ghost nodes in x and y back to zero!
  //Top layer
  {
    const Vector<T,3> wallOrigin {gridOrigin[0] - 3.*deltaX/2.,
                    gridOrigin[1] + gridExtend[1] + deltaX/2.,
                    gridOrigin[2] - deltaX/2.};
    const Vector<T,3> wallExtend {gridExtend[0] + 3.*deltaX,
                    deltaX,
                    gridExtend[2] + deltaX};
 
    IndicatorCuboid3D<T> topWall(wallExtend, wallOrigin);
    sGeometry.rename(1, 0, topWall);
  }
  //Bottom layer
  {
    const Vector<T,3> wallOrigin {gridOrigin[0] - 3.*deltaX/2.,
                    gridOrigin[1] - 3.*deltaX/2.,
                    gridOrigin[2] - deltaX/2.};
    const Vector<T,3> wallExtend {gridExtend[0] + 3.*deltaX,
                    deltaX,
                    gridExtend[2] + deltaX};

    IndicatorCuboid3D<T> bottomWall(wallExtend, wallOrigin);
    sGeometry.rename(1, 0, bottomWall);
  }

  //Front layer
 {
    const Vector<T,3> wallOrigin {gridOrigin[0] - 3.*deltaX/2.,
                    gridOrigin[1] - deltaX/2.,
                    gridOrigin[2] - deltaX/2.};
    const Vector<T,3> wallExtend {deltaX,
                    gridExtend[1] + deltaX,
                    gridExtend[2] + deltaX};

    IndicatorCuboid3D<T> front(wallExtend, wallOrigin);
    sGeometry.rename(1, 0, front);
  }
  //Back layer
  {
    const Vector<T,3> wallOrigin {gridOrigin[0] + gridExtend[0] + deltaX/2.,
                    gridOrigin[1] - deltaX/2.,
                    gridOrigin[2] - deltaX/2.};
    const Vector<T,3> wallExtend {deltaX,
                    gridExtend[1] + deltaX,
                    gridExtend[2] + deltaX};

    IndicatorCuboid3D<T> back(wallExtend, wallOrigin);
    sGeometry.rename(1, 0, back);
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
                     IndicatorSphere3D<T>& indicatorSphere,
		     std::vector<T> const& sphereOrigin,
		     const T radius,
		     const int n) {

	OstreamManager clout(std::cout, "setupRefinement");
	clout << "Setup Refinement ..." << std::endl;

        //Origin of sphere bounding box
	std::vector<T> sphereBoxOrigin = {sphereOrigin[0]-radius,sphereOrigin[1]-radius,sphereOrigin[2]-radius};

	//Heights around wing box for each refinement level
  	const Vector<T,2> hn4 = {0.25*radius,0.3*2.*radius}; //x,y heights in negative direction //Innermost
  	const Vector<T,2> hp4 = {0.25*2.*radius,0.3*2.*radius}; // '' positive

  	const Vector<T,2> hn3 = {1.1*2.*radius,0.5*2.*radius};
	const Vector<T,2> hp3 = {1.1*2.*radius,0.5*2.*radius};

  	const Vector<T,2> hn2 = {2.1*2.*radius,0.9*2.*radius};
  	const Vector<T,2> hp2 = {2.1*2.*radius,0.9*2.*radius};

  	const Vector<T,2> hn1 = {3.1*2.*radius,1.3*2.*radius}; //Outermost
  	const Vector<T,2> hp1 = {3.1*2.*radius,1.3*2.*radius};

	if(n >= 1) {
		// Refinement around the wing box - level 1
	        T coarseDeltaX = coarseGrid.getConverter().getPhysDeltaX();
		Vector<T,3> fineOrigin = {sphereBoxOrigin[0]-hn1[0],sphereBoxOrigin[1]-hn1[1],domainOrigin[2]};
	 	Vector<T,3> fineExtend = {2.*radius+hp1[0]+hn1[0],2.*radius+hp1[1]+hn1[1],domainExtend[2]};

		auto& fineGrid = coarseGrid.refine(fineOrigin, fineExtend,
			false, false, true, false);
		prepareGeometry(fineGrid, domainOrigin, domainExtend, indicatorSphere);

		Vector<T,3> origin	= fineGrid.getOrigin() + Vector<T,3> {0.,0.,0.5*coarseDeltaX};
		Vector<T,3> extend	= fineGrid.getExtend() - Vector<T,3> {0.,0.,0.5*coarseDeltaX};

		Vector<T,3> extendXZ	= {extend[0], 0, extend[2]};
		Vector<T,3> extendYZ	= {0, extend[1], extend[2]};
		coarseGrid.addFineCoupling(fineGrid, origin, extendXZ);
		coarseGrid.addFineCoupling(fineGrid, origin, extendYZ);

		Vector<T,3> extendX	= {extend[0], 0, 0};
		Vector<T,3> extendY	= {0, extend[1], 0};
		coarseGrid.addFineCoupling(fineGrid, origin + extendX, extendYZ);
		coarseGrid.addFineCoupling(fineGrid, origin + extendY, extendXZ);

		Vector<T,3> innerOrigin = origin + Vector<T,3> {coarseDeltaX, coarseDeltaX, 0.};

		Vector<T,3> innerExtendXZ = {extend[0] - 2*coarseDeltaX, 0, extend[2]};
		Vector<T,3> innerExtendYZ = {0, extend[1] - 2*coarseDeltaX, extend[2]};

		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin, innerExtendXZ);
		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin, innerExtendYZ);

		Vector<T,3> innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
		Vector<T,3> innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};

		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin + innerExtendX, innerExtendYZ);
		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin + innerExtendY, innerExtendXZ);

		Vector<T,3> refinedOrigin = origin + Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX,
											-2*coarseDeltaX};
		Vector<T,3> refinedExtend = extend - Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX,
											-4*coarseDeltaX};
		IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
		coarseGrid.getSuperGeometry().reset(refined);


		if(n>=2) {
			// Refinement around the wing box - level 2
			const T deltaX0 = fineGrid.getConverter().getPhysDeltaX();
			coarseDeltaX = deltaX0;
			
			Vector<T,3> fineOrigin2 = {sphereBoxOrigin[0]-hn2[0],sphereBoxOrigin[1]-hn2[1],domainOrigin[2]-deltaX0};
	 		Vector<T,3> fineExtend2 = {2.*radius+hp2[0]+hn2[0],2.*radius+hp2[1]+hn2[1],domainExtend[2]+deltaX0};

			auto& fineGrid2 = fineGrid.refine(fineOrigin2, fineExtend2,
					false, false, true, false);
			prepareGeometry(fineGrid2, domainOrigin, domainExtend, indicatorSphere);

			origin	= fineGrid2.getOrigin() + Vector<T,3> {0.,0.,0.5*coarseDeltaX};
      			extend	= fineGrid2.getExtend() - Vector<T,3> {0.,0.,0.5*coarseDeltaX};

			extendXZ	= {extend[0], 0, extend[2]};
			extendYZ	= {0, extend[1], extend[2]};

			fineGrid.addFineCoupling(fineGrid2, origin, extendXZ);
			fineGrid.addFineCoupling(fineGrid2, origin, extendYZ);

			extendX	= {extend[0], 0, 0};
			extendY	= {0, extend[1], 0};

			fineGrid.addFineCoupling(fineGrid2, origin + extendX, extendYZ);
			fineGrid.addFineCoupling(fineGrid2, origin + extendY, extendXZ);

			innerOrigin = origin + Vector<T,3> {coarseDeltaX, coarseDeltaX, 0.};

			innerExtendXZ = {extend[0] - 2*coarseDeltaX, 0, extend[2]};
			innerExtendYZ = {0, extend[1] - 2*coarseDeltaX, extend[2]};

			fineGrid.addCoarseCoupling(fineGrid2, innerOrigin, innerExtendXZ);
			fineGrid.addCoarseCoupling(fineGrid2, innerOrigin, innerExtendYZ);

			innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
			innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};

			fineGrid.addCoarseCoupling(fineGrid2, innerOrigin + innerExtendX, innerExtendYZ);
			fineGrid.addCoarseCoupling(fineGrid2, innerOrigin + innerExtendY, innerExtendXZ);

			refinedOrigin = origin + Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX, -2*coarseDeltaX};
			refinedExtend = extend - Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX, -4*coarseDeltaX};
			IndicatorCuboid3D<T> refined2(refinedExtend, refinedOrigin);
			fineGrid.getSuperGeometry().reset(refined2);

			if(n>=3) {
			  	// Refinement around the wing box - level 3
				const T deltaX1 = fineGrid2.getConverter().getPhysDeltaX();
				coarseDeltaX = deltaX1;
			
				Vector<T,3> fineOrigin3 = {sphereBoxOrigin[0]-hn3[0],sphereBoxOrigin[1]-hn3[1],domainOrigin[2]-deltaX0-deltaX1};
	 			Vector<T,3> fineExtend3 = {2.*radius+hp3[0]+hn3[0],2.*radius+hp3[1]+hn3[1],domainExtend[2]+deltaX0+deltaX1};

				auto& fineGrid3 = fineGrid2.refine(fineOrigin3, fineExtend3,
						false, false, true, false);
				prepareGeometry(fineGrid3, domainOrigin, domainExtend, indicatorSphere);

				origin	= fineGrid3.getOrigin()+Vector<T,3>{0.,0.,0.5*coarseDeltaX};
				extend	= fineGrid3.getExtend()-Vector<T,3>{0.,0.,0.5*coarseDeltaX};;

				extendXZ	= {extend[0], 0, extend[2]};
				extendYZ	= {0, extend[1], extend[2]};
				fineGrid2.addFineCoupling(fineGrid3, origin, extendXZ);
				fineGrid2.addFineCoupling(fineGrid3, origin, extendYZ);

				extendX	= {extend[0], 0, 0};
				extendY	= {0, extend[1], 0};

				fineGrid2.addFineCoupling(fineGrid3, origin + extendX, extendYZ);
				fineGrid2.addFineCoupling(fineGrid3, origin + extendY, extendXZ);

				innerOrigin = origin + Vector<T,3> {coarseDeltaX, coarseDeltaX, 0};

				innerExtendXZ = {extend[0] - 2*coarseDeltaX, 0, extend[2]};
				innerExtendYZ = {0, extend[1] - 2*coarseDeltaX, extend[2]};
				fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin, innerExtendXZ);
				fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin, innerExtendYZ);

				innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
				innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
				fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin + innerExtendX, innerExtendYZ);
				fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin + innerExtendY, innerExtendXZ);

				refinedOrigin = origin + Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX, -2*coarseDeltaX};
				refinedExtend = extend - Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX, -4*coarseDeltaX};
				IndicatorCuboid3D<T> refined3(refinedExtend, refinedOrigin);
				fineGrid2.getSuperGeometry().reset(refined3);

				if(n>=4) {
					// Refinement around the wing box - level 4 (current innermost)

					const T deltaX2 = fineGrid3.getConverter().getPhysDeltaX();
					coarseDeltaX = deltaX2;
			
					Vector<T,3> fineOrigin4 = {sphereBoxOrigin[0]-hn4[0],sphereBoxOrigin[1]-hn4[1],
						domainOrigin[2]-deltaX0-deltaX1-deltaX2};
	 				Vector<T,3> fineExtend4 = {2.*radius+hp4[0]+hn4[0],2.*radius+hp4[1]+hn4[1],
						domainExtend[2]+deltaX0+deltaX1+deltaX2};

					auto& fineGrid4 = fineGrid3.refine(fineOrigin4, fineExtend4,
				     		 false, false, true, false);
			       		prepareGeometry(fineGrid4, domainOrigin, domainExtend, indicatorSphere);

			    		origin	= fineGrid4.getOrigin() + Vector<T,3>{0.,0.,0.5*coarseDeltaX};
			    		extend	= fineGrid4.getExtend() - Vector<T,3>{0.,0.,0.5*coarseDeltaX};

			    		extendXZ	= {extend[0], 0, extend[2]};
			    		extendYZ	= {0, extend[1], extend[2]};

			    		fineGrid3.addFineCoupling(fineGrid4, origin, extendXZ);
			    		fineGrid3.addFineCoupling(fineGrid4, origin, extendYZ);

			    		extendX	= {extend[0], 0, 0};
			    		extendY	= {0, extend[1], 0};

			    		fineGrid3.addFineCoupling(fineGrid4, origin + extendX, extendYZ);
			    		fineGrid3.addFineCoupling(fineGrid4, origin + extendY, extendXZ);

			    		innerOrigin = origin
			              		+ Vector<T,3> {coarseDeltaX, coarseDeltaX, 0.};

			    		innerExtendXZ = {extend[0] - 2*coarseDeltaX ,0, extend[2]};
			    		innerExtendYZ = {0, extend[1] - 2*coarseDeltaX, extend[2]};
			    		fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin, innerExtendXZ);
			    		fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin, innerExtendYZ);

			    		innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
			    		innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
			    		fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin + innerExtendX, innerExtendYZ);
			    		fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin + innerExtendY, innerExtendXZ);

			    		refinedOrigin = origin
			              		+ Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX,
			                      	-2*coarseDeltaX};
			    		refinedExtend = extend
			              		- Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX,
			                      	-4*coarseDeltaX};
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
		    IndicatorSphere3D<T>& indicatorSphere,
		    const bool& bouzidiOn,
	            const bool& startUpOn) {

	OstreamManager clout(std::cout, "prepareLattice");
	clout << "Prepare lattice ..." << std::endl;

	auto& converter = grid.getConverter();
	auto& sGeometry = grid.getSuperGeometry();
	auto& sLattice  = grid.getSuperLattice();
	const T omega	= converter.getLatticeRelaxationFrequency();

	// Initialize dynamics
	Dynamics<T,DESCRIPTOR>& bulkDynamics = grid.addDynamics(
			std::unique_ptr<Dynamics<T,DESCRIPTOR>>(
				new SmagorinskyBGKdynamics<T,DESCRIPTOR>(
					omega, instances::getBulkMomenta<T,DESCRIPTOR>(),0.1)));

	// Initialize boundary condition types
	//Interp
	sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc =
		grid.getOnLatticeBoundaryCondition();
	createInterpBoundaryCondition3D<T,DESCRIPTOR>(bc);

	//Local
	sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& onbc =
		grid.getOnLatticeBoundaryCondition();
	createLocalBoundaryCondition3D<T,DESCRIPTOR>(onbc);

	//Wall function
//	sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& wfbc =
//		grid.getOnLatticeBoundaryCondition();

	//Bouzidi
	sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc =
		grid.getOffLatticeBoundaryCondition();
	createBouzidiBoundaryCondition3D<T,DESCRIPTOR>(offBc);

	// Define dynamics
	sLattice.defineDynamics(sGeometry, 0, &instances::getNoDynamics<T,DESCRIPTOR>());
//	sLattice.defineDynamics(sGeometry, 2, &instances::getBounceBack<T,DESCRIPTOR>());
//	sLattice.defineDynamics(sGeometry, 2, &bulkDynamics);

	auto bulkIndicator = sGeometry.getMaterialIndicator({1, 2, 3, 4});
	sLattice.defineDynamics(bulkIndicator, &bulkDynamics);

	//sLattice.defineDynamics(sGeometry, 5, &instances::getNoDynamics<T,DESCRIPTOR>());

	// Define boundary conditions
 	onbc.addVelocityBoundary(sGeometry, 2, omega);
  	//bc.addPressureBoundary(sGeometry, 2, omega);
	onbc.addVelocityBoundary(sGeometry, 3, omega);
	bc.addPressureBoundary(sGeometry, 4, omega);

	//Options for blade surface boundary
	if ( bouzidiOn ) {
		// material=5 --> no dynamics + bouzidi zero velocity
		sLattice.defineDynamics( sGeometry,5,&instances::getNoDynamics<T,DESCRIPTOR>() );
		offBc.addZeroVelocityBoundary( sGeometry,5,indicatorSphere );
	}
	else {
		// material=5 --> fullway bounceBack dynamics
		sLattice.defineDynamics( sGeometry, 5, &instances::getBounceBack<T, DESCRIPTOR>() );
	}

	// Initial conditions - characteristic physical velocity and density for inflow
	AnalyticalConst3D<T,T> rhoF {1};
	//Vector<T,3> velocityV {0.};
	Vector<T,3> velocityV {converter.getCharLatticeVelocity(), 0., 0.};

	if (startUpOn) {
		velocityV[0] = 0.;
	}

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

	// No of time steps for smooth start-up
  //int iTmaxStart = converter.getLatticeTime( maxPhysT*0.4 );
	int iTmaxStart = 4000;
	int iTupdate = 10;

	if (iT%iTupdate == 0 && iT <= iTmaxStart) {

		// Smooth start curve, polynomial
		PolynomialStartScale<T,int> StartScale(iTmaxStart, T( 1 ));

		// Creates and sets the Poiseuille inflow profile using functors
		int iTvec[1] = {iT};
		T frac[1] = {};
		StartScale( frac,iTvec );
		std::vector<T> freestreamVelocity(3, 0);
		freestreamVelocity[0] = frac[0]*converter.getCharLatticeVelocity();
		AnalyticalConst3D<T,T> uF(freestreamVelocity);

		//Sides
		sLattice.defineU(sGeometry, 2, uF);
		//Inlet
		sLattice.defineU(sGeometry, 3, uF);

		clout << "step=" << iT << "; Freestream velocity=" << freestreamVelocity[0] << std::endl;
	}
}

// Output results to vtk files
void getVTK(Grid3D<T,DESCRIPTOR>& grid, const std::string& prefix, int iT,
	    SuperLatticeTimeAveragedF3D<T>& sAveragedVel,
	    SuperLatticeTimeAveragedF3D<T>& sAveragedP,
	    SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>& wss,
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
	SuperLatticeRefinementMetricKnudsen3D<T,DESCRIPTOR> quality(sLattice, converter);
	vtmWriter.addFunctor(geometry);
	vtmWriter.addFunctor(velocity);
	vtmWriter.addFunctor(pressure);
	vtmWriter.addFunctor(sAveragedVel);
	vtmWriter.addFunctor(sAveragedP);
	//vtmWriter.addFunctor(knudsen);
	//vtmWriter.addFunctor(quality);

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
	      SuperLatticeYplus3D<T,DESCRIPTOR>& yPlus) {

  OstreamManager clout( std::cout,"getVTK2D" );

	auto& converter = grid.getConverter();
	auto& sLattice  = grid.getSuperLattice();
	auto& sGeometry = grid.getSuperGeometry();

	//1st 2D VTK - x-y along centre span
	//Hyperplane
	auto plane1 = Hyperplane3D<T>().originAt({0.271,0.136,0.158}).normalTo({0,0,1});

	//3D functors
	SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity(sLattice, converter);
	SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure(sLattice, converter);
	SuperLatticeGeometry3D<T,DESCRIPTOR> geometry(sLattice, sGeometry);
	//SuperLatticeKnudsen3D<T,DESCRIPTOR> knudsen(sLattice);
	//SuperLatticeRefinementMetricKnudsen3D<T,DESCRIPTOR> quality(sLattice, converter);

	//Reduce 3D functors to 2D in plane1
	BlockReduction3D2D<T> plane1Velocity(velocity,plane1,BlockDataSyncMode::ReduceAndBcast,BlockDataReductionMode::Discrete);
	BlockReduction3D2D<T> plane1Pressure(pressure,plane1,BlockDataSyncMode::ReduceAndBcast,BlockDataReductionMode::Discrete);
	BlockReduction3D2D<T> plane1Geometry(geometry,plane1,BlockDataSyncMode::ReduceAndBcast,BlockDataReductionMode::Discrete);
	BlockReduction3D2D<T> plane1sAveragedVel(sAveragedVel,plane1,BlockDataSyncMode::ReduceAndBcast,BlockDataReductionMode::Discrete);
	BlockReduction3D2D<T> plane1sAveragedP(sAveragedP,plane1,BlockDataSyncMode::ReduceAndBcast,BlockDataReductionMode::Discrete);
        
	//Might be worth outputting just images instead (easier to just merge images rather than pv files?)
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
   	//auto& sGeometry = grid.getSuperGeometry();
   	auto& converter = grid.getConverter();

    // Writes output on the console
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
}

void getConvergenceStats(Grid3D<T,DESCRIPTOR>& coarseGrid, int iT,
	 util::ValueTracer<T>& tracer) {

	std::vector<T> rhoAv;
	std::vector<T> nCells;
	T rhoAvTot;


	//Collect average density and number of cells for each grid
	coarseGrid.forEachGrid([&](Grid3D<T,DESCRIPTOR>& grid, std::vector<T>& v1,
					std::vector<T>& v2) {
					v1.push_back(grid.getSuperLattice().getStatistics().getAverageRho());
	  			v2.push_back(grid.getSuperGeometry().getStatistics().getNvoxel());
	  			}, rhoAv, nCells);

					//std::cout << rhoAv.size() << rhoAv[0] << rhoAv[1] << rhoAv[2] << rhoAv[3] << rhoAv[4] << endl;
					//std::cout << nCells.size() << nCells[0] << nCells[1] << nCells[2] << nCells[3] << nCells[4] << endl;

	//Compute average over all grids
	rhoAvTot = std::inner_product(std::begin(rhoAv),std::end(rhoAv),
		std::begin(nCells),0.0)/std::accumulate(std::begin(nCells),std::end(nCells),0.0);

	tracer.takeValue(rhoAvTot, true);

	//Write to output at chosen interval
	if(iT % 1 ==0){
		ofstream myfile;
		std::string filename {"tmp/rhoAv.csv"};
		myfile.open(filename, fstream::app);
		cout << std::fixed;
		myfile << iT << "	" << std::setprecision(12)<< rhoAvTot << std::endl;
		myfile.close();
	}

	if ((iT > 100) && tracer.hasConverged()) {
		std::cout << "Simulation converged." << std::endl;
	}
}

void getBladeForce(Grid3D<T,DESCRIPTOR>& grid) {
	auto& sLattice = grid.getSuperLattice();
	auto& superGeometry = grid.getSuperGeometry();
	auto& converter = grid.getConverter();
	SuperLatticePhysDrag3D<T,DESCRIPTOR> drag( sLattice, superGeometry, 5, converter );
	T bladeForce[3];
	int input1[0];
	drag(bladeForce,input1);
	std::cout << "Lift = " << bladeForce[1] << endl;
	std::cout << "Drag = " << bladeForce[0] << endl;
}

int main( int argc, char* argv[] ) {

	std::cout << "Beginning" << endl;

	//Sphere parameters
	const T radius = 0.05*0.04523;
	const std::vector<T> sphereOrigin = {0.18092,0.17935,0.1*0.04523}; //Origin of sphere

	//Domain and simulation parameters
	const int N = 0.05*40;        // resolution of the model (coarse cells per chord)
	const int nRefinement = 2;	//Number of refinement levels (current max = 4)
	const T lDomainPhysx = 0.72368; //Length of domain in physical units (m)
	const T lDomainPhysy = 0.36184;
	const T lDomainPhysz = 0.00905; //0.2*chord
	const T maxPhysT = 10000000.; // max. simulation time in s, SI unit
	const T physL = radius; //Physical reference length (m)

	//Flow conditions
	const T Re = 100000.;       // Reynolds number
	const T Mach = 0.1;
	const T uC = Mach * 1./std::pow(3,0.5); //Lattice characteristic velocity
	const T physNu = 1.468*std::pow(10,-5); //Kinematic viscosity
	const T rho = 1.2;	//Density

	//Options for blade surface boundary condition
	const bool bouzidiOn = true; //true = bouzidi, false = fullway bb

	//Time-loop options
  	const int vtkIter   	 = 10; //Every 10% of max physical time
	const int vtk2DIter      = 1000;
  	const int statIter  	 = 10;
	const int checkIter 	 = 1000;
	const int bladeForceIter = 10;
	const bool startUpOn = false; //Slow start-up
	const std::string checkpoint = "even"; //load even or odd checkpoint

	std::cout << "2" << endl;
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

	// Indicator for sphere geometry
	IndicatorSphere3D<T> testSphere(Vector<T,3>(sphereOrigin[0],sphereOrigin[1],sphereOrigin[2]),radius);

	std::cout << "before grid" << endl;

  	// Construct a background coarse grid
  	Grid3D<T,DESCRIPTOR> coarseGrid(
		coarseDomain,
      		LatticeVelocity<T>(uC),
      		N,
      		PhysCharacteristics,
		false,false,true);

	std::cout << "after grid" << endl;

  	//Overall domain dimensions
  	const Vector<T,3> domainOrigin =
    		coarseGrid.getSuperGeometry().getStatistics().getMinPhysR(0);
  	const Vector<T,3> domainExtend =
    		coarseGrid.getSuperGeometry().getStatistics().getPhysExtend(0);

	std::cout << "before prepareGeometry" << endl;

	// === 2nd Step: Prepare Geometry ===
	prepareGeometry(coarseGrid, domainOrigin, domainExtend, testSphere);
	setupRefinement(coarseGrid, domainOrigin, domainExtend, testSphere, sphereOrigin, radius, nRefinement);

	// === 3rd Step: Prepare Lattice ===
	coarseGrid.forEachGrid(prepareLattice,testSphere,bouzidiOn,startUpOn);
	clout << "Total number of active cells: " << coarseGrid.getActiveVoxelN() << std::endl;

	//Reference to finest grid containing wing
	Grid3D<T,DESCRIPTOR>& wingGrid = coarseGrid.locate(Vector<T,3>(sphereOrigin[0],sphereOrigin[1],sphereOrigin[2]));

	//Instantiate vector of functors for time-averaged quantities, wall shear stress and wall pressure
	std::vector<std::unique_ptr<SuperLatticePhysVelocity3D<T,DESCRIPTOR>>> sVel;
	std::vector<std::unique_ptr<SuperLatticePhysPressure3D<T,DESCRIPTOR>>> sP;

	std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>> sAveragedVel;
	std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>> sAveragedP;

	std::vector<std::unique_ptr<SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>>> wss;
	//std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>> sAveragedP;

	std::vector<std::unique_ptr<SuperLatticeYplus3D<T,DESCRIPTOR>>> yPlus;

	//Pass functor vectors and create new averaged functors for each grid in similar way to 1.4 example
	int i_grid = 0;
	coarseGrid.forEachGrid([&](Grid3D<T,DESCRIPTOR>& grid,
		std::vector<std::unique_ptr<SuperLatticePhysVelocity3D<T,DESCRIPTOR>>>& sVel,
		std::vector<std::unique_ptr<SuperLatticePhysPressure3D<T,DESCRIPTOR>>>& sP,
		std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>>& sAveragedVel,
		std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>>& sAveragedP,
		std::vector<std::unique_ptr<SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>>>& wss,
		std::vector<std::unique_ptr<SuperLatticeYplus3D<T,DESCRIPTOR>>>& yPlus,
		IndicatorSphere3D<T>& indicatorSphere,
		int& i_grid) {
			auto& sGeometry = grid.getSuperGeometry();
			auto& sLattice = grid.getSuperLattice();
			auto& converter = grid.getConverter();
			sVel.push_back(std::unique_ptr<SuperLatticePhysVelocity3D<T,DESCRIPTOR>>(
				new SuperLatticePhysVelocity3D<T,DESCRIPTOR>(sLattice, converter)));
			sAveragedVel.push_back(std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>(
			new SuperLatticeTimeAveragedF3D<T>(*sVel.back())));
			sP.push_back(std::unique_ptr<SuperLatticePhysPressure3D<T,DESCRIPTOR>>(
				new SuperLatticePhysPressure3D<T,DESCRIPTOR>(sLattice, converter)));
			sAveragedP.push_back(std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>(
			new SuperLatticeTimeAveragedF3D<T>(*sP.back())));
			wss.push_back(std::unique_ptr<SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>>(
			new SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>(sLattice,sGeometry,2,converter,indicatorSphere)));
			yPlus.push_back(std::unique_ptr<SuperLatticeYplus3D<T,DESCRIPTOR>>(
			new SuperLatticeYplus3D<T,DESCRIPTOR>(sLattice,converter,sGeometry,indicatorSphere,5)));
			i_grid++;
		},sVel, sP, sAveragedVel, sAveragedP, wss, yPlus, testSphere, i_grid);

	// === 4th Step: Main Loop with Timer ===
	clout << "starting simulation..." << endl;
	Timer<T> timer(
		coarseGrid.getConverter().getLatticeTime(maxPhysT),
		coarseGrid.getSuperGeometry().getStatistics().getNvoxel() );
  	timer.start();

	// Convergence tracer every physical second
	util::ValueTracer<T> rhoTracer(/*coarseGrid.getConverter().getLatticeTime(0.01)*/1000, 1e-5);

  	for ( int iT = 0; iT < 100000/*coarseGrid.getConverter().getLatticeTime( maxPhysT )*/; ++iT ) {

		// Load last checkpoint at startup
		if (iT == 0) {
			coarseGrid.forEachGrid("cylinder3D_"+checkpoint, [&](Grid3D<T,DESCRIPTOR>& grid,
				const std::string& id)
					{grid.getSuperLattice().load(id+".checkpoint");
					 clout << checkpoint + " checkpoint loaded." << std::endl;});
		}

   		 // === 5th Step: Definition of Initial and Boundary Conditions ===
		if (startUpOn) {
    			setBoundaryValues( coarseGrid, iT);
		}

    		// === 6th Step: Collide and Stream Execution ===
    		coarseGrid.collideAndStream();

    		// === 7th Step: Computation and Output of the Results ===

		//Add ensemble to time-averaged functors
		i_grid = 0;
		coarseGrid.forEachGrid([&](Grid3D<T,DESCRIPTOR>& grid,
		        std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>>& sAveragedVel,
			std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>>& sAveragedP,
			int& i_grid) {
				auto& sLattice = grid.getSuperLattice();
		                sLattice.communicate();
				sAveragedVel[i_grid]->addEnsemble();
				sAveragedP[i_grid]->addEnsemble();
				i_grid++;
			},sAveragedVel, sAveragedP,i_grid);


		//Add time-averaged functor vector into getVTK functions
    		if ( iT % vtkIter == 0 ) {
			i_grid = 0;	
      			coarseGrid.forEachGrid("aerofoil3D",
						sAveragedVel,
						sAveragedP,
						wss,
						yPlus,
						i_grid,
						[&](Grid3D<T,DESCRIPTOR>& grid,
          					const std::string& id,
						std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>>& sAveragedVelVector,
						std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>>& sAveragedPVector,
						std::vector<std::unique_ptr<SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>>>& wssVector,
						std::vector<std::unique_ptr<SuperLatticeYplus3D<T,DESCRIPTOR>>>& yPlusVector,

						int& i_grid){
							SuperLatticeTimeAveragedF3D<T>& sAveragedVel = *sAveragedVelVector[i_grid]; 
							SuperLatticeTimeAveragedF3D<T>& sAveragedP = *sAveragedPVector[i_grid]; 
							SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>& wss= *wssVector[i_grid];
							SuperLatticeYplus3D<T,DESCRIPTOR>& yPlus = *yPlusVector[i_grid];
							getVTK(grid, id, iT, sAveragedVel, sAveragedP, wss, yPlus);
							i_grid++;
							});
      						clout << "Get results vtk" << endl;
		}

		if ( iT % vtk2DIter == 0 ) {
			i_grid = 0;
      			coarseGrid.forEachGrid("aerofoil3Dplane1",
						sAveragedVel,
						sAveragedP,
						wss,
						yPlus,
						i_grid,
						[&](Grid3D<T,DESCRIPTOR>& grid,
          					const std::string& id,
						std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>>& sAveragedVelVector,
						std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>>& sAveragedPVector,
                                            	std::vector<std::unique_ptr<SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>>>& wssVector,
                                             	std::vector<std::unique_ptr<SuperLatticeYplus3D<T,DESCRIPTOR>>>& yPlusVector,
	
						int& i_grid){
                       					SuperLatticeTimeAveragedF3D<T>& sAveragedVel = *sAveragedVelVector[i_grid];  
                       					SuperLatticeTimeAveragedF3D<T>& sAveragedP = *sAveragedPVector[i_grid];
							SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>& wss = *wssVector[i_grid];
							SuperLatticeYplus3D<T,DESCRIPTOR>& yPlus = *yPlusVector[i_grid];
							getVTK2D(grid, id, iT, sAveragedVel,sAveragedP,wss,yPlus);
							i_grid++;
							});
      						clout << "Get results vtk" << endl;
		}

		// Save checkpoint
		if ( (iT % checkIter == 0) && (iT != 0)) {
			if (iT % (2 * checkIter) == 0) {
				coarseGrid.forEachGrid("cylinder3D_even", [&](Grid3D<T,DESCRIPTOR>& grid,
					const std::string& id)
					{grid.getSuperLattice().save(id+".checkpoint");
					 clout << "Even checkpoint saved." << std::endl;});
			}
			else {
				coarseGrid.forEachGrid("cylinder3D_odd", [&](Grid3D<T,DESCRIPTOR>& grid,
					const std::string& id)
					{grid.getSuperLattice().save(id+".checkpoint");
					 clout << "Odd checkpoint saved." << std::endl;});
			}
		}

		if ( iT % statIter == 0 ){ 
			getStats(coarseGrid, iT, timer);
			clout << "Get results stats" << endl;
   		}

		if (iT % bladeForceIter == 0) {
			getBladeForce(wingGrid);
		}

		//Residuals and convergence checks
		getConvergenceStats(coarseGrid, iT, rhoTracer);
  }
	timer.stop();
	timer.printSummary();
}
