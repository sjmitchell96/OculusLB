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


// Parameters for the simulation setup
const int N = 20;        // resolution of the model
const T Re = 20.;       // Reynolds number
const T maxPhysT = 10.; // max. simulation time in s, SI unit
const T physL = 0.045; //Physical reference length (m)
const T tau = 0.53; //Lattice relaxation time

const T lDomainPhysx = 0.45; //Length of domain in physical units (m)
const T lDomainPhysy = 0.09;
const T lDomainPhysz = 0.135;

//Characteristics needed by Grid3D
const Characteristics<T> PhysCharacteristics(
		physL,
		0.2,       //Reference velocity (m/s)
		0.2*physL/Re,  //Kinematic viscosity (m2/s)
		1.225);     //Density (kg/m3)


// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( Grid3D<T,DESCRIPTOR>& grid, Vector<T,3> const& origin,
                      Vector<T,3> const& extend) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  auto& converter		= grid.getConverter();
  auto& sGeometry		= grid.getSuperGeometry();
  const T deltaX		= converter.getPhysDeltaX();

  sGeometry.rename( 0,1);

  //Set material number for inflow
  //New scope for easy copy-pasting ;)
  {
    const Vector<T,3> wallOrigin {origin[0] - deltaX/2.,
                    origin[1] + deltaX/2.,
                    origin[2] + deltaX/2.};
    const Vector<T,3> wallExtend {deltaX,
                    extend[1] - deltaX,
                    extend[2] - deltaX};

    IndicatorCuboid3D<T> inflow(wallExtend, wallOrigin);
    sGeometry.rename(1, 3, inflow);
  }
  //Outflow
  {
    const Vector<T,3> wallOrigin {origin[0] + extend[0] - deltaX/2.,
                    origin[1] + deltaX/2.,
                    origin[2] + deltaX/2.};
    const Vector<T,3> wallExtend {deltaX,
                    extend[1] - deltaX,
                    extend[2] - deltaX};

    IndicatorCuboid3D<T> outflow(wallExtend, wallOrigin);
    sGeometry.rename(1, 4, outflow);
  }
  //Top wall
  {
    const Vector<T,3> wallOrigin {origin[0] - deltaX/2.,
                    origin[1] + extend[1] - deltaX/2.,
                    origin[2] - deltaX/2.};
    const Vector<T,3> wallExtend {extend[0] + deltaX,
                    deltaX,
                    extend[2] + deltaX};

    IndicatorCuboid3D<T> topWall(wallExtend, wallOrigin);
    sGeometry.rename(1, 2, topWall);
  }
  //Bottom wall
  {
    const Vector<T,3> wallOrigin {origin[0] - deltaX/2.,
                    origin[1] + deltaX/2.,
                    origin[2] - deltaX/2.};
    const Vector<T,3> wallExtend {extend[0] + deltaX,
                    deltaX,
                    extend[2] + deltaX};

    IndicatorCuboid3D<T> bottomWall(wallExtend, wallOrigin);
    sGeometry.rename(1, 2, bottomWall);
  }
  //Left wall
  {
    const Vector<T,3> wallOrigin {origin[0] - deltaX/2.,
                    origin[1] + deltaX/2.,
                    origin[2] + extend[2] - deltaX/2.};
    const Vector<T,3> wallExtend {extend[0] + deltaX,
                    extend[1] - deltaX,
                    deltaX};

    IndicatorCuboid3D<T> leftWall(wallExtend, wallOrigin);
    sGeometry.rename(1, 2, leftWall);
  }
  //Right wall
  {
    const Vector<T,3> wallOrigin {origin[0] - deltaX/2.,
                    origin[1] + deltaX/2.,
                    origin[2] - deltaX/2.};
    const Vector<T,3> wallExtend {extend[0] + deltaX,
                    extend[1] - deltaX,
                    deltaX};

    IndicatorCuboid3D<T> rightWall(wallExtend, wallOrigin);
    sGeometry.rename(1, 2, rightWall);
  }

  // Set material number for wing

  //STLreader<T> stlReader( "lyonWing_smaller_domain.stl", 1.0 * deltaX, 0.001, 1, true );
  //sGeometry.rename(1,5,stlReader);

  //test cube
  {
    //const Vector<T,3> wallOrigin {0.09,0.045,0.045};
    //const Vector<T,3> wallExtend {0.045,0.003,0.04};
    //IndicatorCuboid3D<T> testCube(wallExtend, wallOrigin);

    STLreader<T> stlReader( "lyonWing.stl", 1.0 * deltaX, 0.001);//, 1, true );
    sGeometry.rename(1, 5, stlReader);
  }

  // Removes all not needed boundary voxels outside the surface
  sGeometry.clean();
  sGeometry.innerClean();
  sGeometry.checkForErrors();

  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void setupRefinement(Grid3D<T,DESCRIPTOR>& coarseGrid,
					 Vector<T,3> const& domainOrigin, Vector<T,3> const& domainExtend) {

	OstreamManager clout(std::cout, "setupRefinement");
	clout << "Setup Refinement ..." << std::endl;

  //Wing bounding box position (could get directly from stlReader...)
  const Vector<T,3> wingOrigin = {0.09,0.045,0.045};
  const Vector<T,3> wingExtend = {0.045,0.003,0.04};
  const T chord = 0.045;
  //Heights around wing box for each refinement level

  const Vector<T,3> hn4 = {-chord/20.,-chord/20.,-chord/20.}; //x,y,z heights in negative direction //Innermost
  const Vector<T,3> hp4 = {chord/10.,chord/20.,chord/20.}; // '' positive

  const Vector<T,3> hn3 = {-chord/10.,-chord/10.,-chord/10.};
  const Vector<T,3> hp3 = {chord/5.,chord/10.,chord/10.};

  const Vector<T,3> hn2 = {-chord/5.,-chord/5.,-chord/5.};
  const Vector<T,3> hp2 = {chord/2.5,chord/5.,chord/5.};

  const Vector<T,3> hn1 = {-chord/2.5,-chord/2.5,-chord/2.5}; //Outermost
  const Vector<T,3> hp1 = {chord/1.,chord/2.5,chord/2.5};

  // Refinement around the wing box - level 1
  const Vector<T,3> fineOrigin = wingOrigin + hn1;
  const Vector<T,3> fineExtend = wingExtend + hp1 - hn1;
  //const Vector<T,3> fineOrigin = {0.045,0.0225,0.0225};
  //const Vector<T,3> fineExtend = {0.09,0.045,0.09};

	auto& fineGrid = coarseGrid.refine(fineOrigin, fineExtend,
			false, false, false, false);
	prepareGeometry(fineGrid, domainOrigin, domainExtend);
	{
		const T coarseDeltaX = coarseGrid.getConverter().getPhysDeltaX();
		const Vector<T,3> origin	= fineGrid.getOrigin();
		const Vector<T,3> extend	= fineGrid.getExtend();

    const Vector<T,3> extendXY  = {extend[0], extend[1], 0};
		const Vector<T,3> extendXZ	= {extend[0], 0, extend[2]};
		const Vector<T,3> extendYZ	= {0, extend[1], extend[2]};

    coarseGrid.addFineCoupling(fineGrid, origin, extendXY);
		coarseGrid.addFineCoupling(fineGrid, origin, extendXZ);
		coarseGrid.addFineCoupling(fineGrid, origin, extendYZ);

		const Vector<T,3> extendX	= {extend[0], 0, 0};
		const Vector<T,3> extendY	= {0, extend[1], 0};
    const Vector<T,3> extendZ = {0, 0, extend[2]};

    coarseGrid.addFineCoupling(fineGrid, origin + extendZ, extendXY);
		coarseGrid.addFineCoupling(fineGrid, origin + extendX, extendYZ);
		coarseGrid.addFineCoupling(fineGrid, origin + extendY, extendXZ);

		const Vector<T,3> innerOrigin = origin
							+ Vector<T,3> {coarseDeltaX, coarseDeltaX, coarseDeltaX};

    const Vector<T,3> innerExtendXY = {extend[0] - 2*coarseDeltaX,
                      extend[1] - 2*coarseDeltaX, 0};
		const Vector<T,3> innerExtendXZ = {extend[0] - 2*coarseDeltaX,
										   0, extend[2] - 2*coarseDeltaX};
		const Vector<T,3> innerExtendYZ = {0, extend[1] - 2*coarseDeltaX,
										   extend[2] - 2*coarseDeltaX};

    coarseGrid.addCoarseCoupling(fineGrid, innerOrigin, innerExtendXY);
		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin, innerExtendXZ);
		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin, innerExtendYZ);

		const Vector<T,3> innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
		const Vector<T,3> innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
    const Vector<T,3> innerExtendZ = {0, 0, extend[2] - 2*coarseDeltaX};

		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin + innerExtendX, innerExtendYZ);
		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin + innerExtendY, innerExtendXZ);
    coarseGrid.addCoarseCoupling(fineGrid, innerOrigin + innerExtendZ, innerExtendXY);

		const Vector<T,3> refinedOrigin = origin
							+ Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX,
											2*coarseDeltaX};
		const Vector<T,3> refinedExtend = extend
							- Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX,
											4*coarseDeltaX};
		IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
		coarseGrid.getSuperGeometry().reset(refined);
	}

  // Refinement around the wing box - level 2
  const Vector<T,3> fineOrigin2 = wingOrigin + hn2;
	const Vector<T,3> fineExtend2 = wingExtend + hp2 - hn2;
	auto& fineGrid2 = fineGrid.refine(fineOrigin2, fineExtend2,
			false, false, false, false);
	prepareGeometry(fineGrid2, domainOrigin, domainExtend);
	{
		const T coarseDeltaX = fineGrid.getConverter().getPhysDeltaX();
		const Vector<T,3> origin	= fineGrid2.getOrigin();
		const Vector<T,3> extend	= fineGrid2.getExtend();


    const Vector<T,3> extendXY  = {extend[0], extend[1], 0};
		const Vector<T,3> extendXZ	= {extend[0], 0, extend[2]};
		const Vector<T,3> extendYZ	= {0, extend[1], extend[2]};

    fineGrid.addFineCoupling(fineGrid2, origin, extendXY);
		fineGrid.addFineCoupling(fineGrid2, origin, extendXZ);
		fineGrid.addFineCoupling(fineGrid2, origin, extendYZ);

		const Vector<T,3> extendX	= {extend[0], 0, 0};
		const Vector<T,3> extendY	= {0, extend[1], 0};
    const Vector<T,3> extendZ = {0, 0, extend[2]};

    fineGrid.addFineCoupling(fineGrid2, origin + extendZ, extendXY);
		fineGrid.addFineCoupling(fineGrid2, origin + extendX, extendYZ);
		fineGrid.addFineCoupling(fineGrid2, origin + extendY, extendXZ);

		const Vector<T,3> innerOrigin = origin
							+ Vector<T,3> {coarseDeltaX, coarseDeltaX, coarseDeltaX};

    const Vector<T,3> innerExtendXY = {extend[0] - 2*coarseDeltaX,
                      extend[1] - 2*coarseDeltaX, 0};
		const Vector<T,3> innerExtendXZ = {extend[0] - 2*coarseDeltaX,
										   0, extend[2] - 2*coarseDeltaX};
		const Vector<T,3> innerExtendYZ = {0, extend[1] - 2*coarseDeltaX,
										   extend[2] - 2*coarseDeltaX};

    fineGrid.addCoarseCoupling(fineGrid2, innerOrigin, innerExtendXY);
		fineGrid.addCoarseCoupling(fineGrid2, innerOrigin, innerExtendXZ);
		fineGrid.addCoarseCoupling(fineGrid2, innerOrigin, innerExtendYZ);

		const Vector<T,3> innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
		const Vector<T,3> innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
    const Vector<T,3> innerExtendZ = {0, 0, extend[2] - 2*coarseDeltaX};

		fineGrid.addCoarseCoupling(fineGrid2, innerOrigin + innerExtendX, innerExtendYZ);
		fineGrid.addCoarseCoupling(fineGrid2, innerOrigin + innerExtendY, innerExtendXZ);
    fineGrid.addCoarseCoupling(fineGrid2, innerOrigin + innerExtendZ, innerExtendXY);

		const Vector<T,3> refinedOrigin = origin
							+ Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX,
											2*coarseDeltaX};
		const Vector<T,3> refinedExtend = extend
							- Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX,
											4*coarseDeltaX};
		IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
		fineGrid.getSuperGeometry().reset(refined);
  }

  // Refinement around the wing box - level 3
  const Vector<T,3> fineOrigin3 = wingOrigin + hn3;
	const Vector<T,3> fineExtend3 = wingExtend + hp3 - hn3;
	auto& fineGrid3 = fineGrid2.refine(fineOrigin3, fineExtend3,
			false, false, false, false);
	prepareGeometry(fineGrid3, domainOrigin, domainExtend);
	{
		const T coarseDeltaX = fineGrid2.getConverter().getPhysDeltaX();
		const Vector<T,3> origin	= fineGrid3.getOrigin();
		const Vector<T,3> extend	= fineGrid3.getExtend();

    const Vector<T,3> extendXY  = {extend[0], extend[1], 0};
		const Vector<T,3> extendXZ	= {extend[0], 0, extend[2]};
		const Vector<T,3> extendYZ	= {0, extend[1], extend[2]};

    fineGrid2.addFineCoupling(fineGrid3, origin, extendXY);
		fineGrid2.addFineCoupling(fineGrid3, origin, extendXZ);
		fineGrid2.addFineCoupling(fineGrid3, origin, extendYZ);

		const Vector<T,3> extendX	= {extend[0], 0, 0};
		const Vector<T,3> extendY	= {0, extend[1], 0};
    const Vector<T,3> extendZ = {0, 0, extend[2]};

    fineGrid2.addFineCoupling(fineGrid3, origin + extendZ, extendXY);
		fineGrid2.addFineCoupling(fineGrid3, origin + extendX, extendYZ);
		fineGrid2.addFineCoupling(fineGrid3, origin + extendY, extendXZ);

		const Vector<T,3> innerOrigin = origin
							+ Vector<T,3> {coarseDeltaX, coarseDeltaX, coarseDeltaX};

    const Vector<T,3> innerExtendXY = {extend[0] - 2*coarseDeltaX,
                      extend[1] - 2*coarseDeltaX, 0};
		const Vector<T,3> innerExtendXZ = {extend[0] - 2*coarseDeltaX,
										   0, extend[2] - 2*coarseDeltaX};
		const Vector<T,3> innerExtendYZ = {0, extend[1] - 2*coarseDeltaX,
										   extend[2] - 2*coarseDeltaX};

    fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin, innerExtendXY);
		fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin, innerExtendXZ);
		fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin, innerExtendYZ);

		const Vector<T,3> innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
		const Vector<T,3> innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
    const Vector<T,3> innerExtendZ = {0, 0, extend[2] - 2*coarseDeltaX};

		fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin + innerExtendX, innerExtendYZ);
		fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin + innerExtendY, innerExtendXZ);
    fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin + innerExtendZ, innerExtendXY);

		const Vector<T,3> refinedOrigin = origin
							+ Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX,
											2*coarseDeltaX};
		const Vector<T,3> refinedExtend = extend
							- Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX,
											4*coarseDeltaX};
		IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
		fineGrid2.getSuperGeometry().reset(refined);
  }

  // Refinement around the wing box - level 4 (current innermost)
  const Vector<T,3> fineOrigin4 = wingOrigin + hn4;
  const Vector<T,3> fineExtend4 = wingExtend + hp4 - hn4;
  auto& fineGrid4 = fineGrid3.refine(fineOrigin4, fineExtend4,
      false, false, false, false);
  prepareGeometry(fineGrid4, domainOrigin, domainExtend);
  {
    const T coarseDeltaX = fineGrid3.getConverter().getPhysDeltaX();
    const Vector<T,3> origin	= fineGrid4.getOrigin();
    const Vector<T,3> extend	= fineGrid4.getExtend();

    const Vector<T,3> extendXY  = {extend[0], extend[1], 0};
    const Vector<T,3> extendXZ	= {extend[0], 0, extend[2]};
    const Vector<T,3> extendYZ	= {0, extend[1], extend[2]};

    fineGrid3.addFineCoupling(fineGrid4, origin, extendXY);
    fineGrid3.addFineCoupling(fineGrid4, origin, extendXZ);
    fineGrid3.addFineCoupling(fineGrid4, origin, extendYZ);

    const Vector<T,3> extendX	= {extend[0], 0, 0};
    const Vector<T,3> extendY	= {0, extend[1], 0};
    const Vector<T,3> extendZ = {0, 0, extend[2]};

    fineGrid3.addFineCoupling(fineGrid4, origin + extendZ, extendXY);
    fineGrid3.addFineCoupling(fineGrid4, origin + extendX, extendYZ);
    fineGrid3.addFineCoupling(fineGrid4, origin + extendY, extendXZ);

    const Vector<T,3> innerOrigin = origin
              + Vector<T,3> {coarseDeltaX, coarseDeltaX, coarseDeltaX};

    const Vector<T,3> innerExtendXY = {extend[0] - 2*coarseDeltaX,
                      extend[1] - 2*coarseDeltaX, 0};
    const Vector<T,3> innerExtendXZ = {extend[0] - 2*coarseDeltaX,
                       0, extend[2] - 2*coarseDeltaX};
    const Vector<T,3> innerExtendYZ = {0, extend[1] - 2*coarseDeltaX,
                       extend[2] - 2*coarseDeltaX};

    fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin, innerExtendXY);
    fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin, innerExtendXZ);
    fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin, innerExtendYZ);

    const Vector<T,3> innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
    const Vector<T,3> innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
    const Vector<T,3> innerExtendZ = {0, 0, extend[2] - 2*coarseDeltaX};

    fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin + innerExtendX, innerExtendYZ);
    fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin + innerExtendY, innerExtendXZ);
    fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin + innerExtendZ, innerExtendXY);

    const Vector<T,3> refinedOrigin = origin
              + Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX,
                      2*coarseDeltaX};
    const Vector<T,3> refinedExtend = extend
              - Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX,
                      4*coarseDeltaX};
    IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
    fineGrid3.getSuperGeometry().reset(refined);
  }

  // -----------------------
	// Refinement at the outlet half
	//const Vector<T,3> outRefineExtend {0.2,
	//							  domainExtend[1],
	//							  domainExtend[2]};
	//const Vector<T,3> outRefineOrigin {domainExtend[0] - outRefineExtend[0],
	//							  domainOrigin[1],
	//							  domainOrigin[2]};
	// add periodicity as well
	//auto& outRefineGrid = coarseGrid.refine(outRefineOrigin, outRefineExtend,
	//		false, false, true, false);
	//prepareGeometry(outRefineGrid, domainOrigin, domainExtend);
	// add couplers manually
	//{
	//	const T coarseDeltaX = coarseGrid.getConverter().getPhysDeltaX();
	//	const Vector<T,3> origin	= outRefineGrid.getOrigin()
	//			+ Vector<T,3> {0., 0., 0.5*coarseDeltaX};
	//	const Vector<T,3> extend	= outRefineGrid.getExtend()
	//			- Vector<T,3> {0., 0., 0.5*coarseDeltaX};
	//	const Vector<T,3> extendYZ	= {0, extend[1], extend[2]};
	//	coarseGrid.addFineCoupling(outRefineGrid, origin, extendYZ);

	//	const Vector<T,3> innerOrigin = origin
	//						+ Vector<T,3> {coarseDeltaX, 0, 0};
	//	coarseGrid.addCoarseCoupling(outRefineGrid, innerOrigin, extendYZ);

	//	const Vector<T,3> refinedOrigin = origin
	//			+ Vector<T,3> {2*coarseDeltaX, 0, -2*coarseDeltaX};
	//	const Vector<T,3> refinedExtend = extend
	//			- Vector<T,3> {2*coarseDeltaX, 0, -4*coarseDeltaX};
	//	IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
	//	coarseGrid.getSuperGeometry().reset(refined);
	//}

	// Refinement at the outlet half
	//const Vector<T,3> outRefineExtend2 {0.1,
	//							  domainExtend[1],
	//							  domainExtend[2] + deltaX0};
	//const Vector<T,3> outRefineOrigin2 {domainExtend[0] - outRefineExtend2[0],
	//							  domainOrigin[1],
	//							  domainOrigin[2] - deltaX0};
	// add periodicity as well
	//auto& outRefineGrid2 = outRefineGrid.refine(outRefineOrigin2, outRefineExtend2,
	//		false, false, true, false);
	//prepareGeometry(outRefineGrid2, domainOrigin, domainExtend);
	// add couplers manually
	//{
	//	const T coarseDeltaX = outRefineGrid.getConverter().getPhysDeltaX();
	//	const Vector<T,3> origin	= outRefineGrid2.getOrigin()
	//			+ Vector<T,3> {0., 0., 0.5*coarseDeltaX};
	//	const Vector<T,3> extend	= outRefineGrid2.getExtend()
	//			- Vector<T,3> {0., 0., 0.5*coarseDeltaX};
	//	const Vector<T,3> extendYZ	= {0, extend[1], extend[2]};
	//	outRefineGrid.addFineCoupling(outRefineGrid2, origin, extendYZ);

	//	const Vector<T,3> innerOrigin = origin
	//						+ Vector<T,3> {coarseDeltaX, 0, 0};
	//	outRefineGrid.addCoarseCoupling(outRefineGrid2, innerOrigin, extendYZ);

	//	const Vector<T,3> refinedOrigin = origin
	//			+ Vector<T,3> {2*coarseDeltaX, 0, -2*coarseDeltaX};
	//	const Vector<T,3> refinedExtend = extend
	//			- Vector<T,3> {2*coarseDeltaX, 0, -4*coarseDeltaX};
	//	IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
	//	outRefineGrid.getSuperGeometry().reset(refined);
	//}

	clout << "Setup Refinement ... OK" << std::endl;
}

// Create lattice structures
void prepareLattice(Grid3D<T,DESCRIPTOR>& grid) {

	OstreamManager clout(std::cout, "prepareLattice");
	clout << "Prepare lattice ..." << std::endl;

	auto& converter = grid.getConverter();
	auto& sGeometry = grid.getSuperGeometry();
	auto& sLattice  = grid.getSuperLattice();
	const T deltaX  = converter.getPhysDeltaX();
	const T omega	= converter.getLatticeRelaxationFrequency();

	// Initialize dynamics
	Dynamics<T,DESCRIPTOR>& bulkDynamics = grid.addDynamics(
			std::unique_ptr<Dynamics<T,DESCRIPTOR>>(
				new SmagorinskyBGKdynamics<T,DESCRIPTOR>(
					omega, instances::getBulkMomenta<T,DESCRIPTOR>(),0.1)));

	// Initialize boundary condition types
	sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc =
		grid.getOnLatticeBoundaryCondition();
	createLocalBoundaryCondition3D<T,DESCRIPTOR>(bc);

	sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc =
		grid.getOffLatticeBoundaryCondition();
	createBouzidiBoundaryCondition3D<T,DESCRIPTOR>(offBc);

	// Define dynamics
	sLattice.defineDynamics(sGeometry, 0, &instances::getNoDynamics<T,DESCRIPTOR>());
//	sLattice.defineDynamics(sGeometry, 2, &instances::getBounceBack<T,DESCRIPTOR>());
//	sLattice.defineDynamics(sGeometry, 2, &bulkDynamics);

	auto bulkIndicator = sGeometry.getMaterialIndicator({1, 2, 3, 4});
	sLattice.defineDynamics(bulkIndicator, &bulkDynamics);
	sLattice.defineDynamics(sGeometry, 5, &instances::getNoDynamics<T,DESCRIPTOR>());

	// Define boundary conditions
	//bc.addVelocityBoundary(sGeometry, 2, omega);
  bc.addPressureBoundary(sGeometry, 2, omega);
	bc.addVelocityBoundary(sGeometry, 3, omega);
	bc.addPressureBoundary(sGeometry, 4, omega);

  STLreader<T> stlReader( "lyonWing_smaller_domain.stl", 1.0 * deltaX, 0.001, 1, true );
	offBc.addZeroVelocityBoundary(sGeometry, 5, stlReader);

	// Initial conditions
	AnalyticalConst3D<T,T> rhoF {1.};
	Vector<T,3> velocityV {0.};
	//Vector<T,3> velocityV {converter.getCharLatticeVelocity(), 0, 0};
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
	int iTmaxStart = converter.getLatticeTime( maxPhysT*0.4 );
	int iTupdate = 30;

	if (iT%iTupdate == 0 && iT <= iTmaxStart) {
		// Smooth start curve, sinus
//		SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

		// Smooth start curve, polynomial
		PolynomialStartScale<T,int> StartScale(iTmaxStart, T( 1 ));

		// Creates and sets the Poiseuille inflow profile using functors
		int iTvec[1] = {iT};
		T frac[1] = {};
		StartScale( frac,iTvec );
		std::vector<T> maxVelocity(3, 0);
		maxVelocity[0] = 1.5*frac[0]*converter.getCharLatticeVelocity();

		T distance2Wall = converter.getConversionFactorLength()/2.;
//		RectanglePoiseuille3D<T> poiseuilleU(sGeometry, 3, maxVelocity,
//				distance2Wall, distance2Wall, distance2Wall);
		Rectangle1DPoiseuille3D<T> poiseuilleU( sGeometry, 3, maxVelocity,
				distance2Wall, 0 );
		sLattice.defineU(sGeometry, 3, poiseuilleU);

		clout << "step=" << iT << "; maxVel=" << maxVelocity[0] << std::endl;
	}
}

// Output results to vtk files
void getVTK(Grid3D<T,DESCRIPTOR>& grid, const std::string& prefix, int iT) {

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
	vtmWriter.addFunctor(knudsen);
	vtmWriter.addFunctor(quality);

	if (iT==0) {
	  vtmWriter.createMasterFile();
	}

	vtmWriter.write(iT);
}

// SPLIT INTO GET-IMAGES AND GET-STATS
// Get gnuplot images
void getImages( Grid3D<T,DESCRIPTOR>& grid, int iT) {

  auto& sLattice  = grid.getSuperLattice();
	auto& sGeometry = grid.getSuperGeometry();
	auto& converter = grid.getConverter();

  OstreamManager clout( std::cout,"getResults" );

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  //SuperLatticeYplus3D<T, DESCRIPTOR> yPlus( sLattice, converter, superGeometry, stlReader, 5 );

  //Write Image
  SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
  BlockReduction3D2D<T> planeReduction( normVel, {0, 0, 1} );
  // write output as JPEG
  heatmap::write(planeReduction, iT);
}

void getStats( Grid3D<T,DESCRIPTOR>& grid, int iT,
               Timer<T> timer) {

    auto& sLattice  = grid.getSuperLattice();
   	auto& sGeometry = grid.getSuperGeometry();
   	auto& converter = grid.getConverter();

    // Writes output on the console
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Drag, lift, pressure drop
    //AnalyticalFfromSuperF3D<T> intpolatePressure( pressure, true );
    //SuperLatticePhysDrag3D<T,DESCRIPTOR> drag( sLattice, sGeometry, 5, converter );

    //std::vector<T> point1V = sGeometry.getStatistics().getCenterPhysR( 5 );
    //std::vector<T> point2V = sGeometry.getStatistics().getCenterPhysR( 5 );
    //T point1[3] = {};
    //T point2[3] = {};
    //for ( int i = 0; i<3; i++ ) {
    //  point1[i] = point1V[i];
    //  point2[i] = point2V[i];
    //}
    //point1[0] = sGeometry.getStatistics().getMinPhysR( 5 )[0] - converter.getConversionFactorLength();
    //point2[0] = sGeometry.getStatistics().getMaxPhysR( 5 )[0] + converter.getConversionFactorLength();

    //T p1, p2;
    //intpolatePressure( &p1,point1 );
    //intpolatePressure( &p2,point2 );

    //clout << "pressure1=" << p1;
    //clout << "; pressure2=" << p2;

    //T pressureDrop = p1-p2;
    //clout << "; pressureDrop=" << pressureDrop;

    //T dragA[3];
    //int input1[0];
    //drag( dragA, input1 );
    //clout << "; drag=" << dragA[0] << "; lift=" << dragA[1] << endl;

    //int input[4] = {};
    //SuperMax3D<T> yPlusMaxF( yPlus, superGeometry, 1 );
    //T yPlusMax[1];
    //yPlusMaxF( yPlusMax,input );
    //clout << "yPlusMax=" << yPlusMax[0] << endl;
  //}
}


int main( int argc, char* argv[] ) {

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

  // Construct a background coarse grid
  Grid3D<T,DESCRIPTOR> coarseGrid(
      coarseDomain,
      RelaxationTime<T>(tau),
      N,
      PhysCharacteristics,
      false, false, true);

  //Overall domain dimensions
  const Vector<T,3> domainOrigin =
    coarseGrid.getSuperGeometry().getStatistics().getMinPhysR(0);
  const Vector<T,3> domainExtend =
    coarseGrid.getSuperGeometry().getStatistics().getPhysExtend(0);

  // === 2nd Step: Prepare Geometry ===

  //Instantiate stlReader and Layer indicators for wing
  //STLreader<T> stlReader( "lyonWing_smaller_domain.stl", 1.0 * converter.getConversionFactorLength(), 0.001, 1, true );

  prepareGeometry(coarseGrid, domainOrigin, domainExtend);

  setupRefinement(coarseGrid, domainOrigin, domainExtend);

  // === 3rd Step: Prepare Lattice ===
  coarseGrid.forEachGrid(prepareLattice);

  clout << "Total number of active cells: " << coarseGrid.getActiveVoxelN() << std::endl;

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  Timer<T> timer(
			coarseGrid.getConverter().getLatticeTime(maxPhysT),
			coarseGrid.getSuperGeometry().getStatistics().getNvoxel() );
  timer.start();

  const int vtkIter   = coarseGrid.getConverter().getLatticeTime( .5 ); //Every 0.5s physical time
  const int imageIter = coarseGrid.getConverter().getLatticeTime( .3 ); //Every 0.3s
  const int statIter  = coarseGrid.getConverter().getLatticeTime( .1 );

  for ( int iT = 0; iT < coarseGrid.getConverter().getLatticeTime( maxPhysT ); ++iT ) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( coarseGrid, iT);

    // === 6th Step: Collide and Stream Execution ===
    coarseGrid.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===

    if ( iT % vtkIter == 0 ) {
      coarseGrid.forEachGrid("aerofoil3D", [&](Grid3D<T,DESCRIPTOR>& grid,
          const std::string& id)
        {getVTK(grid, id, iT); });
      clout << "get results vtk" << endl;
    }

    //if ( iT % imageIter == 0 ) {
    //  getImages(coarseGrid, iT);
    //  clout << "get results images" << endl;
    //}

    if ( iT % statIter == 0 ) {
      getStats(coarseGrid, iT, timer);
      clout << "get results stats" << endl;
    }

  }

  timer.stop();
  timer.printSummary();
}
