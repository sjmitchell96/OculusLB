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
const int N			= 6;			// resolution of the model
const T Re			= 100.;			// Reynolds number
const T maxPhysT	= 100.;		// max. simulation time in s, SI unit
const T D			= 0.1;			// Characteristic physical lenght
const T tau			= 0.502;			// relaxation time

const T lx			= 2.5;
const T ly			= 1.6;
const T lz			= 0.1;

const T cylinderCenterX		= 0.5;
const T cylinderCenterY		= 0.8;

// Characteristics needed by Grid3D
const Characteristics<T> PhysCharacteristics(
		D,
		1.,
		1.*D/Re,
		1.0);

// Assign material numbers for boundary conditions
void prepareGeometry(Grid3D<T,DESCRIPTOR>& grid,
					 Vector<T,3> origin,
					 Vector<T,3> extend) {

	OstreamManager clout(std::cout, "prepareGeometry");
	clout << "Prepare Geometry ..." << std::endl;

	auto& converter		= grid.getConverter();
	auto& sGeometry		= grid.getSuperGeometry();
	const T deltaX		= converter.getPhysDeltaX();

	sGeometry.rename(0, 1);

	// Set material number for y-normal walls
	{
//		const Vector<T,3> wallOrigin = origin - deltaX*3./4.;
		// extend z direction to far away --- really far!
		// Up to level 5 refinement. For level 6 refinement need higher z.
		const Vector<T,3> wallOrigin {origin[0] - deltaX*3./4.,
									  origin[1] - deltaX*3./4.,
									  origin[2] - deltaX*40.};
		const Vector<T,3> wallExtend {extend[0] + deltaX*3./2.,
									  deltaX*3./2.,
									  extend[2] + deltaX*80.};

		// actually right wall
		IndicatorCuboid3D<T> rightWall(wallExtend, wallOrigin);
		sGeometry.rename(1, 2, rightWall);
	}
	{
		const Vector<T,3> wallOrigin {origin[0] - deltaX*3./4.,
									  extend[1] - deltaX*3./4.,
									  origin[2] - deltaX*40.};
		const Vector<T,3> wallExtend {extend[0] + deltaX*3./2.,
									  deltaX*3./2.,
									  extend[2] + deltaX*80.};

		// actually left wall
		IndicatorCuboid3D<T> leftWall(wallExtend, wallOrigin);
		sGeometry.rename(1, 2, leftWall);
	}

	// Set material number for inflow and outflow
	{
		const Vector<T,3> inflowOrigin {origin[0] - deltaX*3./4.,
										deltaX*3./4.,
										origin[2] - deltaX*40.};
		const Vector<T,3> inflowExtend {deltaX*3./2.,
										extend[1] - deltaX*3./2.,
										extend[2] + deltaX*80.};

		IndicatorCuboid3D<T> inflow(inflowExtend, inflowOrigin);
		sGeometry.rename(1, 3, inflow);
	}
	{
		const Vector<T,3> outflowOrigin {extend[0] - deltaX*3./4.,
										deltaX*3./4.,
										origin[2] - deltaX*40.};
		const Vector<T,3> outflowExtend {deltaX*3./2.,
										extend[1] - deltaX*3./2.,
										extend[2] + deltaX*80.};

		IndicatorCuboid3D<T> outflow(outflowExtend, outflowOrigin);
		sGeometry.rename(1, 4, outflow);
	}

	// Set material number for cylinder parallel to Z-axis
	{
		const Vector<T,3> cylinderCenter1 {cylinderCenterX, cylinderCenterY, -deltaX*40.};
		const Vector<T,3> cylinderCenter2 {cylinderCenterX, cylinderCenterY, lz+deltaX*80.};

		IndicatorCylinder3D<T> cylinder(cylinderCenter1, cylinderCenter2, D/2.);
		sGeometry.rename(1, 5, cylinder);
	}

//	while (sGeometry.getExtendedBlockGeometry(0).get(2,2,0)==1) {
//		sGeometry.communicate();
//	}

	// Clean
	sGeometry.clean();
	sGeometry.innerClean();
	sGeometry.checkForErrors();

	clout << "Prepare Geometry ... OK" << std::endl;
}

void setupRefinement(Grid3D<T,DESCRIPTOR>& coarseGrid,
					 Vector<T,3> domainOrigin, Vector<T,3> domainExtend) {

	OstreamManager clout(std::cout, "setupRefinement");
	clout << "Setup Refinement ..." << std::endl;

	// Refinement around the cylinder
	const Vector<T,3> fineExtend = {1.0, 0.4, domainExtend[2]};
	const Vector<T,3> fineOrigin = {0.2, 0.6, domainOrigin[2]};
	auto& fineGrid = coarseGrid.refine(fineOrigin, fineExtend,
			false, false, true, false);
	prepareGeometry(fineGrid, domainOrigin, domainExtend);
	{
		const T coarseDeltaX = coarseGrid.getConverter().getPhysDeltaX();
		const Vector<T,3> origin	= fineGrid.getOrigin()
										+ Vector<T,3> {0., 0., 0.5*coarseDeltaX};
		const Vector<T,3> extend	= fineGrid.getExtend()
										- Vector<T,3> {0., 0., 0.5*coarseDeltaX};
		const Vector<T,3> extendXZ	= {extend[0], 0, extend[2]};
		const Vector<T,3> extendYZ	= {0, extend[1], extend[2]};
		coarseGrid.addFineCoupling(fineGrid, origin, extendXZ);
		coarseGrid.addFineCoupling(fineGrid, origin, extendYZ);

		const Vector<T,3> extendX	= {extend[0], 0, 0};
		const Vector<T,3> extendY	= {0, extend[1], 0};
		coarseGrid.addFineCoupling(fineGrid, origin + extendX, extendYZ);
		coarseGrid.addFineCoupling(fineGrid, origin + extendY, extendXZ);

		const Vector<T,3> innerOrigin = origin
							+ Vector<T,3> {coarseDeltaX, coarseDeltaX, 0.};
		const Vector<T,3> innerExtendXZ = {extend[0] - 2*coarseDeltaX,
										   0, extend[2]};
		const Vector<T,3> innerExtendYZ = {0, extend[1] - 2*coarseDeltaX,
										   extend[2]};
		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin, innerExtendXZ);
		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin, innerExtendYZ);

		const Vector<T,3> innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
		const Vector<T,3> innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin + innerExtendX, innerExtendYZ);
		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin + innerExtendY, innerExtendXZ);

		const Vector<T,3> refinedOrigin = origin
							+ Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX,
											-2*coarseDeltaX};
		const Vector<T,3> refinedExtend = extend
							- Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX,
											-4*coarseDeltaX};
		IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
		coarseGrid.getSuperGeometry().reset(refined);
	}

	// Refinement around the cylinder
	const T deltaX0 = fineGrid.getConverter().getPhysDeltaX();
	const Vector<T,3> fineExtend2 = {0.6, 0.3, domainExtend[2] + deltaX0};
	const Vector<T,3> fineOrigin2 = {0.3, 0.65, domainOrigin[2] - deltaX0};
	auto& fineGrid2 = fineGrid.refine(fineOrigin2, fineExtend2,
			false, false, true, false);
	prepareGeometry(fineGrid2, domainOrigin, domainExtend);
	{
		const T coarseDeltaX = fineGrid.getConverter().getPhysDeltaX();
		const Vector<T,3> origin	= fineGrid2.getOrigin()
				+ Vector<T,3> {0., 0., 0.5*coarseDeltaX};
		const Vector<T,3> extend	= fineGrid2.getExtend()
				- Vector<T,3> {0., 0., 0.5*coarseDeltaX};
		const Vector<T,3> extendXZ	= {extend[0], 0, extend[2]};
		const Vector<T,3> extendYZ	= {0, extend[1], extend[2]};
		fineGrid.addFineCoupling(fineGrid2, origin, extendXZ);
		fineGrid.addFineCoupling(fineGrid2, origin, extendYZ);

		const Vector<T,3> extendX	= {extend[0], 0, 0};
		const Vector<T,3> extendY	= {0, extend[1], 0};
		fineGrid.addFineCoupling(fineGrid2, origin + extendX, extendYZ);
		fineGrid.addFineCoupling(fineGrid2, origin + extendY, extendXZ);

		const Vector<T,3> innerOrigin = origin
							+ Vector<T,3> {coarseDeltaX, coarseDeltaX, 0.};
		const Vector<T,3> innerExtendXZ = {extend[0] - 2*coarseDeltaX,
										   0, extend[2]};
	const Vector<T,3> innerExtendYZ = {0, extend[1] - 2*coarseDeltaX,
										   extend[2]};
		fineGrid.addCoarseCoupling(fineGrid2, innerOrigin, innerExtendXZ);
		fineGrid.addCoarseCoupling(fineGrid2, innerOrigin, innerExtendYZ);

		const Vector<T,3> innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
		const Vector<T,3> innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
		fineGrid.addCoarseCoupling(fineGrid2, innerOrigin + innerExtendX, innerExtendYZ);
		fineGrid.addCoarseCoupling(fineGrid2, innerOrigin + innerExtendY, innerExtendXZ);

		const Vector<T,3> refinedOrigin = origin
				+ Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX, -2*coarseDeltaX};
		const Vector<T,3> refinedExtend = extend
				- Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX, -4*coarseDeltaX};
		IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
		fineGrid.getSuperGeometry().reset(refined);
	}

	// Refinement around the cylinder
	const T deltaX1 = fineGrid2.getConverter().getPhysDeltaX();
	const Vector<T,3> fineExtend3 = {0.2, 0.2, domainExtend[2] + deltaX0 + deltaX1};
	const Vector<T,3> fineOrigin3 = {0.4, 0.7, domainOrigin[2] - deltaX0 - deltaX1};
	auto& fineGrid3 = fineGrid2.refine(fineOrigin3, fineExtend3,
			false, false, true, false);
	prepareGeometry(fineGrid3, domainOrigin, domainExtend);
	{
		const T coarseDeltaX = fineGrid2.getConverter().getPhysDeltaX();
		const Vector<T,3> origin	= fineGrid3.getOrigin()
				+ Vector<T,3> {0., 0., 0.5*coarseDeltaX};
		const Vector<T,3> extend	= fineGrid3.getExtend()
				- Vector<T,3> {0., 0., 0.5*coarseDeltaX};
		const Vector<T,3> extendXZ	= {extend[0], 0, extend[2]};
		const Vector<T,3> extendYZ	= {0, extend[1], extend[2]};
		fineGrid2.addFineCoupling(fineGrid3, origin, extendXZ);
		fineGrid2.addFineCoupling(fineGrid3, origin, extendYZ);

		const Vector<T,3> extendX	= {extend[0], 0, 0};
		const Vector<T,3> extendY	= {0, extend[1], 0};
		fineGrid2.addFineCoupling(fineGrid3, origin + extendX, extendYZ);
		fineGrid2.addFineCoupling(fineGrid3, origin + extendY, extendXZ);

		const Vector<T,3> innerOrigin = origin
							+ Vector<T,3> {coarseDeltaX, coarseDeltaX, 0.};
		const Vector<T,3> innerExtendXZ = {extend[0] - 2*coarseDeltaX,
										   0, extend[2]};
		const Vector<T,3> innerExtendYZ = {0, extend[1] - 2*coarseDeltaX,
										   extend[2]};
		fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin, innerExtendXZ);
		fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin, innerExtendYZ);

		const Vector<T,3> innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
		const Vector<T,3> innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
		fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin + innerExtendX, innerExtendYZ);
		fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin + innerExtendY, innerExtendXZ);

		const Vector<T,3> refinedOrigin = origin
				+ Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX, -2*coarseDeltaX};
		const Vector<T,3> refinedExtend = extend
				- Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX, -4*coarseDeltaX};
		IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
		fineGrid2.getSuperGeometry().reset(refined);
	}

	// Refinement at the outlet half
	const Vector<T,3> outRefineExtend {0.2,
								  domainExtend[1],
								  domainExtend[2]};
	const Vector<T,3> outRefineOrigin {domainExtend[0] - outRefineExtend[0],
								  domainOrigin[1],
								  domainOrigin[2]};
	// add periodicity as well
	auto& outRefineGrid = coarseGrid.refine(outRefineOrigin, outRefineExtend,
			false, false, true, false);
	prepareGeometry(outRefineGrid, domainOrigin, domainExtend);
	// add couplers manually
	{
		const T coarseDeltaX = coarseGrid.getConverter().getPhysDeltaX();
		const Vector<T,3> origin	= outRefineGrid.getOrigin()
				+ Vector<T,3> {0., 0., 0.5*coarseDeltaX};
		const Vector<T,3> extend	= outRefineGrid.getExtend()
				- Vector<T,3> {0., 0., 0.5*coarseDeltaX};
		const Vector<T,3> extendYZ	= {0, extend[1], extend[2]};
		coarseGrid.addFineCoupling(outRefineGrid, origin, extendYZ);

		const Vector<T,3> innerOrigin = origin
							+ Vector<T,3> {coarseDeltaX, 0, 0};
		coarseGrid.addCoarseCoupling(outRefineGrid, innerOrigin, extendYZ);

		const Vector<T,3> refinedOrigin = origin
				+ Vector<T,3> {2*coarseDeltaX, 0, -2*coarseDeltaX};
		const Vector<T,3> refinedExtend = extend
				- Vector<T,3> {2*coarseDeltaX, 0, -4*coarseDeltaX};
		IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
		coarseGrid.getSuperGeometry().reset(refined);
	}

	// Refinement at the outlet half
	const Vector<T,3> outRefineExtend2 {0.1,
								  domainExtend[1],
								  domainExtend[2] + deltaX0};
	const Vector<T,3> outRefineOrigin2 {domainExtend[0] - outRefineExtend2[0],
								  domainOrigin[1],
								  domainOrigin[2] - deltaX0};
	// add periodicity as well
	auto& outRefineGrid2 = outRefineGrid.refine(outRefineOrigin2, outRefineExtend2,
			false, false, false, false);
	prepareGeometry(outRefineGrid2, domainOrigin, domainExtend);
	// add couplers manually
	{
		const T coarseDeltaX = outRefineGrid.getConverter().getPhysDeltaX();
		const Vector<T,3> origin	= outRefineGrid2.getOrigin()
				+ Vector<T,3> {0., 0., 0.5*coarseDeltaX};
		const Vector<T,3> extend	= outRefineGrid2.getExtend()
				- Vector<T,3> {0., 0., 0.5*coarseDeltaX};
		const Vector<T,3> extendYZ	= {0, extend[1], extend[2]};
		outRefineGrid.addFineCoupling(outRefineGrid2, origin, extendYZ);

		const Vector<T,3> innerOrigin = origin
							+ Vector<T,3> {coarseDeltaX, 0, 0};
		outRefineGrid.addCoarseCoupling(outRefineGrid2, innerOrigin, extendYZ);

		const Vector<T,3> refinedOrigin = origin
				+ Vector<T,3> {2*coarseDeltaX, 0, -2*coarseDeltaX};
		const Vector<T,3> refinedExtend = extend
				- Vector<T,3> {2*coarseDeltaX, 0, -4*coarseDeltaX};
		IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
		outRefineGrid.getSuperGeometry().reset(refined);
	}

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
				new BGKdynamics<T,DESCRIPTOR>(
					omega, instances::getBulkMomenta<T,DESCRIPTOR>())));

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
	bc.addVelocityBoundary(sGeometry, 2, omega);
	bc.addVelocityBoundary(sGeometry, 3, omega);
	bc.addPressureBoundary(sGeometry, 4, omega);

	const Vector<T,3> cylinderCenter1 {cylinderCenterX, cylinderCenterY, -deltaX*10};
	const Vector<T,3> cylinderCenter2 {cylinderCenterX, cylinderCenterY, lz+deltaX*10};
	IndicatorCylinder3D<T> cylinder(cylinderCenter1, cylinderCenter2, D/2.);
	offBc.addZeroVelocityBoundary(sGeometry, 5, cylinder);

	// Initial conditions
	AnalyticalConst3D<T,T> rhoF {1.};
//	Vector<T,3> velocityV {0.};
	Vector<T,3> velocityV {converter.getCharLatticeVelocity(), 0, 0};
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
//	int iTmaxStart = converter.getLatticeTime( maxPhysT*0.4 );
	int iTmaxStart = 1000;
	int iTupdate = 10;

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

// Output to vtk files
void getResults(Grid3D<T,DESCRIPTOR>& grid, const std::string& prefix, int iT) {

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

// Measure physical forces and pressure
void takeMeasurements(Grid3D<T,DESCRIPTOR>& grid, int iT, bool print=true) {
	auto& sLattice  = grid.getSuperLattice();
	auto& sGeometry = grid.getSuperGeometry();
	auto& converter = grid.getConverter();

	SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure(sLattice, converter);
	AnalyticalFfromSuperF3D<T> intpolatePressure(pressure, true);
	SuperLatticePhysDrag3D<T,DESCRIPTOR> dragF(sLattice, sGeometry, 5, converter);

	const T point1[3] {cylinderCenterX - D/2., cylinderCenterY, lz/2.};
	const T point2[3] {cylinderCenterX + D/2., cylinderCenterY, lz/2.};

	T pressureInFrontOfCylinder, pressureBehindCylinder;
	intpolatePressure(&pressureInFrontOfCylinder, point1);
	intpolatePressure(&pressureBehindCylinder,    point2);
	const T pressureDrop = pressureInFrontOfCylinder - pressureBehindCylinder;

	const int input[4] {};
	T drag[dragF.getTargetDim()];
	dragF(drag, input);

	if (print) {
		OstreamManager clout(std::cout, "measurement");
		clout << "pressureDrop=" << pressureDrop
		      << "; drag=" << drag[0]
		      << "; lift=" << drag[1]
		      << endl;
	}
}

// Capture the pressure around the cylinder --- middle z
void capturePressure(Grid3D<T,DESCRIPTOR>& coarseGrid, int iT) {
	auto& sLattice = coarseGrid.getSuperLattice();
	auto& converter = coarseGrid.getConverter();
	const T dR = converter.getPhysDeltaX();
	SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
	AnalyticalFfromSuperF3D<T> interpolatePressure(pressure, true);
	const T radiusCylinder = D/2;

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

		myfile << pressureAtPoint << std::endl;
	}

   myfile.close();
}

int main(int argc, char* argv[]) {
	// OLB initialization
	olbInit(&argc, &argv);
	singleton::directories().setOutputDir("./tmp/");
	OstreamManager clout(std::cout, "main");

	// An overall indicator for all domain parts
	Vector<T,3> origin {0.};
	Vector<T,3> extend {lx, ly, lz};
	IndicatorCuboid3D<T> coarseDomain(extend, origin);

	// Construct a background coarse grid
	Grid3D<T,DESCRIPTOR> coarseGrid(
			coarseDomain,
			RelaxationTime<T>(tau),
			N,
			PhysCharacteristics,
			false, false, true);

	const Vector<T,3> domainOrigin =
		coarseGrid.getSuperGeometry().getStatistics().getMinPhysR(0);
	const Vector<T,3> domainExtend =
		coarseGrid.getSuperGeometry().getStatistics().getPhysExtend(0);

	// Prepare geometry ...
	prepareGeometry(coarseGrid, domainOrigin, domainExtend);

	setupRefinement(coarseGrid, domainOrigin, domainExtend);

//	prepareLattice(coarseGrid);
	coarseGrid.forEachGrid(prepareLattice);

	clout << "Total number of active cells: " << coarseGrid.getActiveVoxelN() << std::endl;
	clout << "starting simulation..." << endl;

	Timer<T> timer(
			coarseGrid.getConverter().getLatticeTime(maxPhysT),
			coarseGrid.getSuperGeometry().getStatistics().getNvoxel() );
	timer.start();

	const Vector<T,3> cylinderCenter {cylinderCenterX, cylinderCenterY, lz/2.};
	Grid3D<T,DESCRIPTOR>& cylinderGrid = coarseGrid.locate(cylinderCenter);

	// Convergence tracer every physical second
	util::ValueTracer<T> converge(coarseGrid.getConverter().getLatticeTime(0.1), 1e-5);

	const int statIter	= 10;
	const int vtkIter	= 10;
//	const int recordIter	= 2000;
	const int checkIter	= 1000;

	for (int iT = 0; iT <= coarseGrid.getConverter().getLatticeTime(maxPhysT); ++iT) {

		// Load last checkpoint at startup
		if (iT == 0) {
			coarseGrid.forEachGrid("cylinder3D_even", [&](Grid3D<T,DESCRIPTOR>& grid,
						const std::string& id)
					{grid.getSuperLattice().load(id+".checkpoint");
					 clout << "Checkpoint loaded." << std::endl;});
		}

		if ( iT % vtkIter == 0 ) {
			coarseGrid.forEachGrid("cylinder3D", [&](Grid3D<T,DESCRIPTOR>& grid,
						const std::string& id)
					{getResults(grid, id, iT); });
		}

//		setBoundaryValues(coarseGrid, iT);

//		coarseSLattice.collideAndStream();
		coarseGrid.collideAndStream();

		// Convergence check
		converge.takeValue(
				coarseGrid.getSuperLattice().getStatistics().getAverageEnergy(), true);

		// Update state
		if (iT % statIter == 0) {
			timer.update(iT);
			timer.printStep();

			takeMeasurements(cylinderGrid, iT);
		}

//		if ( iT % vtkIter == 0 ) {
//			coarseGrid.forEachGrid("cylinder3D", [&](Grid3D<T,DESCRIPTOR>& grid,
//						const std::string& id)
//					{getResults(grid, id, iT); });
//		}

//		if ( iT % recordIter == 0 ) {
//			capturePressure(coarseGrid, iT);
//		}

		// Save checkpoint
		if ( (iT % checkIter == 0) && (iT != 0) ) {
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

		if ((iT > 100) && converge.hasConverged()) {
			clout << "Simulation converged." << endl;
			break;
		}
	}

	timer.stop();
	timer.printSummary();

	// Assign material numbers
	return 0;
}
