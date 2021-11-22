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
const int N = 10;        // resolution of the model
const T Re = 50000.;       // Reynolds number
const T maxPhysT = 1.; // max. simulation time in s, SI unit
const T physL = 0.04523; //Physical reference length (m)
const T uC = 0.1 * 1./std::pow(3,0.5); //Lattice characteristic velocity
const T physNu = 1.468*std::pow(10,-5); //Kinematic viscosity

const T lDomainPhysx = 12.*physL; //Length of domain in physical units (m)
const T lDomainPhysy = 6.*physL;
const T lDomainPhysz = 7.*physL;
const int nRefinement = 3;	//Number of refinement levels (current max = 4)
const bool bouzidiOn = true; //true = bouzidi, false = fullway bb
const bool loadCheckpoint = true; //Load checkpoint at startup
const bool saveCheckpoint = true; //Save checkpoints at intervals
const bool startUpOn = true; //Slow start-up

//Characteristics needed by Grid3D
const Characteristics<T> PhysCharacteristics(
		physL,
		Re*physNu/physL,       //Reference velocity (m/s)
		physNu,  //Kinematic viscosity (m2/s)
		1.225);     //Density (kg/m3)


// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( Grid3D<T,DESCRIPTOR>& grid, Vector<T,3> const& origin,
                      Vector<T,3> const& extend, STLreader<T>& stlReader,
											IndicatorLayer3D<T>& wingLayer,
										  STLreader<T>& stlReader2) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  auto& converter		= grid.getConverter();
  auto& sGeometry		= grid.getSuperGeometry();
  const T deltaX		= converter.getPhysDeltaX();

  sGeometry.rename( 0,1);

//Outer domain ---
//Set material number for inflow
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

  //Wing and surrounding 'box' - similar style to cylinder example
  sGeometry.rename(1, 5, wingLayer);
	sGeometry.rename(5, 1, stlReader);

	//Origin and extend of stl outer box
	const T chord = 0.04523;
	const T height = 0.003137;
	const T span = 0.0400;
	const T rShaft = 0.001; //Shaft radius, from cad file
	const Vector<T,3> boxOrigin {3.*chord, 0.5*extend[1]-0.5*height-chord,
									0.5*extend[2]-0.5*span-chord};
	const Vector<T,3> boxExtend {3*chord, 2*chord+height, 2*chord+span};

	const Vector<T,3> layerOrigin {boxOrigin[0] -deltaX, boxOrigin[1]-deltaX,
									boxOrigin[2]-deltaX};
	const Vector<T,3> layerExtend {boxExtend[0]+2*deltaX, boxExtend[1]+2*deltaX,
									boxExtend[2]+2*deltaX};

	//Set outer box of stl file back to fluid
	{
     const Vector<T,3> wallOrigin {layerOrigin[0] - deltaX/2.,
                     layerOrigin[1] + deltaX/2.,
                     layerOrigin[2] + deltaX/2.};
     const Vector<T,3> wallExtend {deltaX,
                     layerExtend[1] - deltaX,
                     layerExtend[2] - deltaX};

     IndicatorCuboid3D<T> inflow(wallExtend, wallOrigin);
     sGeometry.rename(5, 1, inflow);
   }
   //Outflow
   {
     const Vector<T,3> wallOrigin {layerOrigin[0] + layerExtend[0] - deltaX/2.,
                     layerOrigin[1] + deltaX/2.,
                     layerOrigin[2] + deltaX/2.};
     const Vector<T,3> wallExtend {deltaX,
                     layerExtend[1] - deltaX,
                     layerExtend[2] - deltaX};

     IndicatorCuboid3D<T> outflow(wallExtend, wallOrigin);
     sGeometry.rename(5, 1, outflow);
   }
   //Top wall
   {
     const Vector<T,3> wallOrigin {layerOrigin[0] - deltaX/2.,
                     layerOrigin[1] + layerExtend[1] - deltaX/2.,
                     layerOrigin[2] - deltaX/2.};
     const Vector<T,3> wallExtend {layerOrigin[0] + layerExtend[0] + deltaX,
                     layerOrigin[1] + deltaX,
                     layerOrigin[2] + layerExtend[2] + deltaX};

     IndicatorCuboid3D<T> topWall(wallExtend, wallOrigin);
     sGeometry.rename(5, 1, topWall);
   }
   //Bottom wall
   {
     const Vector<T,3> wallOrigin {layerOrigin[0] - deltaX/2.,
                     layerOrigin[1] - deltaX/2.,
                     layerOrigin[2] - deltaX/2.};
     const Vector<T,3> wallExtend {layerExtend[0] + deltaX,
                     deltaX,
                     layerExtend[2] + deltaX};

     IndicatorCuboid3D<T> bottomWall(wallExtend, wallOrigin);
     sGeometry.rename(5, 1, bottomWall);
   }
   //Left wall
   {
     const Vector<T,3> wallOrigin {layerOrigin[0] - deltaX/2.,
                     layerOrigin[1] + deltaX/2.,
                     layerOrigin[2] + layerExtend[2] - deltaX/2.};
     const Vector<T,3> wallExtend {layerExtend[0] + deltaX,
                     layerExtend[1] - deltaX,
                     deltaX};

     IndicatorCuboid3D<T> leftWall(wallExtend, wallOrigin);
     sGeometry.rename(5, 1, leftWall);
   }
   //Right wall
   {
     const Vector<T,3> wallOrigin {layerOrigin[0] - deltaX/2.,
                     layerOrigin[1] + deltaX/2.,
                     layerOrigin[2] - deltaX/2.};
     const Vector<T,3> wallExtend {layerExtend[0] + deltaX,
                     layerExtend[1] - deltaX,
                     deltaX};

     IndicatorCuboid3D<T> rightWall(wallExtend, wallOrigin);
     sGeometry.rename(5, 1, rightWall);
   }

	 //Clean shaft section
	 {
     const Vector<T,3> shaftOrigin {4*chord+0.5*chord-1.1*rShaft,
                     0.5*extend[1]-1.1*rShaft,
                     0.5*extend[2]-0.5*span-chord-0.5*deltaX};
     const Vector<T,3> shaftExtend {2.2*rShaft,
                     2.2*rShaft,
                     chord+deltaX};

     IndicatorCuboid3D<T> shaft(shaftExtend, shaftOrigin);
     sGeometry.rename(5, 1, shaft);
   }
	 //Patch shaft hole (set to boundary)
	 {
			const Vector<T,3> shaftPatchOrigin {4*chord+0.5*chord-1.2*rShaft,
											0.5*extend[1]-1.2*rShaft,
											0.5*extend[2]-0.5*span};
			const Vector<T,3> shaftPatchExtend {2.4*rShaft,
											2.4*rShaft,
											deltaX};

			IndicatorCuboid3D<T> shaftPatch(shaftPatchExtend, shaftPatchOrigin);
			sGeometry.rename(1, 5, shaftPatch);
		}

	 //Reset wing interior to zero... how ... another stl??
	 sGeometry.rename(1,0, stlReader2);

  // Removes all not needed boundary voxels outside the surface
  sGeometry.clean();
  sGeometry.innerClean();
  sGeometry.checkForErrors();

  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void setupRefinement(Grid3D<T,DESCRIPTOR>& coarseGrid,
					 Vector<T,3> const& domainOrigin, Vector<T,3> const& domainExtend,
           STLreader<T>& stlReader, IndicatorLayer3D<T>& wingLayer,
					 STLreader<T>& stlReader2, const int n) {

	OstreamManager clout(std::cout, "setupRefinement");
	clout << "Setup Refinement ..." << std::endl;

	//Wing box dimensions
  const T chord = 0.04523;
	const T height = 0.003137;
	const T span = 0.04000;

	//Wing box positions at resting point (zero pitch)
	const Vector<T,3> wingMin = {4*chord,0.5*domainExtend[1]-0.5*height,0.5*domainExtend[2]-0.5*span};
	const Vector<T,3> wingMax = {4*chord+chord,0.5*domainExtend[1]-0.5*height+height,0.5*domainExtend[2]-0.5*span+span};

  //Heights around wing box for each refinement level
  const Vector<T,3> hn4 = {0.05*chord,0.2*chord,0.05*chord}; //x,y,z heights in negative direction //Innermost
  const Vector<T,3> hp4 = {0.05*chord,0.2*chord,0.05*chord}; // '' positive

  const Vector<T,3> hn3 = {0.1*chord,0.25*chord,0.1*chord};
  const Vector<T,3> hp3 = {0.1*chord,0.25*chord,0.1*chord};

  const Vector<T,3> hn2 = {0.2*chord,0.35*chord,0.2*chord};
  const Vector<T,3> hp2 = {0.5*chord,0.35*chord,0.2*chord};

  const Vector<T,3> hn1 = {0.4*chord,0.55*chord,0.4*chord}; //Outermost
  const Vector<T,3> hp1 = {1.0*chord,0.55*chord,0.4*chord};

	if(n >= 1) {
	  // Refinement around the wing box - level 1
	  Vector<T,3> fineOrigin = wingMin - hn1;
	  Vector<T,3> fineExtend = wingMax + hp1 - (wingMin - hn1);
		T coarseDeltaX = coarseGrid.getConverter().getPhysDeltaX();

		//Check if aligned with coarse grid. If not, round
		for(int i=0;i<3;++i) {
			T remainder = std::fmod(fineOrigin[i],coarseDeltaX);
			if(remainder != 0) {
				fineOrigin[i] -= remainder;
				cout << "Rounding fineOrigin " << i << " by " << remainder << endl;
				cout << "New remainder = " << std::fmod(fineOrigin[i],coarseDeltaX) << endl;
			}
			remainder = std::fmod(fineExtend[i],coarseDeltaX);
			if(remainder != 0) {
				fineExtend[i] -= remainder;
				cout << "Rounding fineExtend " << i << " by " << remainder << endl;
				cout << "New remainder = " << std::fmod(fineExtend[i],coarseDeltaX) << endl;
			}
		}

		auto& fineGrid = coarseGrid.refine(fineOrigin, fineExtend,
				false, false, false, false);
		prepareGeometry(fineGrid, domainOrigin, domainExtend, stlReader, wingLayer, stlReader2);


		Vector<T,3> origin	= fineGrid.getOrigin();
		Vector<T,3> extend	= fineGrid.getExtend();

	  Vector<T,3> extendXY  = {extend[0], extend[1], 0};
		Vector<T,3> extendXZ	= {extend[0], 0, extend[2]};
		Vector<T,3> extendYZ	= {0, extend[1], extend[2]};

	  coarseGrid.addFineCoupling(fineGrid, origin, extendXY);
		coarseGrid.addFineCoupling(fineGrid, origin, extendXZ);
		coarseGrid.addFineCoupling(fineGrid, origin, extendYZ);

		Vector<T,3> extendX	= {extend[0], 0, 0};
		Vector<T,3> extendY	= {0, extend[1], 0};
	  Vector<T,3> extendZ = {0, 0, extend[2]};

	  coarseGrid.addFineCoupling(fineGrid, origin + extendZ, extendXY);
		coarseGrid.addFineCoupling(fineGrid, origin + extendX, extendYZ);
		coarseGrid.addFineCoupling(fineGrid, origin + extendY, extendXZ);

		Vector<T,3> innerOrigin = origin
							+ Vector<T,3> {coarseDeltaX, coarseDeltaX, coarseDeltaX};

	  Vector<T,3> innerExtendXY = {extend[0] - 2*coarseDeltaX,
	                    extend[1] - 2*coarseDeltaX, 0};
		Vector<T,3> innerExtendXZ = {extend[0] - 2*coarseDeltaX,
										   0, extend[2] - 2*coarseDeltaX};
		Vector<T,3> innerExtendYZ = {0, extend[1] - 2*coarseDeltaX,
										   extend[2] - 2*coarseDeltaX};

	  coarseGrid.addCoarseCoupling(fineGrid, innerOrigin, innerExtendXY);
		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin, innerExtendXZ);
		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin, innerExtendYZ);

		Vector<T,3> innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
		Vector<T,3> innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
	  Vector<T,3> innerExtendZ = {0, 0, extend[2] - 2*coarseDeltaX};

		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin + innerExtendX, innerExtendYZ);
		coarseGrid.addCoarseCoupling(fineGrid, innerOrigin + innerExtendY, innerExtendXZ);
	  coarseGrid.addCoarseCoupling(fineGrid, innerOrigin + innerExtendZ, innerExtendXY);

		Vector<T,3> refinedOrigin = origin
							+ Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX,
											2*coarseDeltaX};
		Vector<T,3> refinedExtend = extend
							- Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX,
											4*coarseDeltaX};
		IndicatorCuboid3D<T> refined(refinedExtend, refinedOrigin);
		coarseGrid.getSuperGeometry().reset(refined);


		if(n>=2) {
		  // Refinement around the wing box - level 2
		  Vector<T,3> fineOrigin2 = wingMin - hn2;
			Vector<T,3> fineExtend2 = wingMax + hp2 - (wingMin - hn2);
			coarseDeltaX = fineGrid.getConverter().getPhysDeltaX();

			//Check if aligned with coarse grid. If not, round
			for(int i=0;i<3;++i) {
				T remainder = std::fmod(fineOrigin2[i],coarseDeltaX);
				if(remainder != 0) {
					fineOrigin2[i] -= remainder;
					cout << "Rounding fineOrigin " << i << " by " << remainder << endl;
					cout << "New remainder = " << std::fmod(fineOrigin2[i],coarseDeltaX) << endl;
				}
				remainder = std::fmod(fineExtend2[i],coarseDeltaX);
				if(remainder != 0) {
					fineExtend2[i] -= remainder;
					cout << "Rounding fineExtend " << i << " by " << remainder << endl;
					cout << "New remainder = " << std::fmod(fineExtend2[i],coarseDeltaX) << endl;
				}
			}

			auto& fineGrid2 = fineGrid.refine(fineOrigin2, fineExtend2,
					false, false, false, false);
			prepareGeometry(fineGrid2, domainOrigin, domainExtend, stlReader, wingLayer, stlReader2);

			origin	= fineGrid2.getOrigin();
      extend	= fineGrid2.getExtend();

	    extendXY  = {extend[0], extend[1], 0};
			extendXZ	= {extend[0], 0, extend[2]};
			extendYZ	= {0, extend[1], extend[2]};

	    fineGrid.addFineCoupling(fineGrid2, origin, extendXY);
			fineGrid.addFineCoupling(fineGrid2, origin, extendXZ);
			fineGrid.addFineCoupling(fineGrid2, origin, extendYZ);

			extendX	= {extend[0], 0, 0};
			extendY	= {0, extend[1], 0};
	    extendZ = {0, 0, extend[2]};

	    fineGrid.addFineCoupling(fineGrid2, origin + extendZ, extendXY);
			fineGrid.addFineCoupling(fineGrid2, origin + extendX, extendYZ);
			fineGrid.addFineCoupling(fineGrid2, origin + extendY, extendXZ);

			innerOrigin = origin
								+ Vector<T,3> {coarseDeltaX, coarseDeltaX, coarseDeltaX};

	    innerExtendXY = {extend[0] - 2*coarseDeltaX,
	                      extend[1] - 2*coarseDeltaX, 0};
			innerExtendXZ = {extend[0] - 2*coarseDeltaX,
											   0, extend[2] - 2*coarseDeltaX};
			innerExtendYZ = {0, extend[1] - 2*coarseDeltaX,
											   extend[2] - 2*coarseDeltaX};

	    fineGrid.addCoarseCoupling(fineGrid2, innerOrigin, innerExtendXY);
			fineGrid.addCoarseCoupling(fineGrid2, innerOrigin, innerExtendXZ);
			fineGrid.addCoarseCoupling(fineGrid2, innerOrigin, innerExtendYZ);

			innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
			innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
	    innerExtendZ = {0, 0, extend[2] - 2*coarseDeltaX};

			fineGrid.addCoarseCoupling(fineGrid2, innerOrigin + innerExtendX, innerExtendYZ);
			fineGrid.addCoarseCoupling(fineGrid2, innerOrigin + innerExtendY, innerExtendXZ);
	    fineGrid.addCoarseCoupling(fineGrid2, innerOrigin + innerExtendZ, innerExtendXY);

			refinedOrigin = origin
								+ Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX,
												2*coarseDeltaX};
			refinedExtend = extend
								- Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX,
												4*coarseDeltaX};
			IndicatorCuboid3D<T> refined2(refinedExtend, refinedOrigin);
			fineGrid.getSuperGeometry().reset(refined2);


			if(n>=3) {
			  // Refinement around the wing box - level 3
			  Vector<T,3> fineOrigin3 = wingMin - hn3;
				Vector<T,3> fineExtend3 = wingMax + hp3 - (wingMin - hn3);
				coarseDeltaX = fineGrid2.getConverter().getPhysDeltaX();

				//Check if aligned with coarse grid. If not, round
				for(int i=0;i<3;++i) {
					T remainder = std::fmod(fineOrigin3[i],coarseDeltaX);
					if(remainder != 0) {
						fineOrigin3[i] -= remainder;
						cout << "Rounding fineOrigin " << i << " by " << remainder << endl;
						cout << "New remainder = " << std::fmod(fineOrigin3[i],coarseDeltaX) << endl;

					}
					remainder = std::fmod(fineExtend3[i],coarseDeltaX);
					if(remainder != 0) {
						fineExtend3[i] -= remainder;
						cout << "Rounding fineExtend " << i << " by " << remainder << endl;
						cout << "New remainder = " << std::fmod(fineExtend3[i],coarseDeltaX) << endl;
					}
				}

				auto& fineGrid3 = fineGrid2.refine(fineOrigin3, fineExtend3,
						false, false, false, false);
				prepareGeometry(fineGrid3, domainOrigin, domainExtend, stlReader, wingLayer, stlReader2);

				origin	= fineGrid3.getOrigin();
				extend	= fineGrid3.getExtend();

		    extendXY  = {extend[0], extend[1], 0};
				extendXZ	= {extend[0], 0, extend[2]};
				extendYZ	= {0, extend[1], extend[2]};

		    fineGrid2.addFineCoupling(fineGrid3, origin, extendXY);
				fineGrid2.addFineCoupling(fineGrid3, origin, extendXZ);
				fineGrid2.addFineCoupling(fineGrid3, origin, extendYZ);

				extendX	= {extend[0], 0, 0};
				extendY	= {0, extend[1], 0};
		    extendZ = {0, 0, extend[2]};

		    fineGrid2.addFineCoupling(fineGrid3, origin + extendZ, extendXY);
				fineGrid2.addFineCoupling(fineGrid3, origin + extendX, extendYZ);
				fineGrid2.addFineCoupling(fineGrid3, origin + extendY, extendXZ);

				innerOrigin = origin
									+ Vector<T,3> {coarseDeltaX, coarseDeltaX, coarseDeltaX};

		    innerExtendXY = {extend[0] - 2*coarseDeltaX,
		                      extend[1] - 2*coarseDeltaX, 0};
				innerExtendXZ = {extend[0] - 2*coarseDeltaX,
												   0, extend[2] - 2*coarseDeltaX};
				innerExtendYZ = {0, extend[1] - 2*coarseDeltaX,
												   extend[2] - 2*coarseDeltaX};

		    fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin, innerExtendXY);
				fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin, innerExtendXZ);
				fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin, innerExtendYZ);

				innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
				innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
		    innerExtendZ = {0, 0, extend[2] - 2*coarseDeltaX};

				fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin + innerExtendX, innerExtendYZ);
				fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin + innerExtendY, innerExtendXZ);
		    fineGrid2.addCoarseCoupling(fineGrid3, innerOrigin + innerExtendZ, innerExtendXY);

				refinedOrigin = origin
									+ Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX,
													2*coarseDeltaX};
				refinedExtend = extend
									- Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX,
													4*coarseDeltaX};
				IndicatorCuboid3D<T> refined3(refinedExtend, refinedOrigin);
				fineGrid2.getSuperGeometry().reset(refined3);


				if(n>=4) {
				  // Refinement around the wing box - level 4 (current innermost)
				  Vector<T,3> fineOrigin4 = wingMin - hn4;
				  Vector<T,3> fineExtend4 = wingMax + hp4 - (wingMin - hn4);
					coarseDeltaX = fineGrid3.getConverter().getPhysDeltaX();

					//Check if aligned with coarse grid. If not, round
					for(int i=0;i<3;++i) {
						T remainder = std::fmod(fineOrigin4[i],coarseDeltaX);
						if(remainder != 0) {
							fineOrigin4[i] -= remainder;
							cout << "Rounding fineOrigin " << i << " by " << remainder << endl;
						}
						remainder = std::fmod(fineExtend4[i],coarseDeltaX);
						if(remainder != 0) {
							fineExtend4[i] -= remainder;
							cout << "Rounding fineExtend " << i << " by " << remainder << endl;
						}
					}

				  auto& fineGrid4 = fineGrid3.refine(fineOrigin4, fineExtend4,
				      false, false, false, false);
				  prepareGeometry(fineGrid4, domainOrigin, domainExtend, stlReader, wingLayer, stlReader2);

			    coarseDeltaX = fineGrid3.getConverter().getPhysDeltaX();
			    origin	= fineGrid4.getOrigin();
			    extend	= fineGrid4.getExtend();

			    extendXY  = {extend[0], extend[1], 0};
			    extendXZ	= {extend[0], 0, extend[2]};
			    extendYZ	= {0, extend[1], extend[2]};

			    fineGrid3.addFineCoupling(fineGrid4, origin, extendXY);
			    fineGrid3.addFineCoupling(fineGrid4, origin, extendXZ);
			    fineGrid3.addFineCoupling(fineGrid4, origin, extendYZ);

			    extendX	= {extend[0], 0, 0};
			    extendY	= {0, extend[1], 0};
			    extendZ = {0, 0, extend[2]};

			    fineGrid3.addFineCoupling(fineGrid4, origin + extendZ, extendXY);
			    fineGrid3.addFineCoupling(fineGrid4, origin + extendX, extendYZ);
			    fineGrid3.addFineCoupling(fineGrid4, origin + extendY, extendXZ);

			    innerOrigin = origin
			              + Vector<T,3> {coarseDeltaX, coarseDeltaX, coarseDeltaX};

			    innerExtendXY = {extend[0] - 2*coarseDeltaX,
			                      extend[1] - 2*coarseDeltaX, 0};
			    innerExtendXZ = {extend[0] - 2*coarseDeltaX,
			                       0, extend[2] - 2*coarseDeltaX};
			    innerExtendYZ = {0, extend[1] - 2*coarseDeltaX,
			                       extend[2] - 2*coarseDeltaX};

			    fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin, innerExtendXY);
			    fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin, innerExtendXZ);
			    fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin, innerExtendYZ);

			    innerExtendX = {extend[0] - 2*coarseDeltaX, 0, 0};
			    innerExtendY = {0, extend[1] - 2*coarseDeltaX, 0};
			    innerExtendZ = {0, 0, extend[2] - 2*coarseDeltaX};

			    fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin + innerExtendX, innerExtendYZ);
			    fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin + innerExtendY, innerExtendXZ);
			    fineGrid3.addCoarseCoupling(fineGrid4, innerOrigin + innerExtendZ, innerExtendXY);

			    refinedOrigin = origin
			              + Vector<T,3> {2*coarseDeltaX, 2*coarseDeltaX,
			                      2*coarseDeltaX};
			    refinedExtend = extend
			              - Vector<T,3> {4*coarseDeltaX, 4*coarseDeltaX,
			                      4*coarseDeltaX};
			    IndicatorCuboid3D<T> refined4(refinedExtend, refinedOrigin);
			    fineGrid3.getSuperGeometry().reset(refined4);
  			}
			}
		}
	}

	clout << "Setup Refinement ... OK" << std::endl;
}

// Create lattice structures
void prepareLattice(Grid3D<T,DESCRIPTOR>& grid, STLreader<T>& stlReader) {

	OstreamManager clout(std::cout, "prepareLattice");
	clout << "Prepare lattice ..." << std::endl;

	auto& converter = grid.getConverter();
	auto& sGeometry = grid.getSuperGeometry();
	auto& sLattice  = grid.getSuperLattice();
	//const T deltaX  = converter.getPhysDeltaX();
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

	sLattice.defineDynamics(sGeometry, 5, &instances::getNoDynamics<T,DESCRIPTOR>());

	// Define boundary conditions
  onbc.addVelocityBoundary(sGeometry, 2, omega);
  //bc.addPressureBoundary(sGeometry, 2, omega);
	onbc.addVelocityBoundary(sGeometry, 3, omega);
	bc.addPressureBoundary(sGeometry, 4, omega);


	if ( bouzidiOn ) {
		// material=5 --> no dynamics + bouzidi zero velocity
		sLattice.defineDynamics( sGeometry,5,&instances::getNoDynamics<T,DESCRIPTOR>() );
		offBc.addZeroVelocityBoundary( sGeometry,5,stlReader );
	} else {
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

// Capture the pressure around the wing, at chosen z-planes --- middle z
void getPressure(Grid3D<T,DESCRIPTOR>& grid, int iT, STLreader<T>& wing) {
	auto& sLattice = grid.getSuperLattice();
	auto& converter = grid.getConverter();

	SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
	AnalyticalFfromSuperF3D<T> interpolatePressure(pressure, true);

	const Vector<T,3> wingMin = wing.getMesh().getMin();
	const Vector<T,3> wingMax = wing.getMesh().getMax();

	const T wingSpan = 0.040;
	const T wingChord = 0.04523;
	const T wingHeight = 0.003137;
	const T theta = 0.; //Pitch
	const std::vector<T> cR {0.5*wingChord,0.5*wingHeight};

	//Z-position and name-tag for each x-y plane to read pressure
	const std::vector<T> z {0.7*wingSpan};
	const std::vector<std::string> tags {"0.5span_"};

	for(std::vector<std::string>::size_type i=0; i<tags.size(); ++i) {

		ofstream myfile;
		std::string filename {"tmp/pressure_" + tags[i] +  to_string(iT) + ".csv"};
		myfile.open(filename, fstream::app);

		//X-coordinates of current surface points (relative to local wing origin)
	  std::vector<T> surface_x{0, 0.0000805780000000027, 0.000283419999999999,
			 0.0029097721, 0.0107896067, 0.018852055, 0.0269683804, 0.0350313102,
			 0.0429601641, 0.04501812925, 0.04521576125, 0.04506297569,
			 0.0439259505, 0.0366446105, 0.0285898428, 0.020474961,
			 0.0123952441, 0.0043404288, 0.000492062400000001,
			 0.0001141131, 0.0000025481000000056, 0.009046147, 0.013569221,
			 0.018092295, 0.022615368, 0.009174991, 0.013998084, 0.022910189};

		std::vector<T> surface_y{0.0015917555461, 0.00178170221399999,
			 0.00188643266799999, 0.00216398364599999, 0.00278718292,
			 0.00310168391999999, 0.00308969873, 0.00275129778, 0.002099574308,
			 0.00186987491699999, 0.00171069978199999, 0.001471687579,
			 0.00128696748099999, 0.00113842089, 0.000492467779999999,
			 0.0000890504999999999, 0.0000112311699999998, 0.00026105833,
			 0.000836426680999999, 0.001227102473, 0.001322748836, 0.0015219141547,
			 0.002676053, 0.002932016, 0.003085846, 0.003136704, 0.002684994,
			 0.002951373, 0.003136867};

		 //Option to rotate above points by a specified pitch
		 for(std::vector<T>::size_type j = 0; j < surface_x.size(); ++j) {
			 T xMxrC = surface_x[j] - cR[0];
			 T yMyrC = surface_y[j] - cR[1];

			 surface_x[j] = xMxrC * std::cos(3.141592653589793238463 * theta/180)
			  	- yMyrC * std::sin(3.141592653589793238463 * theta/180) + cR[0];

			 surface_y[j] = yMyrC * std::cos(3.141592653589793238463 * theta/180)
					+ yMyrC * std::cos(3.141592653589793238463 * theta/180) + cR[1];
		 }

		for(std::vector<T>::size_type j = 0; j < surface_x.size(); ++j) {
			const T point[3] {wingMin[0] + surface_x[j], wingMin[1] + surface_y[j],
				 		wingMin[2] + z[i]};
			T pressureAtPoint {};
			interpolatePressure(&pressureAtPoint, point);

			myfile << surface_x[j] << "	" << pressureAtPoint << std::endl;
		}
	  myfile.close();
 	}
}

void getConvergenceStats(Grid3D<T,DESCRIPTOR>& coarseGrid, int iT,
	 util::ValueTracer<T>& tracer) {
	T rhoAv = coarseGrid.getSuperLattice().getStatistics().getAverageRho();

	//Later - upgrade to take average across all grids, not just coarse
	tracer.takeValue(rhoAv, true);

	//Write to output at chosen interval
	if(iT % 100 ==0){
		ofstream myfile;
		std::string filename {"tmp/rhoAv.csv"};
		myfile.open(filename, fstream::app);
		myfile << iT << "	" << rhoAv << std::endl;
		myfile.close();
	}

	if ((iT > 100) && tracer.hasConverged()) {
		std::cout << "Simulation converged." << std::endl;
	}
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
      LatticeVelocity<T>(uC),
      N,
      PhysCharacteristics);

  //Overall domain dimensions
  const Vector<T,3> domainOrigin =
    coarseGrid.getSuperGeometry().getStatistics().getMinPhysR(0);
  const Vector<T,3> domainExtend =
    coarseGrid.getSuperGeometry().getStatistics().getPhysExtend(0);

  // === 2nd Step: Prepare Geometry ===

  //Instantiate stlReader and Layer indicators for wing
  const T finestPhysDx = coarseGrid.getConverter().getPhysDeltaX()*std::pow(0.5,nRefinement); //Assuming 4 refinement levels
  STLreader<T> stlReader( "lyonWingBZ_domain2_0deg_standard.stl", finestPhysDx, 0.001, 1, true );
	//Second STL reader to allow setting interior of wing to zero
	STLreader<T> stlReader2( "lyonWing_domain2_0deg_standard.stl", finestPhysDx, 0.001, 1, true );

	IndicatorLayer3D<T> wingLayer( stlReader, finestPhysDx );

  prepareGeometry(coarseGrid, domainOrigin, domainExtend, stlReader, wingLayer, stlReader2);

  setupRefinement(coarseGrid, domainOrigin, domainExtend, stlReader, wingLayer, stlReader2, nRefinement);

  // === 3rd Step: Prepare Lattice ===
  coarseGrid.forEachGrid(prepareLattice,stlReader);

  clout << "Total number of active cells: " << coarseGrid.getActiveVoxelN() << std::endl;

	//Reference to finest grid containing wing
	Grid3D<T,DESCRIPTOR>& wingGrid = coarseGrid.locate(stlReader.getMesh().getMin());
  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  Timer<T> timer(
			coarseGrid.getConverter().getLatticeTime(maxPhysT),
			coarseGrid.getSuperGeometry().getStatistics().getNvoxel() );
  timer.start();

	// Convergence tracer every physical second
	util::ValueTracer<T> rhoTracer(/*coarseGrid.getConverter().getLatticeTime(0.01)*/1000, 1e-5);

  const int vtkIter   	 = 1000; //Every 10% of max physical time
  const int statIter  	 = 100;
	const int pressureIter = 1000;
	const int checkIter 	 = 1000;

  for ( int iT = 0; iT < coarseGrid.getConverter().getLatticeTime( maxPhysT ); ++iT ) {

		// Load last checkpoint at startup
		if (iT == 0 && loadCheckpoint) {
			coarseGrid.forEachGrid("cylinder3D_even", [&](Grid3D<T,DESCRIPTOR>& grid,
						const std::string& id)
					{grid.getSuperLattice().load(id+".checkpoint");
					 clout << "Checkpoint loaded." << std::endl;});
		}

    // === 5th Step: Definition of Initial and Boundary Conditions ===
		if (startUpOn) {
    	setBoundaryValues( coarseGrid, iT);
		}

    // === 6th Step: Collide and Stream Execution ===
    coarseGrid.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===

    if ( iT % vtkIter == 0 ) {
      coarseGrid.forEachGrid("aerofoil3D", [&](Grid3D<T,DESCRIPTOR>& grid,
          const std::string& id)
        {getVTK(grid, id, iT); });
      clout << "Get results vtk" << endl;
    }

		// Save checkpoint
		if ( (iT % checkIter == 0) && (iT != 0) && saveCheckpoint) {
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

    if ( iT % statIter == 0 ) {
      getStats(coarseGrid, iT, timer);
      clout << "Get results stats" << endl;
    }

		if ( iT % pressureIter == 0 ) {
			getPressure(wingGrid,iT,stlReader2);
			clout << "Get pressure" << endl;
		}

		//Residuals and convergence checks
		getConvergenceStats(coarseGrid, iT, rhoTracer);


  }

  timer.stop();
  timer.printSummary();
}