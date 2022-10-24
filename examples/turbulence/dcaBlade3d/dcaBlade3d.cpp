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
//#define Bouzidi 
#define Grad

#ifdef WALE
#define DESCRIPTOR WALED3Q19Descriptor
#elif defined (Smagorinsky)
#define DESCRIPTOR D3Q19<>
#elif defined (KBC)
#ifdef Grad
#ifdef sponge
#define DESCRIPTOR D3Q27descriptorKBCGradSponge
#else
#define DESCRIPTOR D3Q27descriptorKBCGrad
#endif
#elif defined(sponge)
#define DESCRIPTOR D3Q27descriptorKBCSponge
#else
#define DESCRIPTOR D3Q27descriptorKBC
#endif
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
  const T span    = indicatorBlade.getSpan();

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
                                            bladeOrigin[2]  
                                              + 0.5 * span - 0.5 * deltaX};
  const Vector<T,3> pressureSection1Extend {2. * chord, 2. * chord, deltaX};
  IndicatorCuboid3D<T> pressureSection1(pressureSection1Extend,
                                        pressureSection1Origin);
  sGeometry.rename(5, 7, pressureSection1);

  //Front face
  {
    const Vector<T,3> indiOrigin {origin[0] - deltaX / 2.,
                                  origin[1] - deltaX / 2.,
                                  origin[2] - 50 * deltaX / 2.};
    const Vector<T,3> indiExtend {deltaX,
                                  extend[1] + deltaX,
                                  extend[2] + 100 * deltaX};
    IndicatorCuboid3D<T> ff(indiExtend, indiOrigin);
    sGeometry.rename(1, 3, ff);
  }

  //Upper face
  {
    const Vector<T,3> indiOrigin {origin[0] + deltaX / 2,
                                  origin[1] + extend[1] - deltaX / 2.,
                                  origin[2] - 50 * deltaX / 2.};
    const Vector<T,3> indiExtend {extend[0],
                                  deltaX,
                                  extend[2] + 100 * deltaX};
    IndicatorCuboid3D<T> uf(indiExtend, indiOrigin);
    sGeometry.rename(1, 4, uf);
  }

  //Lower face
  {
    const Vector<T,3> indiOrigin {origin[0] + deltaX / 2,
                                  origin[1] - deltaX / 2.,
                                  origin[2] - 50 * deltaX / 2.};
    const Vector<T,3> indiExtend {extend[0],
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
                                  extend[1] + deltaX / 2.,
                                  extend[2] + 100 * deltaX};
    IndicatorCuboid3D<T> rf(indiExtend, indiOrigin);
    sGeometry.rename(1, 4, rf);
  }

  //Upper rear edge
  {
    const Vector<T,3> indiOrigin {origin[0] + extend[0] - deltaX / 2.,
                                  origin[1] + extend[1] - deltaX / 2.,
                                  origin[2] - 50 * deltaX / 2.};
    const Vector<T,3> indiExtend {deltaX,
                                  deltaX,
                                  extend[2] + 100 * deltaX};
    IndicatorCuboid3D<T> ure(indiExtend, indiOrigin);
    sGeometry.rename(4, 3, ure);
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

  const Vector<T,2> hn6 = {0.018 * chord, 0.018 * chord}; //unused in baseline!
  const Vector<T,2> hp6 = {0.018 * chord, 0.018 * chord}; // '' positive

  const Vector<T,2> hn5 = {0.054 * chord, 0.054 * chord}; 
  const Vector<T,2> hp5 = {0.054 * chord, 0.054 * chord}; // '' positive

  const Vector<T,2> hn4 = {0.126 * chord, 0.126 * chord};
  const Vector<T,2> hp4 = {0.126 * chord, 0.126 * chord};

  const Vector<T,2> hn3 = {0.27 * chord, 1.0 * chord};
  const Vector<T,2> hp3 = {3.0 * chord, 1.0 * chord};

  const Vector<T,2> hn2 = {0.558 * chord, 2.0 * chord}; //Outermost
  const Vector<T,2> hp2 = {14.0 * chord, 2.0 * chord};

  const Vector<T,2> hn1 = {1.134 * chord, 4.0 * chord}; 
  const Vector<T,2> hp1 = {28.0 * chord, 4.0 * chord}; // '' positive

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

          if(n >= 5) {
	          // Refinement around the wing box - level 5 (current innermost)
            const T deltaX3 = fineGrid4.getConverter().getPhysDeltaX();
	            coarseDeltaX = deltaX3;
			
	          Vector<T,3> fineOrigin5 = {bladeBoxOrigin[0] - hn5[0],
		                       bladeBoxOrigin[1] - hn5[1],
				                   domainOrigin[2]-deltaX0-deltaX1-deltaX2-deltaX3};
	          Vector<T,3> fineExtend5 = 
	          {chord + hp5[0] + hn5[0], thickness + hp5[1] + hn5[1],
                 domainExtend[2] + deltaX0 + deltaX1 + deltaX2 + deltaX3};

        	  auto& fineGrid5 = fineGrid4.refine(fineOrigin5, fineExtend5,
		  		             false, false, true, false);
	          prepareGeometry(fineGrid5, domainOrigin, domainExtend,
			          indicatorBlade);
	          origin = fineGrid5.getOrigin() + 
                 Vector<T,3>{0., 0., 0.5 * coarseDeltaX};
      	    extend = fineGrid5.getExtend() - 
              Vector<T,3>{0., 0., 0.5 * coarseDeltaX};
	          extendXZ = {extend[0], 0., extend[2]};
	          extendYZ = {0., extend[1], extend[2]};
            fineGrid4.addFineCoupling(fineGrid5, origin, extendXZ);
	          fineGrid4.addFineCoupling(fineGrid5, origin, extendYZ);
            extendX = {extend[0], 0., 0.};
	          extendY = {0., extend[1], 0.};
            fineGrid4.addFineCoupling(fineGrid5, origin + extendX, extendYZ);
  	        fineGrid4.addFineCoupling(fineGrid5, origin + extendY, extendXZ);

        	  innerOrigin = origin + Vector<T,3> {coarseDeltaX, coarseDeltaX, 0.};
            innerExtendXZ = {extend[0] - 2. * coarseDeltaX , 0., extend[2]};
	          innerExtendYZ = {0., extend[1] - 2. * coarseDeltaX, extend[2]};
	          fineGrid4.addCoarseCoupling(fineGrid5, innerOrigin, innerExtendXZ);
            fineGrid4.addCoarseCoupling(fineGrid5, innerOrigin, innerExtendYZ);

        	  innerExtendX = {extend[0] - 2. * coarseDeltaX, 0., 0.};
	          innerExtendY = {0., extend[1] - 2. * coarseDeltaX, 0.};
	          fineGrid4.addCoarseCoupling(fineGrid5,
			                innerOrigin + innerExtendX,
				        innerExtendYZ);
	          fineGrid4.addCoarseCoupling(fineGrid5,
			                innerOrigin + innerExtendY,
				        innerExtendXZ);

        	  refinedOrigin = origin + Vector<T,3> {2. * coarseDeltaX,
		                                  2. * coarseDeltaX,
				  	       	- 2. * coarseDeltaX};
	          refinedExtend = extend - Vector<T,3> {4. * coarseDeltaX,
		                                  4. * coarseDeltaX,
			                        	- 4. * coarseDeltaX};
	          IndicatorCuboid3D<T> refined5(refinedExtend, refinedOrigin);
	          fineGrid4.getSuperGeometry().reset(refined5);

            if(n >= 6) {
	            // Refinement around the wing box - level 5 (current innermost)
              const T deltaX4 = fineGrid5.getConverter().getPhysDeltaX();
	              coarseDeltaX = deltaX4;
			
	            Vector<T,3> fineOrigin6 = {bladeBoxOrigin[0] - hn6[0],
		                         bladeBoxOrigin[1] - hn6[1],
				                     domainOrigin[2]-deltaX0-deltaX1-deltaX2-deltaX3-deltaX4};
	            Vector<T,3> fineExtend6 = 
	            {chord + hp6[0] + hn6[0], thickness + hp6[1] + hn6[1],
                   domainExtend[2] + deltaX0 + deltaX1 + deltaX2 + deltaX3 + deltaX4};

          	  auto& fineGrid6 = fineGrid5.refine(fineOrigin6, fineExtend6,
		    		             false, false, true, false);
	            prepareGeometry(fineGrid6, domainOrigin, domainExtend,
			            indicatorBlade);
	            origin = fineGrid6.getOrigin() + 
                   Vector<T,3>{0., 0., 0.5 * coarseDeltaX};
      	      extend = fineGrid6.getExtend() - 
                Vector<T,3>{0., 0., 0.5 * coarseDeltaX};
	            extendXZ = {extend[0], 0., extend[2]};
	            extendYZ = {0., extend[1], extend[2]};
              fineGrid5.addFineCoupling(fineGrid6, origin, extendXZ);
	            fineGrid5.addFineCoupling(fineGrid6, origin, extendYZ);
              extendX = {extend[0], 0., 0.};
	            extendY = {0., extend[1], 0.};
              fineGrid5.addFineCoupling(fineGrid6, origin + extendX, extendYZ);
  	          fineGrid5.addFineCoupling(fineGrid6, origin + extendY, extendXZ);

          	  innerOrigin = origin + Vector<T,3> {coarseDeltaX, coarseDeltaX, 0.};
              innerExtendXZ = {extend[0] - 2. * coarseDeltaX , 0., extend[2]};
	            innerExtendYZ = {0., extend[1] - 2. * coarseDeltaX, extend[2]};
	            fineGrid5.addCoarseCoupling(fineGrid6, innerOrigin, innerExtendXZ);
              fineGrid5.addCoarseCoupling(fineGrid6, innerOrigin, innerExtendYZ);

          	  innerExtendX = {extend[0] - 2. * coarseDeltaX, 0., 0.};
	            innerExtendY = {0., extend[1] - 2. * coarseDeltaX, 0.};
	            fineGrid5.addCoarseCoupling(fineGrid6,
			                  innerOrigin + innerExtendX,
				          innerExtendYZ);
	            fineGrid5.addCoarseCoupling(fineGrid6,
			                  innerOrigin + innerExtendY,
				          innerExtendXZ);

          	  refinedOrigin = origin + Vector<T,3> {2. * coarseDeltaX,
		                                    2. * coarseDeltaX,
				    	       	- 2. * coarseDeltaX};
	            refinedExtend = extend - Vector<T,3> {4. * coarseDeltaX,
		                                    4. * coarseDeltaX,
			                          	- 4. * coarseDeltaX};
	            IndicatorCuboid3D<T> refined6(refinedExtend, refinedOrigin);
	            fineGrid5.getSuperGeometry().reset(refined6);
            }
          }
  	    }
      }
    }
  }
  clout << "Setup Refinement ... OK" << std::endl;
}

// Create lattice structures
void prepareLattice(Grid3D<T,DESCRIPTOR>& grid,
		    IndicatorBladeDca3D<T>& indicatorBlade,
		    const bool& bouzidiOn,
        const T& thetaBC) {
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
      #if defined(sponge)
        Dynamics<T,DESCRIPTOR>& bulkDynamics = 
        grid.addDynamics(std::unique_ptr<Dynamics<T,DESCRIPTOR>>(
          new KBCGradSpongeDynamics<T,DESCRIPTOR>(
            omega, instances::getKBCBulkMomenta<T,DESCRIPTOR>())));
      #else
      Dynamics<T,DESCRIPTOR>& bulkDynamics = 
        grid.addDynamics(std::unique_ptr<Dynamics<T,DESCRIPTOR>>(
          new KBCGradDynamics<T,DESCRIPTOR>(
            omega, instances::getKBCBulkMomenta<T,DESCRIPTOR>())));
      #endif
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

  auto bulkIndicator = sGeometry.getMaterialIndicator({1, 2, 3, 4, 6});
  sLattice.defineDynamics(bulkIndicator, &bulkDynamics);

  // Define boundary conditions
  //onbc.addSlipBoundary(sGeometry, 2);// SLIP BC
  //bc.addOutletBoundary(sGeometry, 4, {1, 2, 3, 4, 6});
  //OLD BC:
  //bc.addPressureBoundary(sGeometry, 4, omega);
  //bc.addVelocityBoundary(sGeometry, 3, omega);

  //NEXT RECONFIGURE BCS TO INLET OUTLET


  #if defined(Bouzidi)
    // material=5, 7 --> no dynamics + bouzidi zero velocity
    sLattice.defineDynamics( sGeometry,5,&instances::getNoDynamics<T,DESCRIPTOR>() );
    sLattice.defineDynamics( sGeometry,7,&instances::getNoDynamics<T,DESCRIPTOR>() );
    offBc.addSecondOrderZeroVelocityBoundary( sGeometry,5,indicatorBlade );
    offBc.addSecondOrderZeroVelocityBoundary( sGeometry,7,indicatorBlade );
  #elif defined(Grad)
    sLattice.defineDynamics( sGeometry,5,&instances::getNoDynamics<T,DESCRIPTOR>() );
    sLattice.defineDynamics( sGeometry,7,&instances::getNoDynamics<T,DESCRIPTOR>() );
    offBc.addZeroVelocityGradBoundary( sGeometry,5,indicatorBlade,std::vector<int>{1} );
    offBc.addZeroVelocityGradBoundary( sGeometry,7,indicatorBlade,std::vector<int>{1} );
  #else
    //material=5,7 --> fullway bounceBack dynamics
    sLattice.defineDynamics( sGeometry, 5, &instances::getBounceBack<T, DESCRIPTOR>() );
    sLattice.defineDynamics( sGeometry, 7, &instances::getBounceBack<T, DESCRIPTOR>() );
  #endif

  #ifdef sponge
    //Define and initialise viscosity sponge zones
    //Sponge indicator 1 - y-z outlet
    const T physChord = 0.051;
    const T deltaX = converter.getPhysDeltaX();
    const Vector<T,3> spongeOrigin = {56. * physChord - deltaX /2000, - deltaX / 2,
      - 4 * deltaX};
    const Vector<T,3> spongeExtend = {4. * physChord + deltaX/1000,
      40 * physChord + deltaX,
      0.2 * physChord + 8 * deltaX};
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

    //Sponge indicator 2 - x-z 
    const Vector<T,3> spongeOrigin2 = {0. * physChord - deltaX /2, 36. * physChord - deltaX / 2000,
      - 4 * deltaX};
    const Vector<T,3> spongeExtend2 = {60. * physChord + deltaX,
      4. * physChord + deltaX / 1000,
      0.2 * physChord + 8 * deltaX};
    IndicatorCuboid3D<T> spongeRegion2(spongeExtend2, spongeOrigin2);
    //Orientation
    const Vector<T,3> spongeOrientation2 = {0., 1., 0.};

    outletSponge.addSineSponge(sGeometry, spongeRegion2, spongeOrientation2,
      tauSpongeBase, tauSpongeMax, spongeMaterials);

    sLattice.initialiseSponges();
  #endif

  // Initial conditions - characteristic physical velocity and density for inflow
  AnalyticalConst3D<T,T> rhoF {1.};
  Vector<T,3> velocityV {converter.getCharLatticeVelocity() * cos(thetaBC * M_PI / 180.),
                         converter.getCharLatticeVelocity() * sin(thetaBC * M_PI / 180.), 
                         0.};
  //Rotation tranformation here
  AnalyticalConst3D<T,T> uF(velocityV);

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(Grid3D<T,DESCRIPTOR>& grid, int iT, const T& thetaBC) {

	OstreamManager clout(std::cout, "setBoundaryValues");

	auto& converter	= grid.getConverter();
	auto& sGeometry	= grid.getSuperGeometry();
	auto& sLattice	= grid.getSuperLattice();

  Vector<T,3> inVel {
    converter.getCharLatticeVelocity() * cos(thetaBC * M_PI / 180.),
    converter.getCharLatticeVelocity() * sin(thetaBC * M_PI / 180.),
    0.};
  T inRho = 1.0;
    
  AnalyticalConst3D<T,T> inRhoConst(inRho);
  AnalyticalConst3D<T,T> inVelConst(inVel);

  //sLattice.defineRhoU(sGeometry, 3, inRhoConst, inVelConst);
  sLattice.iniEquilibrium(sGeometry, 3, inRhoConst, inVelConst);
  sLattice.iniEquilibrium(sGeometry, 4, inRhoConst, inVelConst);
}

// Output results to vtk files
void getVTK(Grid3D<T,DESCRIPTOR>& grid, const std::string& prefix, int iT,
	    SuperLatticePhysWallShearStressAndPressure3D<T,DESCRIPTOR>& wssp,
      SuperLatticeYplus3D<T,DESCRIPTOR>& sYplus
	    /*SuperLatticeTimeAveragedF3D<T>& sAveragedWSSP*/) {

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
  vtmWriter.addFunctor(sYplus);
  //vtmWriter.addFunctor(sAveragedWSSP);

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
		                                 std::vector<int>{5,7},
						 converter,0.5*span*chord); //physical span is half geometric!
  T bladeForce[3];
  int input1[0];
  drag(bladeForce,input1);
  //std::cout << "Lift = " << bladeForce[1] << endl;
  //std::cout << "Drag = " << bladeForce[0] << endl;

  ofstream myfile;
  std::string filename {filePath+".csv"};
  myfile.open(filename,fstream::app);
  myfile << theta << "	" << bladeForce[1] << "	" << bladeForce[0] << std::endl;
  myfile.close();
}	

int main( int argc, char* argv[] ) {

  //Blade physical parameters
  const T chord = 0.051; //m
  const T thickness = 0.00382;
  const T span = 0.2 * chord; //Twice as wide, to ensure entire domain is spanned
  const T r1 = 0.1836;  //Upper/lower radius
  const T r2 = 0.00015; //LE/TE radius
  const T xp = 0.02538; //Intersect point
  const T theta = 0.00; //Pitch (+ve = anticlockwise)
  const T thetaBC = 8.00; //Inlet flow angle
  const Vector<T,3> bladeOrigin = {10.5 * chord + chord / 12800., 20. * chord + chord / 12800., - 0.5 * span}; //Origin of blade (make sure it's off-node!)

  //Domain and simulation parameters
  const int N = 30; //14        // resolution of the model (coarse cells per chord)
  const int nRefinement = 6;	//Number of refinement levels (current max = 5)
  const T lDomainPhysx = 60.*chord; //Length of domain in physical units (m)
  const T lDomainPhysy = 40.*chord;
  const T lDomainPhysz = 0.2*chord; //
  const T physL = chord; //Physical reference length (m)

  //Flow conditions
  const T Re = 100000.;       // Reynolds number
  const T Mach = 0.1;
  const T uC = Mach * 1./std::pow(3,0.5); //Lattice characteristic velocity
  const T physuC = Mach * 343.; //Physical characteristic velocity
  const T rho = 1.2;	//Density
  const T physNu = physuC * physL / Re;//m2/s
  const T normFactor = physL / physuC;  //Factor for physical -> convective time (s)

  const T maxNormT = 100; //Max normalised 'convective' time
  const T maxPhysT = maxNormT * normFactor; // max. simulation time in s, SI unit

  //Options for blade surface boundary condition
  const bool bouzidiOn = true; //true = bouzidi, false = fullway bb

  //Time-loop options
  const int vtkIter	   = 500; 
  const int statIter       = 100; 
  const int checkIter 	   = 500;
  const int bladeForceIter = 1;
  const int timeAvgIter    = 500;

  //Checkpoint option
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
  IndicatorBladeDca3D<T> blade(bladeOrigin,chord, thickness, 2. * span, r1, r2, xp,
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
    IndicatorBladeDca3D<T>&, const bool&, const T&)>(prepareLattice),blade,bouzidiOn,thetaBC); 
  clout << "Total number of active cells: " << coarseGrid.getActiveVoxelN() 
	      << std::endl;

  //Reference to finest grid containing wing
  Grid3D<T,DESCRIPTOR>& wingGrid = coarseGrid.locate(
    Vector<T,3>(bladeOrigin[0],bladeOrigin[1],bladeOrigin[2]+span/2.));

	//Functor vectors for 3D VTK
	std::vector<std::unique_ptr<SuperLatticePhysVelocity3D<T,DESCRIPTOR>>> sVel;
	std::vector<std::unique_ptr<SuperLatticePhysPressure3D<T,DESCRIPTOR>>> sP;
	std::vector<std::unique_ptr<SuperLatticeYplus3D<T,DESCRIPTOR>>> sYplus;
	std::vector<std::unique_ptr<SuperLatticePhysWallShearStressAndPressure3D<
    T,DESCRIPTOR>>> wssp;
	//std::vector<std::unique_ptr<SuperLatticeTimeAveragedF3D<T>>> sAveragedWSSP;


  //Helper lambdas for initialisation and management of data
  auto initialiseVTK = [](Grid3D<T,DESCRIPTOR>& grid, auto& sVel, auto& sP, auto& sYplus,
	  auto& wssp,/* auto& sAveragedWSSP,*/
	  IndicatorBladeDca3D<T>& indicatorBlade) {
			auto& sGeometry = grid.getSuperGeometry();
			auto& sLattice = grid.getSuperLattice();
			auto& converter = grid.getConverter();

      typedef typename std::remove_reference<decltype(sVel)>::
        type::value_type sVelType; 
      typedef typename std::remove_reference<decltype(sP)>::
        type::value_type sPtype;
      typedef typename std::remove_reference<decltype(sYplus)>::
        type::value_type sYplusType;
      typedef typename std::remove_reference<decltype(wssp)>::
        type::value_type wsspType;
      //typedef typename std::remove_reference<decltype(sAveragedWSSP)>::
      //  type::value_type taType;

			sVel.push_back(sVelType(new typename sVelType::element_type(
        sLattice, converter)));
			sP.push_back(sPtype(new typename sPtype::element_type(
        sLattice, converter)));
      sYplus.push_back(sYplusType(new typename sYplusType::element_type(
        sLattice, converter, sGeometry, indicatorBlade, 7)));
			wssp.push_back(wsspType(new typename wsspType::element_type(
        sLattice,converter, sGeometry,7,indicatorBlade)));
			//sAveragedWSSP.push_back(taType(new typename taType::element_type(
      //  *wssp.back())));
	};

  auto loadCheckpoint = [](Grid3D<T,DESCRIPTOR>& grid, std::string&& id,
    const std::string& checkpoint, int& i_grid) {
      grid.getSuperLattice().load(id+".checkpoint");
			std::cout << checkpoint + " checkpoint loaded." << std::endl;
      i_grid++;
      id = "dcaBlade3d_"+std::to_string(i_grid)+"_"+checkpoint;
  };

  //auto addTaEnsemble = [](Grid3D<T,DESCRIPTOR>& grid,
	//			auto& sAveragedWSSP, int& i_grid) {
	//				auto& sLattice = grid.getSuperLattice();
	//	      sLattice.communicate();
	//				sAveragedWSSP[i_grid]->addEnsemble();
	//				i_grid++;
	//};
 
  auto writeTaVTK = [](Grid3D<T,DESCRIPTOR>& grid, std::string&& id, int& iT,
	  auto& wsspVector, auto& sYplusVector,/* auto& sAveragedWSSPVector,*/ int& i_grid){
			auto& wssp= *wsspVector[i_grid];
			auto& sYplus= *sYplusVector[i_grid];
			//auto& sAveragedWSSP = *sAveragedWSSPVector[i_grid];
			getVTK(grid, id, iT, wssp, sYplus/*,
        sAveragedWSSP*/);
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
  coarseGrid.forEachGrid(initialiseVTK, sVel, sP, sYplus,
    wssp/*, sAveragedWSSP*/, blade);

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

    setBoundaryValues(coarseGrid,iT,thetaBC);

    // === 7th Step: Computation and Output of the Results ===
		//Add ensemble to time-averaged functors
		//if ( iT % timeAvgIter == 0) {
		//	i_grid = 0;
		//	coarseGrid.forEachGrid(addTaEnsemble,
    //   sAveragedWSSP, i_grid);
		//}

		//Add time-averaged functor vector into getVTK functions
    if ( iT % vtkIter == 0 ) {
			i_grid = 0;	
      coarseGrid.forEachGrid(writeTaVTK,"dcaBlade3d_"+std::to_string(i_grid),
        iT, wssp, sYplus/*, sAveragedWSSP*/, i_grid);
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
