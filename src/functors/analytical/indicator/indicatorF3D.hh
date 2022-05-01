/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Mathias J. Krause, Albert Mink
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

#ifndef INDICATOR_F_3D_HH
#define INDICATOR_F_3D_HH

#include <vector>
#include <cmath>
#include <cassert>
#include <sstream>
#include <algorithm> //std::min

#include "indicatorF3D.h"
#include "indicCalc3D.h"
#include "utilities/vectorHelpers.h"

using namespace std;

namespace olb {

template <typename S>
IndicatorTranslate3D<S>::IndicatorTranslate3D(std::array<S,3> translate, IndicatorF3D<S>& indicator)
  : _translate(translate), _indicator(indicator)
{
}

template< typename S>
bool IndicatorTranslate3D<S>::operator() (bool output[], const S input[] )
{
  const S inputTranslated[3] = {
    input[0] - _translate[0],
    input[1] - _translate[1],
    input[2] - _translate[2] };

  _indicator.operator()(output, inputTranslated);
  return output[0];
}

template <typename S>
IndicatorCircle3D<S>::IndicatorCircle3D(Vector<S,3> center, Vector<S,3> normal, S radius)
  :  _center(center), _normal(normal), _radius2(radius*radius)
{
  _normal.normalize();
  this->_myMin = _center - radius;
  this->_myMax = _center + radius;
}

template <typename S>
IndicatorCircle3D<S>::IndicatorCircle3D(S center0, S center1, S center2,
                                        S normal0, S normal1, S normal2, S radius)
  :  _radius2(radius*radius)
{
  _center = {center0, center1, center2};
  _normal = {normal0, normal1, normal2};
  _normal.normalize();
  this->_myMin = _center - radius;
  this->_myMax = _center + radius;
}

// returns true if x is inside the circle
template <typename S>
bool IndicatorCircle3D<S>::operator()(bool output[], const S input[])
{
  S eps = std::numeric_limits<S>::epsilon();
  IndicatorCylinder3D<S> cylinder(_center, _normal, getRadius(), 5*eps);
  return cylinder(output,input);
}

template <typename S>
Vector<S,3> const& IndicatorCircle3D<S>::getCenter() const
{
  return _center;
}

template <typename S>
Vector<S,3> const& IndicatorCircle3D<S>::getNormal() const
{
  return _normal;
}

template <typename S>
S IndicatorCircle3D<S>::getRadius() const
{
  return sqrt(_radius2);
}



template <typename S>
IndicatorSphere3D<S>::IndicatorSphere3D(Vector<S,3> center, S radius)
  :  _center(center), _radius2(radius*radius)
{
  this->_myMin = _center - radius;
  this->_myMax = _center + radius;
}

template <typename S>
IndicatorSphere3D<S>::IndicatorSphere3D(const IndicatorSphere3D& sphere)
{
  this->_myMin = sphere._myMin;
  this->_myMax = sphere._myMax;
  _center = sphere._center;
  _radius2 = sphere._radius2;
}

// returns true if x is inside the sphere
template <typename S>
bool IndicatorSphere3D<S>::operator()(bool output[], const S input[])
{
  output[0] = (  (_center[0] - input[0]) * (_center[0]-input[0])
                 +(_center[1] - input[1]) * (_center[1]-input[1])
                 +(_center[2] - input[2]) * (_center[2]-input[2]) <= _radius2 );
  return true;
}

template <typename S>
bool IndicatorSphere3D<S>::distance(S& distance, const Vector<S,3>& origin)
{
  distance = sqrt(  (_center[0] - origin[0]) * (_center[0]-origin[0])
                 +(_center[1] - origin[1]) * (_center[1]-origin[1])
                 +(_center[2] - origin[2]) * (_center[2]-origin[2]) ) - sqrt(_radius2);
  return true;
}

template <typename S>
bool IndicatorSphere3D<S>::distance(S& distance, const Vector<S,3>& origin,
                                    const Vector<S,3>& direction, int iC)
{
  // computes pos. distance by solving quadratic equation by a-b-c-formula
  S a = direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2];

  // returns 0 if point is at the boundary of the sphere
  if ( util::approxEqual(a,_radius2) ) {
    distance = S();
    return true;
  }
  // norm of direction
  a = sqrt(a);

  S b = 2.*((origin[0] - _center[0])*direction[0] +
            (origin[1] - _center[1])*direction[1] +
            (origin[2] - _center[2])*direction[2])/a;
  S c = -_radius2 + (origin[0] - _center[0])*(origin[0] - _center[0])
        + (origin[1] - _center[1])*(origin[1] - _center[1])
        + (origin[2] - _center[2])*(origin[2] - _center[2]);

  // discriminant
  S d = b*b - 4.*c;
  if (d < 0) {
    return false;
  }

  S x1 = (- b + sqrt(d)) *0.5;
  S x2 = (- b - sqrt(d)) *0.5;

  // case if origin is inside the sphere
  if ((x1<0.) || (x2<0.)) {
    if (x1>0.) {
      distance = x1;
      return true;
    }
    if (x2>0.) {
      distance = x2;
      return true;
    }
  }
  // case if origin is ouside the sphere
  else {
    distance = min(x1,x2);
    return true;
  }

  return false;
}


template <typename S>
IndicatorLayer3D<S>::IndicatorLayer3D(IndicatorF3D<S>& indicatorF, S layerSize)
  :  _indicatorF(indicatorF), _layerSize(layerSize)
{
  this->_myMin = indicatorF.getMin() - layerSize;
  this->_myMax = indicatorF.getMax() + layerSize;
}

// returns true if x is inside the layer
template <typename S>
bool IndicatorLayer3D<S>::operator()(bool output[], const S input[])
{
  output[0] = false;
  S r[3];
  for (int iX =- 1; iX < 2; ++iX) {
    for (int iY =- 1; iY < 2; ++iY) {
      for (int iZ =- 1; iZ < 2; ++iZ) {
        r[0] = input[0] + iX*_layerSize;
        r[1] = input[1] + iY*_layerSize;
        r[2] = input[2] + iZ*_layerSize;
        _indicatorF(output,r);
        if (output[0]) {
          return true;
        }
      }
    }
  }
  return true;
}

template <typename S>
IndicatorInternal3D<S>::IndicatorInternal3D(IndicatorF3D<S>& indicatorF, S layerSize)
  :  _indicatorF(indicatorF), _layerSize(layerSize)
{
  this->_myMin = indicatorF.getMin() + layerSize;
  this->_myMax = indicatorF.getMax() - layerSize;
}

// returns true if x is inside the layer
template <typename S>
bool IndicatorInternal3D<S>::operator()(bool output[], const S input[])
{
  output[0] = false;
  _indicatorF(output,input);
  if (!output[0])
    return true;

  S r[3];
  for (int iX =- 1; iX < 2; ++iX) {
    for (int iY =- 1; iY < 2; ++iY) {
      for (int iZ =- 1; iZ < 2; ++iZ) {
        r[0] = input[0] + iX*_layerSize;
        r[1] = input[1] + iY*_layerSize;
        r[2] = input[2] + iZ*_layerSize;
        _indicatorF(output,r);
        if (!output[0]) {
          return true;
        }
      }
    }
  }
  return true;
}


/// Indicator function for a cylinder
template <typename S>
IndicatorCylinder3D<S>::IndicatorCylinder3D(Vector<S,3> center1,
    Vector<S,3> center2, S radius)
  :  _center1(center1), _center2(center2), _radius2(radius*radius)
{
  // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0
  init();
}

/// Indicator function for a cylinder
template <typename S>
IndicatorCylinder3D<S>::IndicatorCylinder3D(Vector<S,3> center1,
    Vector<S,3> normal, S radius, S eps)
  :  _center1(center1), _center2(center1), _radius2(radius*radius)
{
  normal.normalize();
  // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0
  _center1 -= .5*eps*normal;
  _center2 = _center1 + eps*normal;

  init();
}

/// Indicator function for a cylinder
template <typename S>
IndicatorCylinder3D<S>::IndicatorCylinder3D(IndicatorCircle3D<S> const& circleF, S eps)
  :  _center1(circleF.getCenter()), _center2(circleF.getCenter()),
     _radius2(circleF.getRadius()*circleF.getRadius())
{
  // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0
  _center1 -= .5*eps*circleF.getNormal();
  _center2 = _center1 + eps*circleF.getNormal();

  init();
}

// returns true if x is inside the cylinder
template <typename S>
bool IndicatorCylinder3D<S>::operator()(bool output[], const S input[])
{
  S X = _I[0]*(input[0]-_center1[0]) + _I[1]*(input[1]-_center1[1]) + _I[2]*(input[2]-_center1[2]);
  S Y = _J[0]*(input[0]-_center1[0]) + _J[1]*(input[1]-_center1[1]) + _J[2]*(input[2]-_center1[2]);
  S Z = _K[0]*(input[0]-_center1[0]) + _K[1]*(input[1]-_center1[1]) + _K[2]*(input[2]-_center1[2]);

  // X^2 + Y^2 <= _radius2
  output[0] = ( Z <= _length && Z >= 0 && X*X + Y*Y <= _radius2 );
  return output[0];
}

template <typename S>
void IndicatorCylinder3D<S>::init()
{
  _length = sqrt( (_center2[0]-_center1[0]) * (_center2[0]-_center1[0])
                  +(_center2[1]-_center1[1]) * (_center2[1]-_center1[1])
                  +(_center2[2]-_center1[2]) * (_center2[2]-_center1[2]) );

  // _K = centre2 - centre1 (normalized)
  _K = {(_center2[0]-_center1[0]) / _length, (_center2[1]-_center1[1]) / _length,
        (_center2[2]-_center1[2]) / _length
       };

  // _I and _J form an orthonormal base with _K
  if ( util::approxEqual(_center2[1],_center1[1]) && util::approxEqual(_center2[0],_center1[0]) ) {
    if ( util::approxEqual(_center2[2],_center1[2]) ) {
      OstreamManager clout = OstreamManager(std::cout,"IndicatorCylinder3D");
      clout << "Warning: in the cylinder, the two centers have the same coordinates" << std::endl;
      clout << _center1 << std::endl;
      clout << _center2 << std::endl;
    }
    _I = {1,0,0};
    _J = {0,1,0};
  } else {
    S normi = sqrt (_K[1]*_K[1]+_K[0]*_K[0]);
    _I = {-_K[1]/normi, _K[0]/normi,0};
    _J = {_K[1]*_I[2] - _K[2]*_I[1], _K[2]*_I[0] - _K[0]*_I[2], _K[0]*_I[1] - _K[1]*_I[0]};
  }

  double r = sqrt(_radius2);
  S maxx, maxy, maxz;
  S minx, miny, minz;

  maxx= _center1[0] + sqrt(_I[0]*_I[0]+_J[0]*_J[0])*r + max(_K[0]*_length, 0.);
  minx= _center1[0] - sqrt(_I[0]*_I[0]+_J[0]*_J[0])*r + min(_K[0]*_length, 0.);

  maxy= _center1[1] + sqrt(_I[1]*_I[1]+_J[1]*_J[1])*r + max(_K[1]*_length, 0.);
  miny= _center1[1] - sqrt(_I[1]*_I[1]+_J[1]*_J[1])*r + min(_K[1]*_length, 0.);

  maxz= _center1[2] + sqrt(_I[2]*_I[2]+_J[2]*_J[2])*r + max(_K[2]*_length, 0.);
  minz= _center1[2] - sqrt(_I[2]*_I[2]+_J[2]*_J[2])*r + min(_K[2]*_length, 0.);

  this->_myMin = {minx, miny, minz};
  this->_myMax = {maxx, maxy, maxz};
}

//SM - Analytical cylinder distance method definitions
//Currently just for infinite cylinder 

template <typename S>
bool IndicatorCylinder3D<S>::distance(S& distance, const Vector<S,3>& origin,
                const Vector<S,3>& direction, int iC)
{
  const S r = sqrt(_radius2);
  const Vector <S, 3> cylOrigin = getCenter1();

  //Ray - cylinder intersection
  //Transform origin and direction  (just translation for now)
  const Vector<S,3> originLocal(_I[0]*(origin[0]-cylOrigin[0]) + _I[1]*(origin[1]-cylOrigin[1]) + _I[2]*(origin[2]-cylOrigin[2]),
                                  _J[0]*(origin[0]-cylOrigin[0]) + _J[1]*(origin[1]-cylOrigin[1]) + _J[2]*(origin[2]-cylOrigin[2]),
                                  _K[0]*(origin[0]-cylOrigin[0]) + _K[1]*(origin[1]-cylOrigin[1]) + _K[2]*(origin[2]-cylOrigin[2]));

  Vector<S,3> directionLocal(_I[0]*direction[0] + _I[1]*direction[1] + _I[2]*direction[2],
                             _J[0]*direction[0] + _J[1]*direction[1] + _J[2]*direction[2],
                             _K[0]*direction[0] + _K[1]*direction[1] + _K[2]*direction[2]);

  //Normalise local direction vector
  S dirMag = directionLocal[0]*directionLocal[0] + directionLocal[1]*directionLocal[1] + directionLocal[2]*directionLocal[2];
  if (dirMag >= 0.) {
    directionLocal[0] = directionLocal[0]/sqrt(dirMag);
    directionLocal[1] = directionLocal[1]/sqrt(dirMag);
    directionLocal[2] = directionLocal[2]/sqrt(dirMag);
  }
  else {
    //std::cout << "ERROR: DIRECTION OF ZERO SPECIFIED" << endl;
    distance = -pow(10,6);
    return true;
  }

  S a = directionLocal[0] * directionLocal[0] + directionLocal[1] * directionLocal[1]; 
  S b = 2. * originLocal[0] * directionLocal[0] + 2. * originLocal[1] * directionLocal[1];
  S c = originLocal[0] * originLocal[0] + originLocal[1] * originLocal[1] - r * r;
  S sol1 = -1.;
  S sol2 = -1.;

  S discriminant = b*b - 4.*a*c;

  if (util::nearZero(discriminant))
    discriminant = 0.0;

  if (discriminant >= 0. && a != 0.) {
    sol1 = (-b-pow(discriminant,0.5))/(2.*a);
    sol2 = (-b+pow(discriminant,0.5))/(2.*a);
    
    //Both positive - outside cylinder - take smallest (incoming intersection)
    if (( sol1 >= 0.) && (sol2 >= 0.)) {
      distance = std::min(sol1, sol2);
    }
     //Both near zero
    else if (util::nearZero(sol1) && util::nearZero(sol2)) {
      distance = 0.;
    }
    //One negative - inside - take positive
    else if (sol1 * sol2 < 0.0 ) {
      distance = std::max(sol1, sol2);
    }
    else {
      std::cout << "UNEXPECTED CYLINDER DISTANCES " << sol1 << " " << sol2 << std::endl;
      distance = -pow(10,6);
    } 
  }
  else {
    std::cout << "No real distance to infinite cylinder found! " << b*b - 4*a*c << endl;
    distance = -pow(10,6);
    }
  return true;
}

template <typename S>
Vector<S,3> const& IndicatorCylinder3D<S>::getCenter1() const
{
  return _center1;
}

template <typename S>
Vector<S,3> const& IndicatorCylinder3D<S>::getCenter2() const
{
  return _center2;
}

template <typename S>
S IndicatorCylinder3D<S>::getRadius() const
{
  return sqrt(_radius2);
}

// cone defined by the centers of the two extremities and the radiuses of the two extremities
// the 2nd radius is optional: if it is not defined, the 2nd center is the vertex of the cone
template <typename S>
IndicatorCone3D<S>::IndicatorCone3D(Vector<S,3> center1, Vector<S,3> center2,
                                    S radius1, S radius2)
  :  _center1(center1), _center2(center2),
     _radius1(radius1), _radius2(radius2)
{
  // _I,_J,_K is the new base where _K is the axe of the cone
  // _K = centre2 - centre1 (normalized)
  _length = sqrt( (_center2[0]-_center1[0]) * (_center2[0]-_center1[0])
                  +(_center2[1]-_center1[1]) * (_center2[1]-_center1[1])
                  +(_center2[2]-_center1[2]) * (_center2[2]-_center1[2]) );
  // _K = centre2 - centre1 (normalized)
  _K = {(_center2[0]-_center1[0]) / _length, (_center2[1]-_center1[1]) / _length,
        (_center2[2]-_center1[2]) / _length
       };

  // _I and _J form an orthonormal base with _K
  if ( util::approxEqual(_center2[1],_center1[1]) && util::approxEqual(_center2[0],_center1[0]) ) {
    if ( util::approxEqual(_center2[2],_center1[2]) ) {
      OstreamManager clout = OstreamManager(std::cout,"IndicatorCone3D");
      clout << "Warning: in the cone, the two centers have the same coordinates" << std::endl;
      clout << _center1 << std::endl;
      clout << _center2 << std::endl;
    }
    _I = {1,0,0};
    _J = {0,1,0};
  } else {
    S normi = sqrt(_K[1]*_K[1] + _K[0]*_K[0]);
    _I = {-_K[1]/normi, _K[0]/normi,0};
    _J = {_K[1]*_I[2] - _K[2]*_I[1], _K[2]*_I[0] - _K[0]*_I[2], _K[0]*_I[1] - _K[1]*_I[0]};
  }

  S maxx, maxy, maxz;
  S minx, miny, minz;

  maxx= _center1[0] + max( sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius1,
                                sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius2 + _K[0]*_length);
  minx= _center1[0] + min(-sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius1,
                               -sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius2 + _K[0]*_length);

  maxy= _center1[1] + max( sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius1,
                                sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius2 + _K[1]*_length);
  miny= _center1[1] + min(-sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius1,
                               -sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius2 + _K[1]*_length);

  maxz= _center1[2] + max( sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius1,
                                sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius2 + _K[2]*_length);
  minz= _center1[2] + min(-sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius1,
                               -sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius2 + _K[2]*_length);


  this->_myMin = {minx, miny, minz};
  this->_myMax = {maxx, maxy, maxz};
}

// returns true if x is inside the cone(Vector<S,3> center1, Vector<S,3> center2, S radius1
template <typename S>
bool IndicatorCone3D<S>::operator()(bool output[], const S input[])
{
  // radius: the radius of the cone at the point x
  S X = _I[0]*(input[0]-_center1[0]) + _I[1]*(input[1]-_center1[1]) + _I[2]*(input[2]-_center1[2]);
  S Y = _J[0]*(input[0]-_center1[0]) + _J[1]*(input[1]-_center1[1]) + _J[2]*(input[2]-_center1[2]);
  S Z = _K[0]*(input[0]-_center1[0]) + _K[1]*(input[1]-_center1[1]) + _K[2]*(input[2]-_center1[2]);
  S radius = _radius1 + (_radius2 - _radius1)*Z / _length;

  output[0] = ( Z <= _length && Z >= 0 && X*X + Y*Y <= radius*radius );
  return true;
}



// Warning : the cuboid is only defined parallel to the plans x=0, y=0 and z=0 !!!
// For other cuboids, please use the parallelepiped version
template <typename S>
IndicatorCuboid3D<S>::IndicatorCuboid3D(Vector<S,3> extend, Vector<S,3> origin)
  : _center(origin+.5*extend),_xLength(extend[0]), _yLength(extend[1]), _zLength(extend[2])
{
  assert(_xLength>0 && _yLength>0 && _zLength>0);

  this->_myMin = origin;
  this->_myMax = origin + extend;
}

template <typename S>
IndicatorCuboid3D<S>::IndicatorCuboid3D(S xLength, S yLength, S zLength, Vector<S,3> center)
  : _center(center), _xLength(xLength), _yLength(yLength), _zLength(zLength)
{
  assert(_xLength>0 && _yLength>0 && _zLength>0);

  this->_myMin = {_center[0] - _xLength/2., _center[1] - _yLength/2., _center[2] - _zLength/2.};
  this->_myMax = {_center[0] + _xLength/2., _center[1] + _yLength/2., _center[2] + _zLength/2.};
}


template <typename S>
bool IndicatorCuboid3D<S>::operator()(bool output[], const S input[])
{
  // returns true if x is inside the cuboid
  output[0] = ( (fabs(_center[0] - input[0]) < _xLength/2. || util::approxEqual(fabs(_center[0] - input[0]),_xLength/2.))
                && (fabs(_center[1] - input[1]) < _yLength/2. || util::approxEqual(fabs(_center[1] - input[1]),_yLength/2.))
                && (fabs(_center[2] - input[2]) < _zLength/2. || util::approxEqual(fabs(_center[2] - input[2]),_zLength/2.)) );
  return output[0];
}



template <typename S>
IndicatorCuboidRotate3D<S>::IndicatorCuboidRotate3D(Vector<S,3> extend, Vector<S,3> origin, S theta, int plane, Vector<S,3> centerRotation)
  : IndicatorCuboid3D<S>(extend, origin), _theta(theta), _plane(plane), _centerRotation(centerRotation)
{
  assert(plane==0 || plane==1 || plane==2);
}

template <typename S>
IndicatorCuboidRotate3D<S>::IndicatorCuboidRotate3D(S xLength, S yLength, S zLength, Vector<S,3> origin, S theta, int plane, Vector<S,3> centerRotation)
  : IndicatorCuboid3D<S>(xLength, yLength, zLength, origin), _theta(theta), _plane(plane), _centerRotation(centerRotation)
{
  assert(plane==0 || plane==1 || plane==2);
}

// do transformation to axis aligned cuboid, then call operator() of basic cuboid.
template <typename S>
bool IndicatorCuboidRotate3D<S>::operator()(bool output[], const S input[])
{
  //initialize for _plane == 2
  int i=0;
  int j=1;
   if (_plane == 1) { // rotation around y axis
    i=0;
    j=2;
  } else if (_plane == 0) {  // rotation around x axis
    i=1;
    j=2;
  }
  // translation to _centerRotation
  S x = input[i] - _centerRotation[i];
  S y = input[j] - _centerRotation[j];
  // rotation of _theta in rad
  S xx = x*cos(_theta) - y*sin(_theta);
  S yy = x*sin(_theta) + y*cos(_theta);
  // change back to standard coordinate system
  x = xx + _centerRotation[i];
  y = yy + _centerRotation[j];
  S newInput[3];
  newInput[_plane] = input[_plane];
  newInput[i] = x;
  newInput[j] = y;
  IndicatorCuboid3D<S>::operator()(output, newInput);
  return output[0];
}

//SM - Indicator functor for DCA blade aligned with Z axis

template<typename S>
IndicatorBladeDca3D<S>::IndicatorBladeDca3D(Vector<S,3> origin, S chord, S thickness, S span, S radius1, S radius2, S xp, S theta)
  :  _origin(origin), _chord(chord), _thickness(thickness), _span(span), _radius1(radius1), _radius2(radius2), _xp(xp), _theta(theta)
{
  init();
}

//Init
template<typename S>
void IndicatorBladeDca3D<S>::init()
{ 
  //Base vectors to define blade-local axes 
  _I = {cos(_theta*M_PI/180),sin(_theta*M_PI/180),0.};
  _J = {-sin(_theta*M_PI/180),cos(_theta*M_PI/180),0.};
  _K = {0.,0.,1.};

  //DCA circle centres in local reference frame
  _xc2 = _chord/2 - _radius2;
  _yc1 = _thickness/2 - _radius1;

  //Basic min/max approx for now - corners of bounding box
  this-> _myMin = {_origin[0]-0.5*_chord,_origin[1]-0.5*_thickness,_origin[2]};
  this-> _myMax = {_origin[0]+0.5*_chord,_origin[1]+0.5*_thickness,_origin[2]+_span};
}

//operator()  
template<typename S>
bool IndicatorBladeDca3D<S>::operator()(bool output[], const S input[]) 
{
  //Transform input position to blade-local coordinates
  S inputLocalX = _I[0]*(input[0]-_origin[0]) + _I[1]*(input[1]-_origin[1]) + _I[2]*(input[2]-_origin[2]);
  S inputLocalY = _J[0]*(input[0]-_origin[0]) + _J[1]*(input[1]-_origin[1]) + _J[2]*(input[2]-_origin[2]);
  S inputLocalZ = _K[0]*(input[0]-_origin[0]) + _K[1]*(input[1]-_origin[1]) + _K[2]*(input[2]-_origin[2]);

  output[0] = false;
  
  if ((inputLocalZ >= 0) && (inputLocalZ <= _span)) {
    if (inputLocalX == -_chord/2 ) {
      if(inputLocalY == 0.) {
        output[0] = true;
	return true;
      }
    }
    else if ((inputLocalX > -_chord/2) && (inputLocalX < -_xp)) {
      if ((inputLocalY >= -sqrt(_radius2*_radius2-(inputLocalX+_xc2)*(inputLocalX+_xc2))) && (inputLocalY <= sqrt(_radius2*_radius2-(inputLocalX+_xc2)*(inputLocalX+_xc2)))) {
        output[0] = true;
	return true;
      }
    }
    else if ((inputLocalX >= -_xp) && (inputLocalX <= _xp)) {
      if ((inputLocalY >= -_yc1-sqrt(_radius1*_radius1-inputLocalX*inputLocalX)) && (inputLocalY <= _yc1+sqrt(_radius1*_radius1-inputLocalX*inputLocalX))) {
        output[0] = true;
	return true;
      }
    }
    else if ((inputLocalX > _xp) && (inputLocalX < _chord/2)) {
      if ((inputLocalY >= -sqrt(_radius2*_radius2-(inputLocalX-_xc2)*(inputLocalX-_xc2))) && (inputLocalY <= sqrt(_radius2*_radius2-(inputLocalX-_xc2)*(inputLocalX-_xc2)))) {
        output[0] = true;
	return true;
      }
    }
    else if (inputLocalX == _chord/2) {
      if (inputLocalY == 0.) {
        output[0] = true;
	return true;
      }
    }
  }
  return true;
}

/* OLD
//Analytical distance functor for DCA blade
//Only works for points outside blade - returns -1 if inside
template<typename S>
bool IndicatorBladeDca3D<S>::distance(S& distance, const Vector<S,3>& origin,
                                     const Vector<S,3>& direction, int iC)
{
  //std::cout << "-----------------------------------------------------" << endl;
  //Convert input origin and direction to blade-local coordinates
  const Vector<S,3> originLocal((origin[0]-_origin[0])*cos(_theta*M_PI/180) + (origin[1]-_origin[1])*sin(_theta*M_PI/180),
                                  -(origin[0]-_origin[0])*sin(_theta*M_PI/180)  + (origin[1]-_origin[1])*cos(_theta*M_PI/180),
                                  origin[2] - _origin[2]);

  Vector<S,3> directionLocal(direction[0]*cos(_theta*M_PI/180) + direction[1]*sin(_theta*M_PI/180),
  			       -direction[0]*sin(_theta*M_PI/180)  + direction[1]*cos(_theta*M_PI/180),
                               direction[2]);

  //std::cout << "DIRECTION LOCAL " << directionLocal[0] << " " << directionLocal[1] << " " << directionLocal[2] << endl;

  //Normalise local direction vector
  S dirMag = directionLocal[0]*directionLocal[0] + directionLocal[1]*directionLocal[1] + directionLocal[2]*directionLocal[2];
  //std::cout << "DIRMAG" << dirMag << endl;
  if (dirMag >= 0.) {
    directionLocal[0] = directionLocal[0]/sqrt(dirMag);
    directionLocal[1] = directionLocal[1]/sqrt(dirMag);
    directionLocal[2] = directionLocal[2]/sqrt(dirMag);
  }
  else {
    //std::cout << "ERROR: DIRECTION OF ZERO SPECIFIED" << endl;
    distance = -1.;
    return true;
  }

  //std::cout << "ORIGIN LOCAL " << originLocal[0] << " " << originLocal[1] << " " << originLocal[2] << endl; 
  //std::cout << "DIRECTION LOCAL " << directionLocal[0] << " " << directionLocal[1] << " " << directionLocal[2] << endl;

  //Lambdas to compute distance from cylinder and plane
  auto cylDistance = [](S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction,
  		       const Vector<S,3>& x1, const Vector<S,3>& x2, const S& r)
  {
    //std::cout << "ORIGIN CYL LAMBDA" << origin[0] << " " << origin[1] << " " << origin[2] << endl;
    //Analytical solution via 3D point-line distance formula
    S s1 = -x1[2]*direction[0]-direction[2]*x2[0]+x1[0]*direction[2]+direction[0]*x2[2];
    S s2 = -x1[1]*direction[2]-direction[1]*x2[2]+x1[2]*direction[1]+direction[2]*x2[1];
    S s3 = -x1[0]*direction[1]-direction[0]*x2[1]+x1[1]*direction[0]+direction[1]*x2[0];
			            
    S s4 = -origin[1]*x2[2]-x1[1]*origin[2]+x1[1]*x2[2]+origin[2]*x2[1]+x1[2]*origin[1]-x1[2]*x2[1];
    S s5 = -origin[2]*x2[0]-x1[2]*origin[0]+x1[2]*x2[0]+origin[0]*x2[2]+x1[0]*origin[2]-x1[0]*x2[2];
    S s6 = -origin[0]*x2[1]-x1[0]*origin[1]+x1[0]*x2[1]+origin[1]*x2[0]+x1[1]*origin[0]-x1[1]*x2[0];
    S s7 = r*r*((x2[0]-x1[0])*(x2[0]-x1[0])+(x2[1]-x1[1])*(x2[1]-x1[1])+(x2[2]-x1[2])*(x2[2]-x1[2]));
								            
    S a = s1*s1 + s2*s2 + s3*s3;
    S b = 2*(s4*s2+s1*s5+s3*s6);
    S c = s4*s4+s5*s5+s6*s6-s7;
    S sol1 = -1.;
    S sol2 = -1.;

    S discriminant = b*b - 4.*a*c;

    if (discriminant >= 0.) {
      sol1 = (-b-pow(discriminant,0.5))/(2.*a);
      sol2 = (-b+pow(discriminant,0.5))/(2.*a);
    }
    else if (discriminant > -pow(10,-15)) {
      discriminant = 0.;
      //std::cout << "Set discriminant to zero!" << endl;
      sol1 = (-b-pow(discriminant,0.5))/(2.*a);
      sol2 = (-b+pow(discriminant,0.5))/(2.*a);
    }
    else {
      //std::cout << "No real distance to infinite cylinder found! " << b*b - 4*a*c << endl;
      distance = -1.;
      return true;
    }
    if (sol1 >= 0. && sol2 >= 0.) {
      distance = min(sol1,sol2);
      //std::cout << "MIN(SOL1,SOL2)= " << distance << endl;
    }
    else {
      //std::cout << "Origin point is inside cylinder " << sol1 << " " << sol2 << endl;
      distance = -1.;
      return true;
    }
    if (origin[2] + distance*direction[2] >= x1[2] - pow(10,-15) && origin[2] + distance*direction[2] <= x2[2] + pow(10,-15))     {
      return true;
    }
    else {
      //std::cout << "Distance is out of finite cylinder bounds " << origin[2] + distance*direction[2] << " " << x1[2] << " " << x2[2] << endl;
      distance = -1.;
      return true;
    }
  };

  //Plane
  auto planeDistance = [](S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction,
  		       const Vector<S,3>& pOrigin, const Vector<S,3>& pNormal)
  { 
    //Ray-plane intersection formula
    S denom = direction[0]*pNormal[0] + direction[1]*pNormal[1] + direction[2]*pNormal[2];
    if (denom != 0.) {
      distance = ((origin[0]-pOrigin[0])*pNormal[0]+(origin[1]-pOrigin[1])*pNormal[1]+(origin[2]-pOrigin[2])*pNormal[2])/denom;
      return true;
    }
    else {
      //std::cout << "Error: Direction perpendicular to plane" << endl;
      distance = -1;
      return true;
    }
  };

  //Compare distances from different fundamental components of blade geometry
  S d1, d2;
  std::vector<S> cylDistances(4);
  //Check if origin is bounded by blade in x-y
  bool isInsideXY = false;

  if (originLocal[0] == -_chord/2 ) {
    if(originLocal[1] == 0.) {
      isInsideXY = true;
    }
  }
  else if ((originLocal[0] > -_chord/2) && (originLocal[0] < -_xp)) {
    if ((originLocal[1] >= -sqrt(_radius2*_radius2-(originLocal[0]+_xc2)*(originLocal[0]+_xc2))) && (originLocal[1] <= sqrt(_radius2*_radius2-(originLocal[0]+_xc2)*(originLocal[0]+_xc2)))) {
      isInsideXY = true;
    }
  }
  else if ((originLocal[0] >= -_xp) && (originLocal[0] <= _xp)) {
    if ((originLocal[1] >= -_yc1-sqrt(_radius1*_radius1-originLocal[0]*originLocal[0])) && (originLocal[1] <= _yc1+sqrt(_radius1*_radius1-originLocal[0]*originLocal[0]))) {
      isInsideXY = true;
    }
  }
  else if ((originLocal[0] > _xp) && (originLocal[0] < _chord/2)) {
    if ((originLocal[1] >= -sqrt(_radius2*_radius2-(originLocal[0]-_xc2)*(originLocal[0]-_xc2))) && (originLocal[1] <= sqrt(_radius2*_radius2-(originLocal[0]-_xc2)*(originLocal[0]-_xc2)))) {
      isInsideXY = true;
    }
  }
  else if (originLocal[0] == _chord/2) {
    if (originLocal[1] == 0.) {
      isInsideXY = true;
      }
  }
  //std::cout << isInsideXY << endl;
  //std::cout << "posx vs -xp= " << originLocal[0]+_xp << endl;
  //std::cout << "posx vs xp=" << originLocal[0]-_xp << endl;
  //std::cout << "posx vs chord" << originLocal[0]+_chord/2 << endl;
  //std::cout << "posx vs chord" << originLocal[0]-_chord/2 << endl;

  if (isInsideXY && originLocal[2] <= 0.) {
    planeDistance(d1, originLocal, directionLocal, Vector<S,3>{0.,0.,0.}, Vector<S,3>{0.,0.,-1.});
    distance = d1;
    //std::cout << "D1" << endl;
    return true; 
  }
  else if (isInsideXY && originLocal[2] >= _span) {
   planeDistance(d2, originLocal, directionLocal, Vector<S,3>{0.,0.,_span}, Vector<S,3>{0.,0.,1.});
    distance = d2;
    //std::cout << "D2" << endl;
    return true; 
  }
  else {
    cylDistance(cylDistances[0], originLocal, directionLocal, Vector<S,3>{-_xc2,0.,0.}, Vector<S,3>{-_xc2,0.,_span},_radius2);
    S intersectCyl3x = originLocal[0] + cylDistances[0]*directionLocal[0];
    if (intersectCyl3x < -_chord/2. - pow(10,-15) || intersectCyl3x > -_xp + pow(10,-15.)) {
      cylDistances[0] = -1.; //Distance invalid as x intersect is outside of the cylinder bounds
      //std::cout << "d3 invalid " << intersectCyl3x + _chord/2 << " "  << intersectCyl3x + _xp << endl;
    }
    cylDistance(cylDistances[1], originLocal, directionLocal, Vector<S,3>{0.,-_yc1,0.},Vector<S,3>{0.,-_yc1,_span},_radius1);
    S intersectCyl4x = originLocal[0] + cylDistances[1]*directionLocal[0];
    if (intersectCyl4x <= -_xp || intersectCyl4x >= _xp) {
      cylDistances[1] = -1.; //Distance invalid as x intersect is outside of the cylinder bounds
      //std::cout << "d4 invalid " << intersectCyl4x + _xp << " " << intersectCyl4x - _xp << endl;
    }
    cylDistance(cylDistances[2], originLocal, directionLocal, Vector<S,3>{0.,_yc1,0.}, Vector<S,3>{0.,_yc1,_span},_radius1); 
    S intersectCyl5x = originLocal[0] + cylDistances[2]*directionLocal[0];
    if (intersectCyl5x <= -_xp || intersectCyl5x >= _xp) {
      cylDistances[2] = -1.; //Distance invalid as x intersect is outside of the cylinder bounds
      //std::cout << "d5 invalid " << intersectCyl5x + _xp << " " << intersectCyl5x - _xp << endl;
    }
    cylDistance(cylDistances[3], originLocal, directionLocal, Vector<S,3>{_xc2,0.,0.}, Vector<S,3>{_xc2,0.,_span},_radius2);
    S intersectCyl6x = originLocal[0] + cylDistances[3]*directionLocal[0];
    if (intersectCyl6x < _xp - pow(10,-15) || intersectCyl6x > _chord/2. + pow(10,-15.)) {
      cylDistances[3] = -1.; //Distance invalid as x intersect is outside of the cylinder bounds
      //std::cout << "d6 invalid " << intersectCyl6x - _chord/2 << " " << intersectCyl6x - _xp << endl;
    }
    //std::cout << cylDistances[0] << " " <<  cylDistances[1] << " " << cylDistances[2] << " " <<  cylDistances[3] << endl;
    if (cylDistances[0]<0. && cylDistances[1]<0. && cylDistances[2]<0. && cylDistances[3]<0.) {
      //std::cout << "NO POSITIVE CYL DISTANCE FOUND" << endl;
      distance = -1.;
      return true;
    } else {
      //Get min valid distance
      distance = -1.;
      for (int i=0;i<4;++i) {
        if (distance < 0. && cylDistances[i] >= 0.) {
	  distance = cylDistances[i];
	}
	if (cylDistances[i] > 0. && cylDistances[i]<distance) {
	  distance = cylDistances[i];
	}
      }
      return true;
    }
  }
}
*/


//Analytical distance functor for DCA blade
template<typename S>
bool IndicatorBladeDca3D<S>::distance(S& distance, const Vector<S,3>& origin,
                                     const Vector<S,3>& direction, int iC)
{
  //Convert input origin and direction to blade-local coordinates
  const Vector<S,3> originLocal(_I[0]*(origin[0]-_origin[0]) + _I[1]*(origin[1]-_origin[1]) + _I[2]*(origin[2]-_origin[2]),
                                  _J[0]*(origin[0]-_origin[0]) + _J[1]*(origin[1]-_origin[1]) + _J[2]*(origin[2]-_origin[2]),
                                  _K[0]*(origin[0]-_origin[0]) + _K[1]*(origin[1]-_origin[1]) + _K[2]*(origin[2]-_origin[2]));

  Vector<S,3> directionLocal(_I[0]*direction[0] + _I[1]*direction[1] + _I[2]*direction[2],
                             _J[0]*direction[0] + _J[1]*direction[1] + _J[2]*direction[2],
                             _K[0]*direction[0] + _K[1]*direction[1] + _K[2]*direction[2]);

  //Normalise local direction vector
  S dirMag = directionLocal[0]*directionLocal[0] + directionLocal[1]*directionLocal[1] + directionLocal[2]*directionLocal[2];
  if (dirMag >= 0.) {
    directionLocal[0] = directionLocal[0]/sqrt(dirMag);
    directionLocal[1] = directionLocal[1]/sqrt(dirMag);
    directionLocal[2] = directionLocal[2]/sqrt(dirMag);
  }
  else {
    //std::cout << "ERROR: DIRECTION OF ZERO SPECIFIED" << endl;
    distance = -pow(10,6);
    return true;
  }
  //Lambdas to compute distance from cylinder and plane
  auto cylDistance = [](Vector<S,2>& distance, const Vector<S,3>& origin, const Vector<S,3>& direction,
  		       const Vector<S,3>& x1, const Vector<S,3>& x2, const S& r)
  {
    //std::cout << "ORIGIN CYL LAMBDA" << origin[0] << " " << origin[1] << " " << origin[2] << endl;
    //Analytical solution via 3D point-line distance formula
    S s1 = -x1[2]*direction[0]-direction[2]*x2[0]+x1[0]*direction[2]+direction[0]*x2[2];
    S s2 = -x1[1]*direction[2]-direction[1]*x2[2]+x1[2]*direction[1]+direction[2]*x2[1];
    S s3 = -x1[0]*direction[1]-direction[0]*x2[1]+x1[1]*direction[0]+direction[1]*x2[0];
			            
    S s4 = -origin[1]*x2[2]-x1[1]*origin[2]+x1[1]*x2[2]+origin[2]*x2[1]+x1[2]*origin[1]-x1[2]*x2[1];
    S s5 = -origin[2]*x2[0]-x1[2]*origin[0]+x1[2]*x2[0]+origin[0]*x2[2]+x1[0]*origin[2]-x1[0]*x2[2];
    S s6 = -origin[0]*x2[1]-x1[0]*origin[1]+x1[0]*x2[1]+origin[1]*x2[0]+x1[1]*origin[0]-x1[1]*x2[0];
    S s7 = r*r*((x2[0]-x1[0])*(x2[0]-x1[0])+(x2[1]-x1[1])*(x2[1]-x1[1])+(x2[2]-x1[2])*(x2[2]-x1[2]));
								            
    S a = s1*s1 + s2*s2 + s3*s3;
    S b = 2*(s4*s2+s1*s5+s3*s6);
    S c = s4*s4+s5*s5+s6*s6-s7;
    S sol1 = -1.;
    S sol2 = -1.;

    S discriminant = b*b - 4.*a*c;

    if (discriminant >= 0. && a > 0.) {
      sol1 = (-b-pow(discriminant,0.5))/(2.*a);
      sol2 = (-b+pow(discriminant,0.5))/(2.*a);
      distance[0] = sol1;
      distance[1] = sol2;
      return true;
    }
    //else if (discriminant > -pow(10,-15) && a > 0.) {
    //  discriminant = 0.;
      //std::cout << "Set discriminant to zero!" << endl;
    //  sol1 = (-b-pow(discriminant,0.5))/(2.*a);
    //  sol2 = (-b+pow(discriminant,0.5))/(2.*a);
    //  distance[0] = sol1;
    //  distance[1] = sol2;
    //}
    else {
      //std::cout << "No real distance to infinite cylinder found! " << b*b - 4*a*c << endl;
      distance[0] = -pow(10,6);
      distance[1] = -pow(10,6);
      return true;
    }
    return true;
  };
  
  //Plane
  auto planeDistance = [](S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction,
  		       const Vector<S,3>& pOrigin, const Vector<S,3>& pNormal)
  { 
    //Ray-plane intersection formula
    S denom = direction[0]*pNormal[0] + direction[1]*pNormal[1] + direction[2]*pNormal[2];
    if (denom != 0.) {
      distance = ((origin[0]-pOrigin[0])*pNormal[0]+(origin[1]-pOrigin[1])*pNormal[1]+(origin[2]-pOrigin[2])*pNormal[2])/denom;
      return true;
    }
    else {
      //std::cout << "Error: Direction perpendicular to plane" << endl;
      distance = -pow(10,6);
      return true;
    }
  };

  //Lambda to compare if point is inside X-Y bounds of blade
  auto isInsideXY = [this](bool& isInside, const Vector<S,3>& origin)
  { 
    isInside = false; 
    if (origin[0] == -_chord/2 ) {
      if(origin[1] == 0.) {
        isInside = true;
      }
    }
    else if ((origin[0] > -_chord/2) && (origin[0] < -_xp)) {
      if ((origin[1] >= -sqrt(_radius2*_radius2-(origin[0]+_xc2)*(origin[0]+_xc2))) && (origin[1] <= sqrt(_radius2*_radius2-(origin[0]+_xc2)*(origin[0]+_xc2)))) {
        isInside = true;
      }
    }
    else if ((origin[0] >= -_xp) && (origin[0] <= _xp)) {
      if ((origin[1] >= -_yc1-sqrt(_radius1*_radius1-origin[0]*origin[0])) && (origin[1] <= _yc1+sqrt(_radius1*_radius1-origin[0]*origin[0]))) {
        isInside = true;
      }
    }
    else if ((origin[0] > _xp) && (origin[0] < _chord/2)) {
      if ((origin[1] >= -sqrt(_radius2*_radius2-(origin[0]-_xc2)*(origin[0]-_xc2))) && (origin[1] <= sqrt(_radius2*_radius2-(origin[0]-_xc2)*(origin[0]-_xc2)))) {
        isInside = true;
      }
    }
    else if (origin[0] == _chord/2) {
      if (origin[1] == 0.) {
        isInside = true;
      }
    }
  return true;
  };

  bool isInside;
  S d1, d2;
  Vector<S,2> d3;
  Vector<S,2> d4;
  Vector<S,2> d5;
  Vector<S,2> d6;

  planeDistance(d1, originLocal, directionLocal, Vector<S,3>{0.,0.,0.}, Vector<S,3>{0.,0.,-1.});
  Vector<S,3> intersect1;
  for(int i=0; i<3; ++i) {
    intersect1[i] = originLocal[i] + d1 * directionLocal[i];
  }
  isInsideXY(isInside,intersect1);
  if(!isInside) {
    d1 = -pow(10,6);
  }
  planeDistance(d2, originLocal, directionLocal, Vector<S,3>{0.,0.,_span}, Vector<S,3>{0.,0.,1.});
  Vector<S,3> intersect2;
  for(int i=0; i<3; ++i) {
    intersect2[i] = originLocal[i] + d2 * directionLocal[i];
  }
  isInsideXY(isInside,intersect2);
  if(!isInside) {
    d2 = -pow(10,6);
  }
  cylDistance(d3, originLocal, directionLocal, Vector<S,3>{-_xc2,0.,0.}, Vector<S,3>{-_xc2,0.,_span},_radius2);
  Vector<S,2> intersectCyl3x;
  for (int i=0;i<2;++i){
    intersectCyl3x[i] = originLocal[0] + d3[i] * directionLocal[0];
    if (intersectCyl3x[i] < -_chord/2. - pow(10,-15) || intersectCyl3x[i] > -_xp + pow(10,-15.)) {
      d3[i] = -pow(10,6); //Distance invalid as x intersect is outside of the cylinder bounds
    }
  }
  cylDistance(d4, originLocal, directionLocal, Vector<S,3>{0.,-_yc1,0.},Vector<S,3>{0.,-_yc1,_span},_radius1);
  Vector<S,2> intersectCyl4x;
  for (int i=0;i<2;++i){
    intersectCyl4x[i] = originLocal[0] + d4[i] * directionLocal[0];
    if (intersectCyl4x[i] < -_xp || intersectCyl4x[i] > _xp) {
      d4[i] = -pow(10,6); //Distance invalid as x intersect is outside of the cylinder bounds
    }
  }
  cylDistance(d5, originLocal, directionLocal, Vector<S,3>{0.,_yc1,0.}, Vector<S,3>{0.,_yc1,_span},_radius1); 
  Vector<S,2> intersectCyl5x;
  for (int i=0;i<2;++i){
    intersectCyl5x[i] = originLocal[0] + d5[i] * directionLocal[0];
    if (intersectCyl5x[i] < -_xp || intersectCyl5x[i] > _xp) {
      d5[i] = -pow(10,6); //Distance invalid as x intersect is outside of the cylinder bounds
    }
  }
  cylDistance(d6, originLocal, directionLocal, Vector<S,3>{_xc2,0.,0.}, Vector<S,3>{_xc2,0.,_span},_radius2);
  Vector<S,2> intersectCyl6x;
  for (int i=0;i<2;++i){
    intersectCyl6x[i] = originLocal[0] + d6[i] * directionLocal[0];
    if (intersectCyl6x[i] < _xp - pow(10,-15) || intersectCyl6x[i] > _chord/2. + pow(10,-15.)) {
      d6[i] = -pow(10,6); //Distance invalid as x intersect is outside of the cylinder bounds
    }
  }
  std::vector<S> absDistances;
  absDistances.assign({abs(d1),abs(d2),abs(d3[0]),abs(d3[1]),abs(d4[0]),abs(d4[1]),abs(d5[0]),abs(d5[1]),abs(d6[0]),abs(d6[1])});
  distance = *std::min_element(std::begin(absDistances),std::end(absDistances));
  //std::cout << "---" << std::endl;  
  //std::cout << "ORIGIN LOCAL" << originLocal[0] << " " << originLocal[1] << " " << originLocal[2] << endl;;
  //std::cout << "DIRECTION LOCAL" << directionLocal[0] << " " << directionLocal[1] << " " << directionLocal[2] << endl;
  //std::cout << "DISTANCE" << distance << endl;
  return true;
}

template<typename S> 
Vector<S,3> const&  IndicatorBladeDca3D<S>::getOrigin() const 
{
  return _origin;
}

template<typename S> 
S IndicatorBladeDca3D<S>::getChord() const 
{
  return _chord;
}

template<typename S> 
S IndicatorBladeDca3D<S>::getThickness() const 
{
  return _thickness;
}

template<typename S> 
S IndicatorBladeDca3D<S>::getSpan() const 
{
  return _span;
}

template<typename S> 
S IndicatorBladeDca3D<S>::getTheta() const 
{
  return _theta;
}

////// creator functions /////
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorCircle3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCircle3D");

  Vector<S,3> center;
  Vector<S,3> normal;
  S radius = 1;

  //  params.setWarningsOn(false);
  params.setWarningsOn(true);

  std::stringstream xmlCenter1( params.getAttribute("center") );
  xmlCenter1 >> center[0] >> center[1] >> center[2];
  std::stringstream xmlCenter2( params.getAttribute("normal") );
  xmlCenter2 >> normal[0] >> normal[1] >> normal[2];
  std::stringstream xmlRadius( params.getAttribute("radius") );
  xmlRadius >> radius;

  return std::make_shared<IndicatorCircle3D<S>>(center, normal, radius);
}

template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorSphere3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorSphere3D");

  Vector<S,3> center;
  S radius = 1;

  //  params.setWarningsOn(false);
  params.setWarningsOn(true);

  std::stringstream xmlCenter1( params.getAttribute("center") );
  xmlCenter1 >> center[0] >> center[1] >> center[2];
  std::stringstream xmlRadius( params.getAttribute("radius") );
  xmlRadius >> radius;

  return std::make_shared<IndicatorSphere3D<S>>(center, radius);
}

// creator function for a cylinder3d
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorCylinder3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCylinder3D");

  Vector<S,3> center1;
  Vector<S,3> center2(S(1),S(1),S(1));
  S radius = 1;

  //  params.setWarningsOn(false);
  //  params.setWarningsOn(true);

  std::stringstream xmlCenter1( (params).getAttribute("center1") );
  xmlCenter1 >> center1[0] >> center1[1] >> center1[2];
  std::stringstream xmlCenter2( (params).getAttribute("center2") );
  xmlCenter2 >> center2[0] >> center2[1] >> center2[2];
  std::stringstream xmlRadius( (params).getAttribute("radius") );
  xmlRadius >> radius;

  /// for debugging purpose
//  print(center1, "center1: ");
//  print(center2, "center2: ");
//  print(radius, "radius: ");

  return std::make_shared<IndicatorCylinder3D<S>>(center1, center2, radius);
}

template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorCone3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCone3D");

  Vector<S,3> center1;
  Vector<S,3> center2(S(1), S(1), S(1));
  S radius1 = S(0);
  S radius2 = S(1);

  //  params.setWarningsOn(false);
  params.setWarningsOn(true);

  std::stringstream xmlCenter1( params.getAttribute("center1") );
  xmlCenter1 >> center1[0] >> center1[1] >> center1[2];
  std::stringstream xmlCenter2( params.getAttribute("center2") );
  xmlCenter2 >> center2[0] >> center2[1] >> center2[2];
  std::stringstream xmlRadius1( params.getAttribute("radius1") );
  xmlRadius1 >> radius1;
  std::stringstream xmlRadius2( params.getAttribute("radius2") );
  xmlRadius2 >> radius2;

  return std::make_shared<IndicatorCone3D<S>>(center1, center2, radius1, radius2);
}

template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorCuboid3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCuboid3D");

  Vector<S,3> origin;
  Vector<S,3> extend(S(1),S(1),S(1));

  //  params.setWarningsOn(false);
  params.setWarningsOn(true);

  std::stringstream xmlOrigin( params.getAttribute("origin") );
  xmlOrigin >> origin[0] >> origin[1] >> origin[2];
  std::stringstream xmlExtend( params.getAttribute("extend") );
  xmlExtend >> extend[0] >> extend[1] >> extend[2];

  return std::make_shared<IndicatorCuboid3D<S>>(extend, origin);
}

// Create Union with XML - file
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorUnion3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorUnion3D");

  //  clout << "xml position: " << params.getName() << std::endl;
  //  params.print(2);

  std::shared_ptr<IndicatorF3D<S>> output = createIndicatorF3D<S>(**params.begin());
  for (auto it = params.begin()+1; it != params.end(); ++it) {
    output = output + createIndicatorF3D<S>(**it);
  }
  return output;
}

// Create Without with XML - file
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorWithout3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorWithout3D");

  //  clout << "xml position: " << params.getName() << std::endl;
  //  params.print(2);

  std::shared_ptr<IndicatorF3D<S>> output = createIndicatorF3D<S>(**params.begin());
  for (auto it = params.begin()+1; it != params.end(); ++it) {
    output = output - createIndicatorF3D<S>(**it);
  }
  return output;
}

// Create Intersection with XML - file
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorIntersection3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorIntersection3D");

  //  clout << "xml position: " << params.getName() << std::endl;
  //  params.print(2);

  std::shared_ptr<IndicatorF3D<S>> output = createIndicatorF3D<S>(**params.begin());
  for (auto it = params.begin()+1; it != params.end(); ++it) {
    output = output * createIndicatorF3D<S>(**it);
  }
  return output;
}

// Create Geometry
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorF3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorF3D");

  //  clout << "XML element: "<< params.getName() << std::endl;
  //  params.print(2);

  std::string actualName = params.getName();
  if ( actualName == "IndicatorCircle3D" ) {
    return createIndicatorCircle3D<S>(params);
  } else if ( actualName == "IndicatorSphere3D" ) {
    return createIndicatorSphere3D<S>(params);
  } else if ( actualName == "IndicatorCylinder3D" ) {
    return createIndicatorCylinder3D<S>(params);
  } else if ( actualName == "IndicatorCone3D" ) {
    return createIndicatorCone3D<S>(params);
  } else if ( actualName == "IndicatorCuboid3D" ) {
    return createIndicatorCuboid3D<S>(params);
  } else if ( actualName == "IndicatorUnion3D" ) {
    return createIndicatorUnion3D<S>(params);
  } else if ( actualName == "IndicatorWithout3D" ) {
    return createIndicatorWithout3D<S>(params);
  } else if ( actualName == "IndicatorIntersection3D" ) {
    return createIndicatorIntersection3D<S>(params);
  } else {
    auto firstChild = params.begin(); // get iterator of childTree
    return createIndicatorF3D<S>( **firstChild );
  }
}


} // namespace olb

#endif
