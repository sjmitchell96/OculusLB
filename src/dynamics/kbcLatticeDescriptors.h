/*  This file
 *
 *  Copyright 
 *  
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA 02110-1301, USA.
*/

/** \file
 *  -- header file
 */

#ifndef KBC_LATTICE_DESCRIPTORS_H
#define KBC_LATTICE_DESCRIPTORS_H

#include "latticeDescriptors.h"

namespace olb {

namespace descriptors {

namespace tag {

struct KBC : public CATEGORY, public DESCRIPTOR_TAG { };

} //namespace tag

using D3Q27descriptorKBC = D3Q27<tag::KBC>;
using D3Q27descriptorKBCGrad = D3Q27<tag::KBC,VELOCITY>;

namespace kbc_data {

using utilities::Fraction;

template<unsigned D, unsigned Q>
constexpr Fraction t[Q] = {};

template<>
constexpr Fraction t<3, 27>[27] = {
  {8,27},

  {2, 27}, {2, 27}, {2, 27},
  {1, 54}, {1, 54}, {1, 54},
  {1, 54}, {1, 54}, {1, 54},
  {1, 216}, {1, 216}, {1, 216}, {1, 216},

  {2, 27}, {2, 27}, {2, 27},
  {1, 54}, {1, 54}, {1, 54},
  {1, 54}, {1, 54}, {1, 54},
  {1, 216}, {1, 216}, {1, 216}, {1, 216}
};


}; //namespace kbc_data

template<typename T, unsigned D, unsigned Q>
constexpr T t(unsigned iPop, tag::KBC)
{
  return kbc_data::t<D,Q>[iPop].template as<T>();
}

} //namespace descriptors

} //namespace olb

#endif