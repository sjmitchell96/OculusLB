/*  This file is 
 *
 *  
 *  E-mail contact: 
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

/** \file
 * Template specializations for some computationally intensive LB
 * functions of the header file kbcLbHelpers.h - D3Q27
 */

#ifndef KBC_LB_HELPERS_D3Q27_H
#define KBC_LB_HELPERS_D3Q27_H

#include <cmath>

namespace olb {

// Efficient specialisation for D3Q27 lattice

template<typename T>
struct kbcLbHelpers<T, descriptors::D3Q27descriptorKBC> {

  static T kbcCollision(Cell<T, descriptors::D3Q27descriptorKBC>& cell, T rho, T u[3], const T& beta) {

  typedef descriptors::D3Q27descriptorKBC L;
   
  //std::cout << "SPECIAL KBC COLLIDE" << std::endl;

  //Compute feq
    //Subexpressions
  T subexpr1 = std::sqrt(1. + 3. * u[0] * u[0]);
  T subexpr2 = std::sqrt(1. + 3. * u[1] * u[1]);
  T subexpr3 = std::sqrt(1. + 3. * u[2] * u[2]);

  T aux = 2. - subexpr1;
  T auy = 2. - subexpr2;
  T auz = 2. - subexpr3;

  T bux = (2. * u[0] + subexpr1) / (1. - u[0]);
  T buy = (2. * u[1] + subexpr2) / (1. - u[1]);
  T buz = (2. * u[2] + subexpr3) / (1. - u[2]);

  T buxInv = 1. / bux;
  T buyInv = 1. / buy;
  T buzInv = 1. / buz;

  T buxInvBuy = buxInv * buy; 
  T buxInvBuyInv = buxInv * buyInv;
  T buyInvBuzInv = buyInv * buzInv;
  T buyInvBuz = buyInv * buz;
  T buxBuy = bux * buy;

  T rhoPhi = rho * aux * auy * auz;

  T W0 =  descriptors::t<T,L>(0);
  T W1 =  descriptors::t<T,L>(1);
  T W4 =  descriptors::t<T,L>(4);
  T W10 = descriptors::t<T,L>(10);

  T rhoPhiW1 = rhoPhi * W1; //weight products
  T rhoPhiW4 = rhoPhi * W4; 
  T rhoPhiW10 = rhoPhi * W10; 

  T fEq[27];
  fEq[0] = rhoPhi * W0;
  fEq[1] = rhoPhiW1 * buxInv;
  fEq[2] = rhoPhiW1 * buyInv;
  fEq[3] = rhoPhiW1 * buzInv;
  fEq[4] = rhoPhiW4 * buxInvBuyInv;
  fEq[5] = rhoPhiW4 * buxInvBuy;
  fEq[6] = rhoPhiW4 * buxInv * buzInv;
  fEq[7] = rhoPhiW4 * buxInv * buz;
  fEq[8] = rhoPhiW4 * buyInvBuzInv;
  fEq[9] = rhoPhiW4 * buyInvBuz;
  fEq[10] = rhoPhiW10 * buxInvBuyInv * buzInv;
  fEq[11] = rhoPhiW10 * buxInvBuyInv * buz;
  fEq[12] = rhoPhiW10 * buxInvBuy * buzInv;
  fEq[13] = rhoPhiW10 * buxInvBuy * buz;
  fEq[14] = rhoPhiW1 * bux;
  fEq[15] = rhoPhiW1 * buy;
  fEq[16] = rhoPhiW1 * buz;
  fEq[17] = rhoPhiW4 * buxBuy;
  fEq[18] = rhoPhiW4 * bux * buyInv;
  fEq[19] = rhoPhiW4 * bux * buz;
  fEq[20] = rhoPhiW4 * bux * buzInv;
  fEq[21] = rhoPhiW4 * buy * buz;
  fEq[22] = rhoPhiW4 * buy * buzInv;
  fEq[23] = rhoPhiW10 * buxBuy * buz;
  fEq[24] = rhoPhiW10 * buxBuy * buzInv;
  fEq[25] = rhoPhiW10 * bux * buyInvBuz;
  fEq[26] = rhoPhiW10 * bux * buyInvBuzInv;

  T df[27];
  df[0] = cell[0] - fEq[0] + W0;
  df[1] = cell[1] - fEq[1] + W1; 
  df[2] = cell[2] - fEq[2] + W1; 
  df[3] = cell[3] - fEq[3] + W1; 
  df[4] = cell[4] - fEq[4] + W4; 
  df[5] = cell[5] - fEq[5] + W4; 
  df[6] = cell[6] - fEq[6] + W4; 
  df[7] = cell[7] - fEq[7] + W4; 
  df[8] = cell[8] - fEq[8] + W4; 
  df[9] = cell[9] - fEq[9] + W4; 
  df[10] = cell[10] - fEq[10] + W10; 
  df[11] = cell[11] - fEq[11] + W10; 
  df[12] = cell[12] - fEq[12] + W10; 
  df[13] = cell[13] - fEq[13] + W10; 
  df[14] = cell[14] - fEq[14] + W1; 
  df[15] = cell[15] - fEq[15] + W1; 
  df[16] = cell[16] - fEq[16] + W1; 
  df[17] = cell[17] - fEq[17] + W4; 
  df[18] = cell[18] - fEq[18] + W4; 
  df[19] = cell[19] - fEq[19] + W4; 
  df[20] = cell[20] - fEq[20] + W4; 
  df[21] = cell[21] - fEq[21] + W4; 
  df[22] = cell[22] - fEq[22] + W4; 
  df[23] = cell[23] - fEq[23] + W10; 
  df[24] = cell[24] - fEq[24] + W10; 
  df[25] = cell[25] - fEq[25] + W10; 
  df[26] = cell[26] - fEq[26] + W10; 

  //Moments computations

  T rhoDeltaM002 = df[3] + df[6] + df[7] + df[8] + df[9] +
                   df[10] + df[11] + df[12] + df[13] + df[16] +
                   df[19] + df[20] + df[21] + df[22] + df[23] +
                   df[24] + df[25] + df[26];

  T rhoDeltaM200 = df[1] + df[4] + df[5] + df[6] + df[7] +
                   df[10] + df[11] + df[12] + df[13] + df[14] +
                   df[17] + df[18] + df[19] + df[20] + df[23] +
                   df[24] + df[25] + df[26];

  T rhoDeltaM020 = df[2] + df[4] + df[5] + df[8] + df[9] +
                   df[10] + df[11] + df[12] + df[13] + df[15] +
                   df[17] + df[18] + df[21] + df[22] + df[23] +
                   df[24] + df[25] + df[26];

  T rhoDeltaPiXY = df[4] - df[5] + df[10] + df[11] - df[12] - df[13] +                 
                   df[17] - df[18] + df[23] + df[24] - df[25] - df[26];

  T rhoDeltaPiXZ = df[6] - df[7] + df[10] - df[11] + df[12] - df[13] +                 
                   df[19] - df[20] + df[23] - df[24] + df[25] - df[26];

  T rhoDeltaPiYZ = df[8] - df[9] + df[10] - df[11] - df[12] + df[13] +                 
                   df[21] - df[22] + df[23] - df[24] - df[25] + df[26];

  T invRho = 1. / rho;
  T const1 = 1. / 6.;
  T const2 = 1. / 4.;

  T dNxz = invRho * (rhoDeltaM200 - rhoDeltaM002);
  T dNyz = invRho * (rhoDeltaM020 - rhoDeltaM002);

  //delta s
  T ds1 = const1 * (2. * dNxz - dNyz);
  T ds2 = const1 * (2. * dNyz - dNxz);
  T ds3 = - const1 * (dNxz + dNyz);
  T ds4 = const2 * invRho * rhoDeltaPiXY;
  T ds5 = - ds4;
  T ds6 = const2 * invRho * rhoDeltaPiXZ;
  T ds7 = - ds6;
  T ds8 = const2 * invRho * rhoDeltaPiYZ;
  T ds9 = - ds8;

  T ds14 = ds1;
  T ds15 = ds2;
  T ds16 = ds3;
  T ds17 = ds4;
  T ds18 = ds5;
  T ds19 = ds6;
  T ds20 = ds7;
  T ds21 = ds8;
  T ds22 = ds9;

  T dh[27];
  dh[0] = df[0];
  dh[1] = df[1] - ds1; 
  dh[2] = df[2] - ds2; 
  dh[3] = df[3] - ds3; 
  dh[4] = df[4] - ds4; 
  dh[5] = df[5] - ds5; 
  dh[6] = df[6] - ds6; 
  dh[7] = df[7] - ds7; 
  dh[8] = df[8] - ds8; 
  dh[9] = df[9] - ds9; 
  dh[10] = df[10];
  dh[11] = df[11];
  dh[12] = df[12];
  dh[13] = df[13];
  dh[14] = df[14] - ds14; 
  dh[15] = df[15] - ds15; 
  dh[16] = df[16] - ds16; 
  dh[17] = df[17] - ds17; 
  dh[18] = df[18] - ds18; 
  dh[19] = df[19] - ds19; 
  dh[20] = df[20] - ds20; 
  dh[21] = df[21] - ds21; 
  dh[22] = df[22] - ds22; 
  dh[23] = df[23];
  dh[24] = df[24];
  dh[25] = df[25];
  dh[26] = df[26];

  //Stabiliser (gamma)
  T dhInvFeq[27];
  dhInvFeq[0] = dh[0] /fEq[0]; 
  dhInvFeq[1] = dh[1] /fEq[1]; 
  dhInvFeq[2] = dh[2] /fEq[2]; 
  dhInvFeq[3] = dh[3] /fEq[3]; 
  dhInvFeq[4] = dh[4] /fEq[4]; 
  dhInvFeq[5] = dh[5] /fEq[5]; 
  dhInvFeq[6] = dh[6] /fEq[6]; 
  dhInvFeq[7] = dh[7] /fEq[7]; 
  dhInvFeq[8] = dh[8] /fEq[8]; 
  dhInvFeq[9] = dh[9] /fEq[9]; 
  dhInvFeq[10] = dh[10] /fEq[10]; 
  dhInvFeq[11] = dh[11] /fEq[11]; 
  dhInvFeq[12] = dh[12] /fEq[12]; 
  dhInvFeq[13] = dh[13] /fEq[13]; 
  dhInvFeq[14] = dh[14] /fEq[14]; 
  dhInvFeq[15] = dh[15] /fEq[15]; 
  dhInvFeq[16] = dh[16] /fEq[16]; 
  dhInvFeq[17] = dh[17] /fEq[17]; 
  dhInvFeq[18] = dh[18] /fEq[18]; 
  dhInvFeq[19] = dh[19] /fEq[19]; 
  dhInvFeq[20] = dh[20] /fEq[20]; 
  dhInvFeq[21] = dh[21] /fEq[21]; 
  dhInvFeq[22] = dh[22] /fEq[22]; 
  dhInvFeq[23] = dh[23] /fEq[23]; 
  dhInvFeq[24] = dh[24] /fEq[24]; 
  dhInvFeq[25] = dh[25] /fEq[25]; 
  dhInvFeq[26] = dh[26] /fEq[26]; 

  T entProd1 = ds1 * dhInvFeq[1] + ds2 * dhInvFeq[2] +
               ds3 * dhInvFeq[3] + ds4 * dhInvFeq[4] + ds5 * dhInvFeq[5] +
               ds6 * dhInvFeq[6] + ds7 * dhInvFeq[7] + ds8 * dhInvFeq[8] +
               ds9 * dhInvFeq[9] +
               ds14 * dhInvFeq[14] +
               ds15 * dhInvFeq[15] + ds16 * dhInvFeq[16] + ds17 * dhInvFeq[17] +
               ds18 * dhInvFeq[18] + ds19 * dhInvFeq[19] + ds20 * dhInvFeq[20] +
               ds21 * dhInvFeq[21] + ds22 * dhInvFeq[22];

  T entProd2 = dh[0] * dhInvFeq[0] + dh[1] * dhInvFeq[1] + dh[2] * dhInvFeq[2] +
               dh[3] * dhInvFeq[3] + dh[4] * dhInvFeq[4] + dh[5] * dhInvFeq[5] +
               dh[6] * dhInvFeq[6] + dh[7] * dhInvFeq[7] + dh[8] * dhInvFeq[8] +
               dh[9] * dhInvFeq[9] + dh[10] * dhInvFeq[10] + dh[11] * dhInvFeq[11] +
               dh[12] * dhInvFeq[12] + dh[13] * dhInvFeq[13] + dh[14] * dhInvFeq[14] +
               dh[15] * dhInvFeq[15] + dh[16] * dhInvFeq[16] + dh[17] * dhInvFeq[17] +
               dh[18] * dhInvFeq[18] + dh[19] * dhInvFeq[19] + dh[20] * dhInvFeq[20] +
               dh[21] * dhInvFeq[21] + dh[22] * dhInvFeq[22] + dh[23] * dhInvFeq[23] +
               dh[24] * dhInvFeq[24] + dh[25] * dhInvFeq[25] + dh[26] * dhInvFeq[26];

  T invBeta = 1. / beta;
  
  T gamma = invBeta - (2. - invBeta) * (entProd1 / entProd2);
  //std::cout << "gamma" << gamma << std::endl;

  //Collide
  cell[0] -= beta * (gamma * dh[0]);
  cell[1] -= beta * (2. * ds1 + gamma * dh[1]);
  cell[2] -= beta * (2. * ds2 + gamma * dh[2]);
  cell[3] -= beta * (2. * ds3 + gamma * dh[3]);
  cell[4] -= beta * (2. * ds4 + gamma * dh[4]);
  cell[5] -= beta * (2. * ds5 + gamma * dh[5]);
  cell[6] -= beta * (2. * ds6 + gamma * dh[6]);
  cell[7] -= beta * (2. * ds7 + gamma * dh[7]);
  cell[8] -= beta * (2. * ds8 + gamma * dh[8]);
  cell[9] -= beta * (2. * ds9 + gamma * dh[9]);
  cell[10] -= beta * (gamma * dh[10]);
  cell[11] -= beta * (gamma * dh[11]);
  cell[12] -= beta * (gamma * dh[12]);
  cell[13] -= beta * (gamma * dh[13]);
  cell[14] -= beta * (2. * ds14 + gamma * dh[14]);
  cell[15] -= beta * (2. * ds15 + gamma * dh[15]);
  cell[16] -= beta * (2. * ds16 + gamma * dh[16]);
  cell[17] -= beta * (2. * ds17 + gamma * dh[17]);
  cell[18] -= beta * (2. * ds18 + gamma * dh[18]);
  cell[19] -= beta * (2. * ds19 + gamma * dh[19]);
  cell[20] -= beta * (2. * ds20 + gamma * dh[20]);
  cell[21] -= beta * (2. * ds21 + gamma * dh[21]);
  cell[22] -= beta * (2. * ds22 + gamma * dh[22]);
  cell[23] -= beta * (gamma * dh[23]);
  cell[24] -= beta * (gamma * dh[24]);
  cell[25] -= beta * (gamma * dh[25]);
  cell[26] -= beta * (gamma * dh[26]);

  T uSqr = util::normSqr<T,3>(u);
  return uSqr;
  }

}; //struct


} //namespace olb

#endif