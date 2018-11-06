/*

 * FFequilibrium.h
 *
 *  Created on: 31/07/2018
 *      Author: Carlos Trujillo
 *
 *This file is part of the "Free Fluids" application
 *Copyright (C) 2008-2018  Carlos Trujillo Gonzalez

 *This program is free software; you can redistribute it and/or
 *modify it under the terms of the GNU General Public License version 3
 *as published by the Free Software Foundation

 *This program is distributed in the hope that it will be useful,
 *but WITHOUT ANY WARRANTY; without even the implied warranty of
 *MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *GNU General Public License for more details.

 *You should have received a copy of the GNU General Public License
 *along with this program; if not, write to the Free Software
 *Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
#ifndef FFEQUILIBRIUM_H
#define FFEQUILIBRIUM_H

#if (defined(_WIN32) || defined(__WIN32__))
  #define CALLCONV __cdecl//Definition of the calling convention, option is __stdcall
  /*** You should define EXPORTS at the gcc command line only when building the shared library ***/
  #ifdef EXPORTS
    #define EXP_IMP __declspec(dllexport)/*used to declare functions as export*/
  #else
    #ifdef IMPORTS
      #define EXP_IMP __declspec(dllimport)/*used to declare functions as import*/
    #else
      #define EXP_IMP
    #endif
  #endif
#else/*If we are not in Windows, we give EXP_IMP and CALLCONV null value*/
  #define EXP_IMP
  #define CALLCONV
#endif

#include <stdbool.h>
#include "FFbasic.h"
#include "FFeosPure.h"
#include "FFactivity.h"
#include "FFeosMix.h"

#ifdef __cplusplus
extern "C"
{
#endif
//Calculation of the minimum tangent plane distance and the corresponding composition
EXP_IMP void CALLCONV FF_StabilityCheck(FF_PTXfeed *data,int *useOptimizer,double *tpd,double tpdX[]);

//Mixture bubble temperature calculation, given P, comoposition, eos and mixing rule
EXP_IMP void CALLCONV FF_BubbleT(FF_MixData *mix,const double *P,const double x[],const double *bTguess, double *bT,double y[],
                                 double substPhiL[],double substPhiG[]);
//Mixture dew temperature calculation, given P, comoposition, eos and mixing rule
EXP_IMP void CALLCONV FF_DewT(FF_MixData *mix,const double *P, const double y[],const double *dTguess, double *dT,double x[],
                              double substPhiL[],double substPhiG[]);
//Temperature envelope of a binary mixture
void CALLCONV FF_TemperatureEnvelope(FF_MixData *mix,const double *P, const int *nPoints, double c[],double bT[],double y[],double dT[],double x[]);
//Mixture bubble pressure calculation, given T, comoposition, eos and mixing rule
EXP_IMP void CALLCONV FF_BubbleP(FF_MixData *mix,const double *T, const double c[],const double *bPguess, double *bP,
                                 double y[],double substPhiL[],double substPhiG[]);
//Mixture dew pressure calculation, given T, comoposition, eos and mixing rule
EXP_IMP void CALLCONV FF_DewP(FF_MixData *mix,const double *T, const double y[],const double *dPguess, double *dP,double x[],double substPhiL[],double substPhiG[]);
//Pressure envelope of a binary mixture
EXP_IMP void CALLCONV FF_PressureEnvelope(FF_MixData *mix,const double *T, const int *nPoints, double c[], double bP[],double y[],double dP[],double x[]);

//VL flash calculation, given T, P, feed composition, eos and mixing rule
EXP_IMP void CALLCONV FF_VLflashPT(FF_MixData *mix,const double *T,const double *P,const double f[],
                                   double x[],double y[],double substPhiL[],double substPhiG[],double *beta);
//Determines the Gibbs energy of 1 mol of a mix, given T,P,total composition and the part of it assigned to one phase
EXP_IMP double CALLCONV FF_TwoPhasesGibbs(unsigned nVar, const double coef[], double grad[], FF_PTXfeed *data);
//VL flash calculation, given T,P, composition, and thermo model to use. By global optimization of residual Gibbs energy
EXP_IMP void CALLCONV FF_TwoPhasesFlashPTGO(FF_PTXfeed *data, double x[],double y[],double substPhiL[],double substPhiG[],double *beta);
//void CALLCONV FF_VLflashPTGO();
//Mixture VL flash, given P,T, composition, and thermo model to use. By global optimization simulated annealing of residual Gibbs energy
void CALLCONV FF_TwoPhasesFlashPTSA(FF_PTXfeed *data, double x[],double y[],double substPhiL[],double substPhiG[],double *beta,double *Gr);
#ifdef __cplusplus
}
#endif
#endif // FFEQUILIBRIUM_H
